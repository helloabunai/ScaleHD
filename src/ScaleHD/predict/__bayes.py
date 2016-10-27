#/usr/bin/python
__version__ = 0.01
__author__ = 'alastair.maxwell@glasgow.ac.uk'

##
## Generic imports
import gc
import os
import sys
import csv
import socket
import numpy as np
import pandas as pd
import logging as log

##
## Package backend
from ..__backend import Colour as clr
from . import GenotypePrediction as gtp

##
## R <-> Python interface
import rpy2.robjects as robj
from rpy2.robjects import pandas2ri
import rpy2.robjects.packages as rpackages
from rpy2.robjects.vectors import StrVector

class BayesianLikelihood:
    def __init__(self, bayesian_path, data_pair, likely_matrix, raw_matrix):
        """
        Class which utilises bayesian statistics to determine a genotype sample
        Currently an alternative approach to SVM/DE/FOD utilised in __prediction; depending on the
        strength of the performance, it may replace/be included along-side as a QoL metric for confidence.
        Work in progress.. more details as they are determined
        :param bayesian_path: path to output folder for bayes results
        :param data_pair: fw/rv data files for this current sample
        :param likely_matrix: likelihood heatmap for HD data
        :param raw_matrix: read count distributions
        """

        ##
        ## Instantiate variables
        self.bayesian_path = bayesian_path
        self.data_pair = data_pair
        self.polyglu_matrix = likely_matrix
        self.raw_matrix = raw_matrix

        ##
        ## Get input distributions for this sample
        self.forward_distribution = gtp.scrape_distro(self.data_pair[0])
        self.reverse_distribution = gtp.scrape_distro(self.data_pair[1])

        ##
        ## Prepare input distributions as dataframes for R functions
        self.forward_frame = self.create_frame(self.forward_distribution)
        self.reverse_frame = self.create_frame(self.reverse_distribution)

        ##
        ## Check we have an internet connection (for R package install/CRAN)
        self.is_connected = True
        if not self.check_connection():
            log.info('{}{}{}{}'.format(clr.yellow, 'shd__ ', clr.end, 'No internet connection!!'))
            log.info('{}{}{}{}'.format(clr.yellow, 'shd__ ', clr.end, 'Hopefully you are not missing any R packages.'))

        ##
        ## Installed required R Package (if missing)
        ## Note: R spam is annoying when loading packages etc so..
        ## re-route it to DEVNULL
        self.install_r_packages()

        ##
        ## Since we can't assume there is an internet connection
        ## and the above may have failed, check that the required R packages are installed..
        self.check_r_packages()

        ##
        ## We need to make python objects of our R-genotyping functions
        ## Do so here..
        self.get_uniques = None
        self.weight_matrix = None
        self.slippage_handler = None
        self.genotype_caller = None
        self.generate_python_r_functions()

        ##
        ## Run the genotyping..
        self.main()

    def create_frame(self, input_distro):
        """
        docstring...
        :return:
        """
        ##
        ## Remove CAG101-200 (since bayesian model lacks)
        ## Add a nolabel to the end
        ccg_split = np.split(input_distro, 20)
        new_distribution = []
        for ccg_contig in ccg_split:
            half = ccg_contig[0:100]
            new_distribution = np.hstack((new_distribution, half))
        labelled_newdistro = np.hstack((new_distribution, ["NoLabel"]))

        ##TODO >>> d = {'a': robjects.IntVector((1,2,3)), 'b': robjects.IntVector((4,5,6))}
        ##TODO >>> dataf = robject.DataFrame(d)

        ##
        ## Get a header (laziness..)
        with open(self.raw_matrix, 'r') as csvfi:
            reader = csv.reader(csvfi)
            header = next(reader)
            csvfi.close()
            header = header[1:]

        df = pd.DataFrame(labelled_newdistro, header).transpose()
        return df

    def check_connection(self):
        """
        docstr...
        :return:
        """
        try:
            host_socket = socket.gethostbyname('www.google.com') #if google is down the world has ended
            socket.create_connection((host_socket,80),2)
            return True
        except:
            self.is_connected = False
            pass
        socket.socket.close()
        return False

    @staticmethod
    def install_r_packages():
        """
        docstr..
        :return:
        """
        ##
        ## Route stdout to NULL
        stdout = sys.stdout
        f = open(os.devnull,'w')
        sys.stdout = f

        ##
        ## Activate R environment interface
        ## Install required packages if missing (assuming internet connection)
        pandas2ri.activate()
        rutils = rpackages.importr('utils')
        rutils.chooseCRANmirror(ind=1)
        packnames = ['proxy']
        install_targets = [x for x in packnames if not rpackages.isinstalled(x)]
        if len(install_targets) > 0: rutils.install_packages(StrVector(install_targets))

        ##
        ## Return stdout to stdout
        sys.stdout = stdout

    def check_r_packages(self):
        """
        docstr..
        :return:
        """
        packnames = ['proxy']
        install_targets = [x for x in packnames if not rpackages.isinstalled(x)]
        if len(install_targets) > 0 and self.is_connected == False:
            raise IOError('{}{}'.format('You are missing R packages and have no internet connection: ', install_targets))

    def generate_python_r_functions(self):
        """
        docstring...
        :return:
        """

        ##
        ## R comment
        self.weight_matrix = robj.r(
            '''
            weight_matrix <- function(pp, rcounts, p=0.1)
            {
                require("proxy")
                point = which(pp == min(pp), arr.ind = TRUE)
                edistmat =p*dist(-point[1]:((nrow(pp)-point[1])-1), -point[2]:((ncol(pp)-point[2])-1))
                e = ((edistmat/(max(edistmat)))*pp[point])*rcounts
                return(list("alleleA"=which(e == min(e), arr.ind = TRUE), "alleleB"=point))
            }
            ''')

        ##
        ## R comment
        self.slippage_handler = robj.r(
            '''
            ccg_cag_bayes_slippage <- function(p_mat, read_data, p=0.1)
            {
                pp = NULL
                rcounts = p_mat

                for(y in 1:length(p_mat[,1]))
                {
                    nind = sapply(1:ncol(p_mat), function(x) paste(row.names(p_mat)[y], colnames(p_mat)[x], sep=""))
                    read_counts = read_data[c(nind)]
                    rcounts[y,] = read_counts
                    ps = p_mat[y,]*read_counts
                    pp = rbind(pp, ps)
                }

                return((weight_matrix(pp, rcounts, p)))
            }
            ''')

        ##
        ## R comment
        self.genotype_caller = robj.r(
            '''
            run_bayes_ccg_cag_calling <- function(p_mat, data, p=0.1)
            {
                tags = vector()

                for(x in 1:length(data[,1]))
                {
                    read_data = data[x,]
                    bayes_m = ccg_cag_bayes_slippage(p_mat, read_data, p)
                    tags[x] = paste("CAG", bayes_m$alleleA[1], "CCG", bayes_m$alleleA[2], "-", "CAG", bayes_m$alleleB[1], "CCG",
                    bayes_m$alleleB[2], sep="")
                }

                return(tags)
            }
            ''')

    def main(self):
        """
        Main function of this class
        Fill out docstring later on..
        :return:
        """

        ##
        ## Set up DataFrame objects for the likelihood matrix
        robj.globalenv["p_mat"] = robj.DataFrame.from_csvfile(self.polyglu_matrix, header=True, row_names=1)
        robj.globalenv["raw_data"] = robj.DataFrame.from_csvfile(self.raw_matrix, header=True, sep=",", row_names=1)

        ##
        ## Build a DataFrame object for our input distribution for this sample...
        #calls = self.genotype_caller(robj.globalenv["p_mat"], robj.globalenv["data"])
        #print calls

        print self.forward_frame

        ##
        ## Build environment frames for our input distributions
        for distribution in [self.forward_frame, self.reverse_frame]:
            robj.globalenv["data"] = distribution
            prediction = self.genotype_caller(robj.globalenv["p_mat"], robj.globalenv["data"])
            print ">> Sample Prediction: ", prediction
