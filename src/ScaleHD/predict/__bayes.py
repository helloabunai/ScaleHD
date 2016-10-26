#/usr/bin/python
__version__ = 0.01
__author__ = 'alastair.maxwell@glasgow.ac.uk'

##
## Generic imports
import os
import pandas

class BayesianLikelihood:
    def __init__(self, bayesian_path, raw_matrix, likely_matrix):
        """
        Class which utilises bayesian statistics to determine a genotype sample
        Currently an alternative approach to SVM/DE/FOD utilised in __prediction; depending on the
        strength of the performance, it may replace/be included along-side as a QoL metric for confidence.
        Work in progres.. more details as they are determined
        :param bayesian_path: path to output folder for bayes results
        :param raw_matrix: data matrix included with the package; reads for CAG1-200/CCG1-20
        """
        self.bayesian_path = bayesian_path
        self.raw_matrix = raw_matrix

        ##
        ## Read log likelihood heatmap ..
        ## Remove duplicates from DataFrame
        self.polyglu_probabilities = pandas.read_table(likely_matrix)
        self.polyglu_matrix = self.remove_bias()

        #self.main()

    def remove_bias(self):
        """
        Ensure we only have unique samples within our input matrix
        :return: unique_polyglu (matrix of only unique samples)
        """

        ##
        ## Read the matrix as is into DataFrame and drop duplicates (would skew bayes)
        polyglu_matrix = pandas.read_csv(self.raw_matrix,header=1,sep=',')
        unique_polyglu = polyglu_matrix.drop_duplicates()
        return unique_polyglu

    def bayes_genotype(self, p_mat, read_data):

        pp = None
        for index, row in self.polyglu_probabilities.iterrows():
            nind = None ## returns CAGxCCG1-20; x+1 for x in len(0,200)
            test = []
            read_counts = None ## for current read_data, gets read count at index from nind (dict)
            ps = None
            pp = None

        return 'hi'

    def main(self):

        tags = []
        likelihoods = None
        for index, row in self.polyglu_matrix.iterrows():

            ##
            ## Get the current row
            read_data = row
            bayes_m = self.bayes_genotype(self.polyglu_probabilities,read_data)