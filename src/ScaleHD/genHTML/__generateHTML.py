#!/bin/bash/python
__version__ = 0.321
__author__ = 'alastair.maxwell@glasgow.ac.uk'

## imports
import os
import pysam
import numpy as np
import pkg_resources
from .. import align
from .. import predict
from shutil import move
from fadapa import Fadapa
from shutil import copyfile
from tempfile import mkstemp
from os import fdopen, remove
from ..__backend import mkdir_p


class genHTML:
    def __init__(self, scalehdResults=None, shdVersion=None, jobLabel=None, outputPath=None):

        """
        docstring todo
        """

        ## class attributes
        self.instance_objects = scalehdResults
        self.WEB_BASE = os.path.dirname(os.path.abspath(__file__))
        self.TEMPLATES_BASE = os.path.join(self.WEB_BASE, 'templates')
        self.OUTPUT_ROOT = outputPath
        self.HTML_FILE = '{}/{}{}'.format(outputPath, jobLabel, 'HTMLResults.html')

        ## Collate ScaleHD results in meaningful way for HTML
        self.SAMPLES = []; self.get_processed_samples()

        ## Base HTML template for SHD output
        ## items to be filled with generated HTML strings from templates
        ## before being placed into their respective HTML locations
        base_template = os.path.join(self.TEMPLATES_BASE, 'base.html')
        cag_summary = self.get_summary('CAG'); ccg_summary = self.get_summary('CCG')
        allele_dict = self.implement_summaries(jobLabel, cag_summary, ccg_summary) #instance allele summary
        styling_str = self.get_styling() #CSS
        script_str = self.get_javascript() #Javascript
        version_str = shdVersion # Get ScaleHD version
        instancelabel_str = jobLabel # Get ScaleHD jobLabel
        lists_str = self.get_lists_html() # Get sample list of which SHD processed
        alleletable_str = self.get_alleletable() # Get summary table for all alleles (landing page)
        analysis_str = self.get_seqdata() # Get SubStage information for each sample

        ## WRITE EVERYTHING COLLECTED TO BASE HTML TEMPLATE
        ## THIS IS THE FINAL OUTPUT CURATION STAGE
        base_template = os.path.join(self.TEMPLATES_BASE, 'base.html')
        f = open(base_template, 'r')
        output = ''

        for line in f:
            line = line.format(
            CSS=styling_str,
            SAMPLE_LIST=lists_str, shd_version=version_str,
            instance_label=instancelabel_str,
            CAG_TITLE = allele_dict['CAG_TITLE'], CAG_DESCR=allele_dict['CAG_DESCR'], CAG_LABELS = allele_dict['CAG_LABELS'], CAG_VALUES = allele_dict['CAG_VALUES'], CAG_X = allele_dict['CAG_X'], CAG_Y = allele_dict['CAG_Y'],
            CCG_TITLE = allele_dict['CCG_TITLE'], CCG_DESCR=allele_dict['CCG_DESCR'], CCG_LABELS = allele_dict['CCG_LABELS'], CCG_VALUES = allele_dict['CCG_VALUES'], CCG_X = allele_dict['CCG_X'], CCG_Y = allele_dict['CCG_Y'],
            ALLELETABLE = alleletable_str, SEQDATA=analysis_str, JAVASCRIPT=script_str
            )
            output = '{0}{1}'.format(output, line)
        f.close()

        ## Write to final output HMTL
        with open(self.HTML_FILE,'w') as outfi:
            outfi.write(output)

    def get_processed_samples(self):
        for individual in self.instance_objects:
            self.SAMPLES.append(individual.get_label())

    def get_summary(self, triplet):

        ## CAG
        cag_blank = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
        if triplet == 'CAG':
            for individual in self.instance_objects:
                try:
                    primary_cag = individual.get_primaryallele().get_cag()
                    secondary_cag = individual.get_secondaryallele().get_cag()
                    for allele in [primary_cag, secondary_cag]:
                        index = 0
                        while index < len(cag_blank):
                            if index == allele:
                                cag_blank[index] += 1
                            index += 1
                except AttributeError: ##skip over samples that have no genotype (i.e. failed)
                    pass

            return cag_blank

        ## CCG
        ccg_blank = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
        if triplet == 'CCG':
            for individual in self.instance_objects:
                try:
                    primary_ccg = individual.get_primaryallele().get_ccg()
                    secondary_ccg = individual.get_secondaryallele().get_ccg()
                    for allele in [primary_ccg, secondary_ccg]:
                        index = 0
                        while index < len(ccg_blank):
                            if index == allele:
                                ccg_blank[index] += 1
                            index += 1
                except AttributeError:
                    pass

            return ccg_blank

    def implement_summaries(self, jobLabel, cag_summary, ccg_summary):

        ## Target
        allele_dict = {}
        cag_labels = ['0', '1', '2', '3', '4', '5', '6', '7', '8', '9', '10', '11', '12', '13', '14', '15', '16', '17', '18', '19', '20', '21', '22', '23', '24', '25', '26', '27', '28', '29', '30', '31', '32', '33', '34', '35', '36', '37', '38', '39', '40', '41', '42', '43', '44', '45', '46', '47', '48', '49', '50', '51', '52', '53', '54', '55', '56', '57', '58', '59', '60', '61', '62', '63', '64', '65', '66', '67', '68', '69', '70', '71', '72', '73', '74', '75', '76', '77', '78', '79', '80', '81', '82', '83', '84', '85', '86', '87', '88', '89', '90', '91', '92', '93', '94', '95', '96', '97', '98', '99', '100', '101', '102', '103', '104', '105', '106', '107', '108', '109', '110', '111', '112', '113', '114', '115', '116', '117', '118', '119', '120', '121', '122', '123', '124', '125', '126', '127', '128', '129', '130', '131', '132', '133', '134', '135', '136', '137', '138', '139', '140', '141', '142', '143', '144', '145', '146', '147', '148', '149', '150', '151', '152', '153', '154', '155', '156', '157', '158', '159', '160', '161', '162', '163', '164', '165', '166', '167', '168', '169', '170', '171', '172', '173', '174', '175', '176', '177', '178', '179', '180', '181', '182', '183', '184', '185', '186', '187', '188', '189', '190', '191', '192', '193', '194', '195', '196', '197', '198', '199', '200']
        ccg_labels = ['0', '1', '2', '3', '4', '5', '6', '7', '8', '9', '10', '11', '12', '13', '14', '15', '16', '17', '18', '19', '20']

        ## CAG Summary
        allele_dict['CAG_TITLE'] = 'CAG allele distribution for {}'.format(jobLabel)
        allele_dict['CAG_DESCR'] = '# of alleles present'
        allele_dict['CAG_LABELS'] = str(cag_labels)
        allele_dict['CAG_VALUES'] = str(cag_summary)
        allele_dict['CAG_X'] = 'CAG Repeat size'
        allele_dict['CAG_Y'] = 'Allele count'

        ## CCG Summary
        allele_dict['CCG_TITLE'] = 'CCG allele distribution for {}'.format(jobLabel)
        allele_dict['CCG_DESCR'] = '# of alleles present'
        allele_dict['CCG_LABELS'] = str(ccg_labels)
        allele_dict['CCG_VALUES'] = str(ccg_summary)
        allele_dict['CCG_X'] = 'CCG Repeat size'
        allele_dict['CCG_Y'] = 'Allele count'

        return allele_dict

    def get_styling(self):
        gridism_path = os.path.join(self.WEB_BASE, 'gridism.css')
        scalehd_path = os.path.join(self.WEB_BASE, 'scalehd.css')
        css_string = ''

        ## gridism styling
        f = open(gridism_path, 'r')
        for line in f:
            css_string += line
        f.close()
        ## scalehd styling
        f = open(scalehd_path, 'r')
        for line in f:
            css_string += line
        f.close()

        css_string = '<style type=\"text/css\">' + css_string + '</style>'
        return css_string

    def get_javascript(self):
        jquery_path = os.path.join(self.WEB_BASE, 'jquery.js')
        scalehd_path = os.path.join(self.WEB_BASE, 'scalehd.js')
        seqview_path = os.path.join(self.WEB_BASE, 'msa.js')
        chart_path = os.path.join(self.WEB_BASE, 'chartBase.js')
        chartBox_path = os.path.join(self.WEB_BASE, 'chartBoxPlot.js')
        chartZoom_path = os.path.join(self.WEB_BASE, 'chartZoom.js')
        js_string = ''
        
        ## jquery scripts
        f = open(jquery_path, 'r')
        for line in f:
            js_string += line
        f.close()
        ## scalehd scripts
        f = open(scalehd_path, 'r')
        for line in f:
            js_string += line
        f.close()
        ## MSA sequence viewer scripts
        f = open(seqview_path, 'r')
        for line in f:
            js_string += line
        f.close()
        ## chart.js scripts
        f = open(chart_path, 'r')
        for line in f:
            js_string += line
        f.close()
        ## chart.js boxplot extension
        f = open(chartBox_path, 'r')
        for line in f:
            js_string += line
        f.close()
        ## chart.js zoom extension
        f = open(chartZoom_path, 'r')
        for line in f:
            js_string += line
        f.close()

        js_string = '<script type=\"text/javascript\">' + js_string + '</script>'
        return js_string

    def get_lists_html(self):

        """
        docstring todo
        """

        list_template = os.path.join(self.TEMPLATES_BASE, 'list.html')
        return_str = ''

        f = open(list_template, 'r')

        for sequence in self.SAMPLES:

            ## Get SHD object for this sample
            targetObject = None
            for x in self.instance_objects:
                if x.get_label() == sequence:
                    targetObject = x

            ## Check whether an exception was raised in SHD
            unicode_indicator = '&#9658'
            if targetObject.get_exceptionraised() != 'N/A':
                unicode_indicator = '&#10060;'

            for line in f:
                line = line.format(UNICODE_INDICATOR=unicode_indicator, ID=sequence)
                return_str = '{0}{1}'.format(return_str, line)

            f.seek(0)

        f.close()
        return return_str

    def get_alleletable(self):

        tablerow_template = os.path.join(self.TEMPLATES_BASE, 'alleletable.html')
        return_str = ''

        f = open(tablerow_template, 'r')
        for sequence in self.SAMPLES:

            targetObject = None
            ## Select object from instance results
            for x in self.instance_objects:
                if x.get_label() == sequence:
                    targetObject = x

            sample_id = ''; primary_genotype = ''; primary_structure = ''; primary_confidence = ''
            secondary_genotype = ''; secondary_structure = ''; secondary_confidence = ''
            exceptions = targetObject.get_exceptionraised()
            if exceptions != 'N/A':
                primary_genotype = 'Fail'; primary_structure = 'Fail'; primary_confidence = 'Fail'
                secondary_genotype = 'Fail'; secondary_structure = 'Fail'; secondary_confidence = 'Fail'

            try:
                pri = targetObject.get_primaryallele(); sec = targetObject.get_secondaryallele()
                primary_genotype = '{}-{}'.format(pri.get_cag(), pri.get_ccg()); primary_structure = pri.get_intervening(); primary_confidence = pri.get_alleleconfidence()
                secondary_genotype = '{}-{}'.format(sec.get_cag(), sec.get_ccg()); secondary_structure = sec.get_intervening(); secondary_confidence = sec.get_alleleconfidence()
            except AttributeError:
                pass ## skip unprocessed samples / shd fail

            for line in f:
                line = line.format(
                SAMPLE_ID = sequence, A1_GENOTYPE = primary_genotype, A1_STRUCTURE = primary_structure, A1_CONFIDENCE = primary_confidence,
                A2_GENOTYPE = secondary_genotype, A2_STRUCTURE = secondary_structure, A2_CONFIDENCE = secondary_confidence
                )
                return_str = '{0}{1}'.format(return_str, line)
            f.seek(0)
        f.close()

        return return_str

    def get_seqdata(self):

        """
        docstring todo
        """

        seqdata_template = os.path.join(self.TEMPLATES_BASE, 'sequencedata.html')
        return_str = ''

        f = open(seqdata_template, 'r')
        for sequence in self.SAMPLES:

            targetObject = None
            ## Select object from instance results
            for x in self.instance_objects:
                if x.get_label() == sequence:
                    targetObject = x

            ## HEADLINE INFORMATION
            primary_fwmap = ''; primary_fwmap_pcnt = ''; primary_fwmap_purge = ''; primary_rvmap = ''; primary_rvmap_pcnt = ''; primary_rvmap_purge = ''
            secondary_fwmap = ''; secondary_fwmap_pcnt = ''; secondary_fwmap_purge = ''; secondary_rvmap = ''; secondary_rvmap_pcnt = ''; secondary_rvmap_purge = ''
            exceptions = targetObject.get_exceptionraised(); passfail=''
            if exceptions == 'N/A': passfail = 'completed'
            else: passfail = 'incomplete (exception raised on {})'.format(exceptions)

            try:
                primary_allele = targetObject.get_primaryallele()
                primary_fwmap = primary_allele.get_fwalncount(); primary_fwmap_pcnt = primary_allele.get_fwalnpcnt(); primary_fwmap_purge = primary_allele.get_fwalnrmvd()
                primary_rvmap = primary_allele.get_rvalncount(); primary_rvmap_pcnt = primary_allele.get_rvalnpcnt(); primary_rvmap_purge = primary_allele.get_rvalnrmvd()

                secondary_allele = targetObject.get_secondaryallele()
                secondary_fwmap = secondary_allele.get_fwalncount(); secondary_fwmap_pcnt = secondary_allele.get_fwalnpcnt(); secondary_fwmap_purge = secondary_allele.get_fwalnrmvd()
                secondary_rvmap = secondary_allele.get_rvalncount(); secondary_rvmap_pcnt = secondary_allele.get_rvalnpcnt(); secondary_rvmap_purge = secondary_allele.get_rvalnrmvd()

            except AttributeError:
                pass ## skip over samples that failed genotyping

            ## replace with results from ScaleHD for the current sample
            sample_seqqc = self.get_sampleQC(sequence)
            sample_seqaln = self.get_sampleALN(sequence)
            sample_gtype = self.get_sampleGTYPE(sequence)

            for line in f:
                line = line.format(
                ID=sequence, PASSFAIL=passfail,
                A1_FWMAP=primary_fwmap, A1_FWMAP_PCNT=primary_fwmap_pcnt, A1_FWMAP_PURGE=primary_fwmap_purge, A1_RVMAP=primary_rvmap, A1_RVMAP_PCNT=primary_rvmap_pcnt, A1_RVMAP_PURGE=primary_rvmap_purge,
                A2_FWMAP=secondary_fwmap, A2_FWMAP_PCNT=secondary_fwmap_pcnt, A2_FWMAP_PURGE=secondary_fwmap_purge, A2_RVMAP=secondary_rvmap, A2_RVMAP_PCNT=secondary_rvmap_pcnt, A2_RVMAP_PURGE=secondary_rvmap_purge,
                SEQ_QC=sample_seqqc, SEQ_ALN=sample_seqaln, GTYPE=sample_gtype)
                return_str = '{0}{1}'.format(return_str, line)
            f.seek(0)
        f.close()

        return return_str

    def get_sampleQC(self, currSample):

        targetObject = None
        ## Select object from instance results
        for x in self.instance_objects:
            if x.get_label() == currSample:
                targetObject = x

        ##
        ## If the sample failed, there won't be any data to collect
        ## so just return a simple string to be placed in the data's stead
        if targetObject.get_exceptionraised() == 'SeqQC':
            return '<p> No Quality Control results! ScaleHD workflow failed!</p>'

        ##################################################################
        ## Check for Trimming report! scrape data and format if present ##
        ##################################################################
        forwardTrimReport = targetObject.get_trimreport()[0]; forwardTrimString = ''
        reverseTrimReport = targetObject.get_trimreport()[1]; reverseTrimString = ''
        try:
            with open(forwardTrimReport, 'r') as infi:
                forwardTrimString=infi.read()
            forwardTrimString = self.format_trimming(forwardTrimString)
        except Exception, e:
            forwardTrimString = 'We could not find/process a forward trimming report! Exception raised: {}'.format(e)
        try:
            with open(reverseTrimReport, 'r') as infi:
                reverseTrimString=infi.read()
            reverseTrimString = self.format_trimming(reverseTrimString)
        except Exception, e:
            reverseTrimString = 'We could not find/process a reverse trimming report! Exception raised: {}'.format(e)

        ################################################################
        ## Check for FastQC report! scrape data and format if present ##
        ################################################################
        forwardFQCReport = targetObject.get_fqcreport()[0]; forwardFQCString = ''
        try:
            forwardFQCString = self.format_fastqc(forwardFQCReport, currSample)
        except Exception, e:
            forwardFQCString = 'We could not find/process a FastQC report! Exception raised: {}'.format(e)

        #########################################################################
        ## Get FastQC graph data and append it to attribute tags in seqqc.html ##
        #########################################################################
        forwardFQCReport = targetObject.get_fqcreport()[0]; fastqc_graphdata = {}
        try:
            fastqc_graphdata = self.format_fastqc_graphs(forwardFQCReport, currSample)
        except Exception, e:
            fastqc_graphdata = {
            'PBSQ_TITLE': 'FastQC failure :(', 'PBSQ_LABELS': '', 'PBSQ_VALUES': '', 'PBSQ_MEANVAL': '', 'PBSQ_DESCR': '', 'PBSQ_X': '', 'PBSQ_Y': '',
            'PBNC_TITLE': 'FastQC failure :(', 'PBNC_LABELS': '', 'PBNC_VALUES': '', 'PBNC_DESCR': '', 'PBNC_X': '', 'PBNC_Y': '',
            'SQLD_TITLE': 'FastQC failure :(', 'SQLD_LABELS': '', 'SQLD_VALUES': '', 'SQLD_DESCR': '', 'SQLD_X': '', 'SQLD_Y': ''
            }

        ###################################################################
        ## Apply scraped and formatted data into HTML template for SeqQC ##
        ###################################################################
        qc_template = os.path.join(self.TEMPLATES_BASE, 'seqqc.html')
        f = open(qc_template, 'r')
        qc_return = ''
        for line in f:
            line = line.format(ID=currSample, FORWARD_TRIM=forwardTrimString, REVERSE_TRIM=reverseTrimString, FASTQC=forwardFQCString,
            PBSQ_TITLE=fastqc_graphdata['PBSQ_TITLE'], PBSQ_LABELS=fastqc_graphdata['PBSQ_LABELS'], PBSQ_VALUES=fastqc_graphdata['PBSQ_VALUES'], PBSQ_MEANVAL=fastqc_graphdata['PBSQ_MEANVAL'], PBSQ_DESCR=fastqc_graphdata['PBSQ_DESCR'], PBSQ_X=fastqc_graphdata['PBSQ_X'], PBSQ_Y=fastqc_graphdata['PBSQ_Y'],
            PBNC_TITLE=fastqc_graphdata['PBNC_TITLE'], PBNC_LABELS=fastqc_graphdata['PBNC_LABELS'], PBNC_VALUES=fastqc_graphdata['PBNC_VALUES'], PBNC_DESCR=fastqc_graphdata['PBNC_DESCR'], PBNC_X=fastqc_graphdata['PBNC_X'], PBNC_Y=fastqc_graphdata['PBNC_Y'],
            SQLD_TITLE=fastqc_graphdata['SQLD_TITLE'], SQLD_LABELS=fastqc_graphdata['SQLD_LABELS'], SQLD_VALUES=fastqc_graphdata['SQLD_VALUES'], SQLD_DESCR=fastqc_graphdata['SQLD_DESCR'], SQLD_X=fastqc_graphdata['SQLD_X'], SQLD_Y=fastqc_graphdata['SQLD_Y'],
            )
            qc_return = '{0}{1}'.format(qc_return, line)
        f.close()

        return qc_return

    def format_trimming(self, rawData):

        ## Trimming template
        trim_template = os.path.join(self.TEMPLATES_BASE, 'trim.html')

        ## First get technical summary; replace python newline with HTML break, remove start/end breaks
        techSummary = rawData.split('=== Summary ===')[0]
        techSummary = techSummary.replace('\n', '<br />').lstrip('<br />').rstrip('<br />')

        ## Trimming Summary; replace python newline with HTML break, remove start/end breaks
        trimSummary = rawData.split('=== Summary ===')[1].split('=== Adapter 1 ===')[0]
        trimSummary = trimSummary.replace('\n', '<br />').lstrip('<br />').rstrip('<br />')
        trimSummary = trimSummary.replace('<br /><br />', '<br />')

        ## Adapter summary; replace python newline with HTML break, remove start/end breaks
        adapterSummary = rawData.split('=== Adapter 1 ===')[1].split('Overview of removed sequences')[0]
        adapterSummary = adapterSummary.replace('\n', '<br />').lstrip('<br />').rstrip('<br />')
        adapterSummary = adapterSummary.replace(';', '<br />')

        trim_return = ''
        f = open(trim_template, 'r')
        for line in f:
            line = line.format(TECHNICAL=techSummary, SUMMARY=trimSummary, ADAPTER=adapterSummary)
            trim_return = '{0}{1}'.format(trim_return, line)
        f.close()

        ## Return formatted trimming report
        return trim_return

    def format_fastqc(self, rawDataPath, currSample):

        ## FastQC templates
        fastqc_template = os.path.join(self.TEMPLATES_BASE, 'fastqc.html')

        ## just fuckin lump it all in there for now and figure out what you want to format next
        fqc_object = Fadapa(rawDataPath)

        ## Module status data
        module_summary = fqc_object.summary()
        module_stats = module_summary[1][0]; module_pbsq = module_summary[2][0]; module_ptsq = module_summary[3][0];
        module_psqs = module_summary[4][0]; module_pbsc = module_summary[5][0]; module_psgcc = module_summary[6][0];
        module_pbnc = module_summary[7][0]; module_seqlendist = module_summary[8][0]; module_seqdup = module_summary[9][0];
        module_overrep = module_summary[10][0]; module_adapter = module_summary[11][0]

        ## Basic statistics data
        basic_stats = fqc_object.clean_data('Basic Statistics')
        file_name = basic_stats[1][1]; file_type = basic_stats[2][1]; encoding = basic_stats[3][1]
        total_sequences = basic_stats[4][1]; poor_quality = basic_stats[5][1]; seq_len = basic_stats[6][1]
        gc_pcnt = basic_stats[7][1]

        ## FastQC html template file with data inserted
        fqc_return = ''
        f = open(fastqc_template, 'r')
        for line in f:
            line = line.format(
            MODULE_STATS = module_stats, MODULE_PBSQ = module_pbsq, MODULE_PTSQ = module_ptsq,
            MODULE_PSQS = module_psqs, MODULE_PBSC = module_pbsc, MODULE_PSGCC = module_psgcc,
            MODULE_PBNC = module_pbnc, MODULE_SEQLENDIST = module_seqlendist, MODULE_SEQDUP = module_seqdup,
            MODULE_OVERREP = module_overrep, MODULE_ADAPTER = module_adapter,
            FQC_FILENAME = file_name, FQC_FILETYPE = file_type, FQC_ENCODING = encoding,
            FQC_TOTALSEQ = total_sequences, FQC_POORQUAL = poor_quality, FQC_SEQLEN = seq_len,
            FQC_GCPCNT = gc_pcnt
            )
            fqc_return = '{0}{1}'.format(fqc_return, line)
        f.close()

        ## return formatted FastQC report
        return fqc_return

    def format_fastqc_graphs(self, rawDataPath, currSample):

        ## Object to get data from FastQC output
        fqc_object = Fadapa(rawDataPath)

        ## Target output dictionary
        fastqc_graphdata = {}

        ## Unextracted data
        fqc_pbsq_data = fqc_object.clean_data('Per base sequence quality')
        fqc_pbnc_data = fqc_object.clean_data('Per base N content')
        fqc_seqlen_data = fqc_object.clean_data('Sequence Length Distribution')

        ##
        ## Per Base Pair Sequence Quality
        ## min = item[5], q1 = item[3], median = item[2], q3 = item[4], max = item[6]
        pbsq_labels = []; pbsq_values = []; pbsq_means = []
        for item in fqc_pbsq_data[1:]:
            pbsq_labels.append(item[0]) ## label for bin
            pbsq_means.append(int(float(item[1]))) ## sample running mean
            bin_values = [item[5],item[3],item[2],item[4],item[6]]
            bin_values = [0.0 if x=='NaN' else x for x in bin_values] ## replace NaN with 0
            bin_values = [int(float(x)) for x in bin_values] ## convert str of float->float->int
            pbsq_values.append(bin_values)
        fastqc_graphdata['PBSQ_TITLE'] = 'FastQC Per base sequence quality'
        fastqc_graphdata['PBSQ_LABELS'] = str(pbsq_labels)
        fastqc_graphdata['PBSQ_VALUES'] = str(pbsq_values)
        fastqc_graphdata['PBSQ_MEANVAL'] = str(pbsq_means)
        fastqc_graphdata['PBSQ_DESCR'] = 'Per base sequence quality'
        fastqc_graphdata['PBSQ_X'] = 'Position in read (BP)'
        fastqc_graphdata['PBSQ_Y'] = 'PHRED quality score'

        ##
        ## Per Base Pair N Content
        fastqc_graphdata['PBNC_TITLE'] = 'FastQC Per base N content for {}'.format(currSample)
        pbnc_labels = []; pbnc_values = []
        for item in fqc_pbnc_data[1:]:
            pbnc_labels.append(item[0]); pbnc_values.append(item[1])
        fastqc_graphdata['PBNC_LABELS'] = str(pbnc_labels)
        fastqc_graphdata['PBNC_VALUES'] = str(pbnc_values)
        fastqc_graphdata['PBNC_DESCR'] = 'N content per base'
        fastqc_graphdata['PBNC_X'] = 'Position in read (BP)'
        fastqc_graphdata['PBNC_Y'] = 'Percentage content (%)'

        ##
        ## Sequence Length Distribution
        fastqc_graphdata['SQLD_TITLE'] = 'FastQC Sequence length distribution for {}'.format(currSample)
        dist_labels = []; dist_values = []
        for item in fqc_seqlen_data[1:]:
            dist_labels.append(item[0]); dist_values.append(item[1])
        fastqc_graphdata['SQLD_LABELS'] = str(dist_labels)
        fastqc_graphdata['SQLD_VALUES'] = str(dist_values)
        fastqc_graphdata['SQLD_DESCR'] = 'Sequence length population'
        fastqc_graphdata['SQLD_X'] = 'Sequence length (BP)'
        fastqc_graphdata['SQLD_Y'] = 'Population (#)'

        return fastqc_graphdata

    def get_sampleALN(self, currSample):

        targetObject = None
        ## Select object from instance results
        for x in self.instance_objects:
            if x.get_label() == currSample:
                targetObject = x

        ##
        ## If the sample failed, there won't be any data to collect
        ## so just return a simple string to be placed in the data's stead
        if targetObject.get_exceptionraised() in ['SeqALN','SeqRE-ALN','DSP','Genotype','SNPCalling']:
            return '<p> No Genotype results! ScaleHD workflow failed/incomplete!</p>'

        ##
        ## Primary allele alignment map
        pri_assembly_object = pysam.AlignmentFile(targetObject.get_primaryallele().get_fwassembly(), 'rb')
        pri_contig = targetObject.get_primaryallele().get_reflabel()
        pri_sequences = ''; counter = 1; pri_reads = None; pri_err_string = ''
        ## if atypical, then the labelling format generated in custom XML/FA is different
        if targetObject.get_primaryallele().get_allelestatus() == 'Atypical':
            pri_contig = '{}_CAG{}_CCG{}_CCT{}'.format(pri_contig, targetObject.get_primaryallele().get_cag(), targetObject.get_primaryallele().get_ccg(), targetObject.get_primaryallele().get_cct())

        try:
            pri_reads = pri_assembly_object.fetch(reference=pri_contig)
        except ValueError:
            ## Atypical fuckery where the specified contig does not exist
            ## find best match of available and take those reads (DSP wrong allele size?)
            original = targetObject.get_primaryallele().get_reflabel(); present_references = pri_assembly_object.references; similar_contigs = []
            for item in present_references: similar_contigs.append((item, align.similar(sec_contig,item)))
            pri_contig = sorted(similar_contigs, key=lambda a: a[1], reverse=True)[0][0]
            pri_reads = pri_assembly_object.fetch(reference=pri_contig)
            pri_err_string = '<p>ScaleHD was unable to extract reads for the contig: {}. Extracted data is from the best contig match: {}</p>'.format(targetObject.get_primaryallele().get_reflabel(), pri_contig)
        for read in pri_reads:
            target_sequence = read.query_alignment_sequence
            if counter <= 50:
                pri_sequences += ">{}\n{}\n".format(counter, target_sequence)
            counter += 1

        ##
        ## Secondary allele alignment map
        sec_assembly_object = pysam.AlignmentFile(targetObject.get_secondaryallele().get_fwassembly(), 'rb')
        sec_contig = targetObject.get_secondaryallele().get_reflabel()
        sec_sequences = ''; counter = 1; sec_reads = None; sec_err_string = ''
        ## if atypical, then the labelling format generated in custom XML/FA is different
        if targetObject.get_secondaryallele().get_allelestatus() == 'Atypical':
            sec_contig = '{}_CAG{}_CCG{}_CCT{}'.format(sec_contig, targetObject.get_secondaryallele().get_cag(), targetObject.get_secondaryallele().get_ccg(), targetObject.get_secondaryallele().get_cct())

        try:
            sec_reads = sec_assembly_object.fetch(reference=sec_contig)
        except ValueError:
            ## Atypical fuckery where the specified contig does not exist
            ## find best match of available and take those reads (DSP wrong allele size?)
            original = targetObject.get_secondaryallele().get_reflabel(); present_references = sec_assembly_object.references; similar_contigs = []
            for item in present_references: similar_contigs.append((item, align.similar(sec_contig,item)))
            sec_contig = sorted(similar_contigs, key=lambda a: a[1], reverse=True)[0][0]
            sec_reads = sec_assembly_object.fetch(reference=sec_contig)
            sec_err_string = '<p>ScaleHD was unable to extract reads for the contig: {}. Extracted data is from the best contig match: {}</p>'.format(targetObject.get_secondaryallele().get_reflabel(), sec_contig)
        for read in sec_reads:
            target_sequence = read.query_alignment_sequence
            if counter <= 50:
                sec_sequences += ">{}\n{}\n".format(counter, target_sequence)
            counter += 1

        ###################################################################
        ## Apply scraped and formatted data into HTML template for SeqQC ##
        ###################################################################
        aln_template = os.path.join(self.TEMPLATES_BASE, 'seqALN.html')
        f = open(aln_template, 'r')
        aln_return = ''
        for line in f:
            line = line.format(ID = currSample, PRI_ERR_STRING = pri_err_string, PRI_CONTIG = pri_contig, SEC_ERR_STRING = sec_err_string, SEC_CONTIG = sec_contig, PRI_SEQUENCES = pri_sequences, SEC_SEQUENCES = sec_sequences)
            aln_return = '{0}{1}'.format(aln_return, line)
        f.close()

        return aln_return

    def get_sampleGTYPE(self, currSample):

        targetObject = None
        ## Select object from instance results
        for x in self.instance_objects:
            if x.get_label() == currSample:
                targetObject = x
        gtype_data = {}

        ##
        ## If the sample failed, there won't be any data to collect
        ## so just return a simple string to be placed in the data's stead
        if targetObject.get_exceptionraised() in ['SeqALN','SeqRE-ALN','DSP','Genotype','SNPCalling']:
            return '<p> No Genotype results! ScaleHD workflow failed/incomplete!</p>'

        ##
        ## Allele objects because we'll be using them a lot durrr
        primary_allele = targetObject.get_primaryallele()
        secondary_allele = targetObject.get_secondaryallele()

        ##################################
        ## Summary genotype information ##
        ##################################
        pri_cag = primary_allele.get_cag(); pri_ccg = primary_allele.get_ccg(); pri_structurelabel = primary_allele.get_allelestatus();
        pri_structure = primary_allele.get_reflabel(); pri_intervening = primary_allele.get_intervening(); pri_slippage = primary_allele.get_backwardsslippage()
        pri_mosaicism = primary_allele.get_somaticmosaicism(); pri_confidence = primary_allele.get_alleleconfidence()

        sec_cag = secondary_allele.get_cag(); sec_ccg = secondary_allele.get_ccg(); sec_structurelabel = secondary_allele.get_allelestatus();
        sec_structure = secondary_allele.get_reflabel(); sec_intervening = secondary_allele.get_intervening(); sec_slippage = secondary_allele.get_backwardsslippage()
        sec_mosaicism = secondary_allele.get_somaticmosaicism(); sec_confidence = secondary_allele.get_alleleconfidence()

        ###############################################
        ## Data for CCG distribution for this sample ##
        ###############################################
        ccg_labels = ['0', '1', '2', '3', '4', '5', '6', '7', '8', '9', '10', '11', '12', '13', '14', '15', '16', '17', '18', '19', '20']
        pri_rvarray = primary_allele.get_rvarray(); sec_rvarray = secondary_allele.get_rvarray()
        allele_super = [0] + np.asarray([a + b for a, b in zip(pri_rvarray,sec_rvarray)]).tolist() ## 0 added to offset label indexing from 0

        gtype_data['CCGDIST_TITLE'] = 'CCG Distribution for {}'.format(currSample)
        gtype_data['CCGDIST_DESCR'] = '# of reads present'
        gtype_data['CCGDIST_LABELS'] = str(ccg_labels)
        gtype_data['CCGDIST_VALUES'] = str(allele_super)
        gtype_data['CCGDIST_X'] = 'CCG Repeat size'
        gtype_data['CCGDIST_Y'] = 'Number of reads'

        ############################################################
        ##          Homozygous / Heterozygous CAG graphs          ##
        ## I am lazy so i'm just gonna overlay the datasets hurrr ##
        ############################################################
        cag_labels = ['0', '1', '2', '3', '4', '5', '6', '7', '8', '9', '10', '11', '12', '13', '14', '15', '16', '17', '18', '19', '20', '21', '22', '23', '24', '25', '26', '27', '28', '29', '30', '31', '32', '33', '34', '35', '36', '37', '38', '39', '40', '41', '42', '43', '44', '45', '46', '47', '48', '49', '50', '51', '52', '53', '54', '55', '56', '57', '58', '59', '60', '61', '62', '63', '64', '65', '66', '67', '68', '69', '70', '71', '72', '73', '74', '75', '76', '77', '78', '79', '80', '81', '82', '83', '84', '85', '86', '87', '88', '89', '90', '91', '92', '93', '94', '95', '96', '97', '98', '99', '100']
        pri_fwarray = primary_allele.get_fwarray(); sec_fwarray = secondary_allele.get_fwarray()
        pri_split = predict.split_cag_target(pri_fwarray); sec_split = predict.split_cag_target(sec_fwarray)
        pri_target = pri_split['CCG{}'.format(pri_ccg)].tolist(); sec_target = sec_split['CCG{}'.format(sec_ccg)].tolist()
        pri_target = pri_target[0:100]; sec_target = sec_target[0:100]

        gtype_data['CAGDIST_TITLE'] = 'CAG Distribution for {}'.format(currSample)
        gtype_data['PRI_DIST_DESCR'] = 'CCG{}'.format(pri_ccg)
        gtype_data['SEC_DIST_DESCR'] = 'CCG{}'.format(sec_ccg)
        gtype_data['CAGDIST_LABELS'] = str(cag_labels)
        gtype_data['CAGDIST_PRI_VALUES'] = str(pri_target)
        gtype_data['CAGDIST_SEC_VALUES'] = str(sec_target)
        gtype_data['CAGDIST_X'] = 'CAG Repeat size'
        gtype_data['CAGDIST_Y'] = 'Number of reads'

        #################
        ## SNP Calling ##
        #################
        pri_snp = primary_allele.get_variantcall(); pri_score = primary_allele.get_variantscore()
        sec_snp = secondary_allele.get_variantcall(); sec_score = secondary_allele.get_variantscore()

        ############################
        ## ScaleHD analysis flags ##
        ############################
    	shd_exception = targetObject.get_exceptionraised(); shd_homozygous = targetObject.get_homozygoushaplotype(); shd_neighbours = targetObject.get_neighbouringpeaks();
        shd_diminished = targetObject.get_diminishedpeaks(); shd_novelatypical = targetObject.get_novel_atypical_structure(); shd_alignmentwarn = targetObject.get_alignmentwarning()
        shd_atypicalalignmentwarn = targetObject.get_atypical_alignmentwarning(); shd_ccgrewrite = targetObject.get_atypical_ccgrewrite(); shd_zygrewrite = targetObject.get_atypical_zygrewrite()
        shd_ccguncertain = targetObject.get_ccguncertainty(); shd_cctuncertain = targetObject.get_cctuncertainty(); shd_svmfail = targetObject.get_svm_failure();
        shd_diffconfuse = targetObject.get_differential_confusion(); shd_missedexpansion = targetObject.get_missed_expansion(); shd_heuristicfilter = targetObject.get_heuristicfilter();
        shd_peakinspection = targetObject.get_peakinspection_warning(); shd_lowdistreads = targetObject.get_distribution_readcount_warning(); shd_lowpeakreads = targetObject.get_fatalreadallele()

        ###################################################################
        ## Apply scraped and formatted data into HTML template for SeqQC ##
        ###################################################################
        gtype_template = os.path.join(self.TEMPLATES_BASE, 'seqGTYPE.html')
        f = open(gtype_template, 'r')
        gtype_return = ''
        for line in f:
            line = line.format(
            ID=currSample,
            A1_CAG = pri_cag, A1_CCG = pri_ccg, A1_STRUCTURELABEL = pri_structurelabel, A1_STRUCTURE = pri_structure, A1_INTERVENING = pri_intervening, A1_SLIPPAGE = pri_slippage, A1_MOSAICISM = pri_mosaicism, A1_CONFIDENCE = pri_confidence,
            A2_CAG = sec_cag, A2_CCG = sec_ccg, A2_STRUCTURELABEL = sec_structurelabel, A2_STRUCTURE = sec_structure, A2_INTERVENING = sec_intervening, A2_SLIPPAGE = sec_slippage, A2_MOSAICISM = sec_mosaicism, A2_CONFIDENCE = sec_confidence,
            CCGDIST_TITLE = gtype_data['CCGDIST_TITLE'], CCGDIST_DESCR = gtype_data['CCGDIST_DESCR'], CCGDIST_LABELS = gtype_data['CCGDIST_LABELS'],
            CCGDIST_VALUES = gtype_data['CCGDIST_VALUES'], CCGDIST_X = gtype_data['CCGDIST_X'], CCGDIST_Y = gtype_data['CCGDIST_Y'],
            CAGDIST_TITLE = gtype_data['CAGDIST_TITLE'], PRI_DIST_DESCR = gtype_data['PRI_DIST_DESCR'], SEC_DIST_DESCR=gtype_data['SEC_DIST_DESCR'],
            CAGDIST_LABELS = gtype_data['CAGDIST_LABELS'], CAGDIST_PRI_VALUES = gtype_data['CAGDIST_PRI_VALUES'], CAGDIST_SEC_VALUES = gtype_data['CAGDIST_SEC_VALUES'], CAGDIST_X = gtype_data['CAGDIST_X'], CAGDIST_Y = gtype_data['CAGDIST_Y'],
            A1_SNP = pri_snp, A1_CALLSCORE = pri_score, A2_SNP = sec_snp, A2_CALLSCORE = sec_score,
            SHDFLAG_EXCEPTION = shd_exception, SHDFLAG_HZYG = shd_homozygous, SHDFLAG_NEIGHBOUR = shd_neighbours, SHDFLAG_DIMINISH = shd_diminished, SHDFLAG_NOVELATYP = shd_novelatypical, SHDFLAG_ALNWARN = shd_alignmentwarn,
            SHDFLAG_ATYPALIGNWARN = shd_atypicalalignmentwarn, SHDFLAG_CCGREWR = shd_ccgrewrite, SHDFLAG_CCGZYG_REWR = shd_zygrewrite, SHDFLAG_CCGUNCERTAIN = shd_ccguncertain, SHDFLAG_CCTUNCERTAIN = shd_cctuncertain,
            SHDFLAG_SVMFAIL = shd_svmfail, SHDFLAG_DIFFCONFUSE = shd_diffconfuse, SHDFLAG_MISSEDEXP = shd_missedexpansion, SHDFLAG_FILTERPASS = shd_heuristicfilter, SHDFLAG_PEAKINSPECT = shd_peakinspection,
            SHDFLAG_LOWDISTREADS = shd_lowdistreads, SHDFLAG_LOWPEAKREADS = shd_lowpeakreads
            )
            gtype_return = '{0}{1}'.format(gtype_return, line)
        f.close()
        return gtype_return
