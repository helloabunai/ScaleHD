#!/bin/bash/python
__version__ = 0.321
__author__ = 'alastair.maxwell@glasgow.ac.uk'

## imports
import os

class genHTML:
    def __init__(self, shdVersion=None, jobLabel=None, sampleData=None):

        """
        docstring todo
        """

        ## Globals?
        self.TEMPLATES_BASE = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'templates')
        self.DEV_OUTPUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'OUTPUT.html') ##todo change to specified output path from ScaleHD
        self.SAMPLES = ['Sample1', 'Sample2', 'Sample3', 'Sample4', 'Sample5', 'Sample6', 'Sample7', 'Sample8', 'Sample9', 'Sample10'] ##todo change to self.instance_objects from ScaleHD

        ## Base HTML template for SHD output
        base_template = os.path.join(self.TEMPLATES_BASE, 'base.html')
        version_str = shdVersion; instancelabel_str = jobLabel # Get ScaleHD version and jobLabel
        lists_str = self.get_lists_html() # Get sample list of which SHD processed
        analysis_str = self.get_seqdata() # Get SubStage information for each sample

        ## WRITE EVERYTHING COLLECTED TO BASE HTML TEMPLATE
        ## THIS IS THE FINAL OUTPUT CURATION STAGE
        base_template = os.path.join(self.TEMPLATES_BASE, 'base.html')
        f = open(base_template, 'r')
        output = ''

        for line in f:
            line = line.format(SAMPLE_LIST=lists_str, shd_version=version_str, instance_label=instancelabel_str, SEQDATA=analysis_str)
            output = '{0}{1}'.format(output, line)
        f.close()

        ## Write to final output HMTL
        with open(self.DEV_OUTPUT,'w') as outfi:
            outfi.write(output)

    def get_lists_html(self):

        """
        docstring todo
        """

        list_template = os.path.join(self.TEMPLATES_BASE, 'list.html')
        return_str = ''

        f = open(list_template, 'r')

        for sequence in self.SAMPLES:
            for line in f:
                line = line.format(ID=sequence)
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

            ## replace with results from ScaleHD for the current sample
            test_seqqc = "SEQ_QC_TEST"
            test_seqaln = "SEQ_ALN_TEST"
            test_gtype = "GTYPE_TEST"

            for line in f:
                line = line.format(ID=sequence, SEQ_QC=test_seqqc, SEQ_ALN=test_seqaln, GTYPE=test_gtype)
                return_str = '{0}{1}'.format(return_str, line)

            f.seek(0)

        f.close()

        return return_str
