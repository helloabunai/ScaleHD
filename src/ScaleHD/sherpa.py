#/usr/bin/python
__version__ = 0.01
__author__ = 'alastair.maxwell@glasgow.ac.uk'

import os
import sys
import argparse
import pkg_resources
import logging as log
import glob
from multiprocessing import cpu_count

##
## Backend junk
from __backend import ConfigReader
from __backend import Colour as clr
from __backend import initialise_libraries
from __backend import sanitise_inputs
from __backend import extract_data
from __backend import sequence_pairings
from __backend import sanitise_outputs

##
## Package stages
from . import seq_qc
from . import align
from . import predict

##
## Globals
VERBOSE = False
LOGGING = True
THREADS = cpu_count()
DEF_OUT = os.path.join(os.path.expanduser('~'),'ScaleHD')

class ScaleHD:
	def __init__(self):

		"""
		ScaleHD -- automated Huntington Disease genotyping pipeline
		~~haha fill this out~~
		"""

		##
		## Package data
		self.generic_descriptor = pkg_resources.resource_filename(__name__, 'train/long_descr.rst')
		self.collapsed_ccg_zygosity = pkg_resources.resource_filename(__name__, 'train/polyglu_zygosity.csv')
		self.training_data = {'GenericDescriptor': self.generic_descriptor, 'CollapsedCCGZygosity': self.collapsed_ccg_zygosity}

		##
		## Argument parser from CLI
		self.parser = argparse.ArgumentParser(prog='scalehd', description='ScaleHD: Automated DNA micro-satellite genotyping.')
		self.parser.add_argument('-v', '--verbose', help='Verbose output mode. Setting this flag enables verbose output. Default: off.', action='store_true')
		input_group = self.parser.add_mutually_exclusive_group(required=True)
		input_group.add_argument('-i', '--input', help='Input data. Path to a single *.sam input file.', nargs=1)
		input_group.add_argument('-b', '--batch', help='Input batch. Folder of multiple .sam files. For FastQ processing, use -c.', nargs=1)
		input_group.add_argument('-c', '--config', help='Pipeline config. Specify a directory to your ArgumentConfig.xml file.', nargs=1)
		self.parser.add_argument('-t', '--threads', help='Thread utilisation. Typically only alters third party alignment performance. Default: system max.', type=int, choices=xrange(1, THREADS+1), default=THREADS)
		self.parser.add_argument('-o', '--output', help='Output path. Specify a directory you wish output to be directed towards.', metavar='output', nargs=1, required=True)
		self.args = self.parser.parse_args()

		##
		## Set verbosity for CLI output
		if self.args.verbose:
			log.basicConfig(format='%(message)s', level=log.DEBUG)
			log.info('{}{}{}{}'.format(clr.bold, 'shd__ ', clr.end, 'ScaleHD: Automated DNA micro-satellite genotyping.'))
			log.info('{}{}{}{}'.format(clr.bold, 'shd__ ', clr.end, 'alastair.maxwell@glasgow.ac.uk\n'))
		else:
			log.basicConfig(format='%(message)s')

		##
		## Check inputs, generate outputs
		if sanitise_inputs(self.args):
			log.error('{}{}{}{}'.format(clr.red, 'shd__ ', clr.end, 'Error with specified input(s) configuration. Exiting.'))
			sys.exit(2)
		self.instance_rundir = sanitise_outputs(self.args.output)

		##
		## Set up config dictionary of all params.
		## if -c used, read from XML. Else, use 'defaults' in set_params().
		script_path = os.path.dirname(__file__)
		if not self.args.config:
			self.instance_params = self.set_prediction_params()
		else:
			self.configfile = self.args.config[0]
			self.instance_params = ConfigReader(script_path, self.configfile)

		##
		## Check third party libraries before continuing
		if initialise_libraries(self.instance_params):
			log.error('{}{}{}{}'.format(clr.red, 'shd__ ', clr.end, 'Detected missing library from system/$PATH. Exiting.'))
			sys.exit(2)
		else:
			log.info('{}{}{}{}'.format(clr.green, 'shd__ ', clr.end, 'Required libraries present, assuming OK!\n'))

		##
		## Depending on input mode, direct flow of functions
		## -i == single file, pass to class
		## -b == multiple files, loop files to class
		## -c == config, do as config parsed flags
		if not self.args.config: self.sam_workflow()
		else: self.instance_workflow()
		log.info('{}{}{}{}'.format(clr.green, 'shd__ ', clr.end, 'ScaleHD pipeline completed; exiting.'))

	@staticmethod
	def set_prediction_params():

		param_dict = {'quality_control':'False',
					  'sequence_alignment':'False',
					  'genotype_prediction':'True',
					  'decision_function_shape':'ovr',
					  'probability_estimate':'True',
					  'max_iteration':'-1',}
		return param_dict

	def sam_workflow(self):
		##TODO rework sam-flow after retooling -config flow
		print self.instance_rundir

	def instance_workflow(self):

		##
		## Config generics
		instance_inputdata = self.instance_params.config_dict['@data_dir']

		##
		## Pre-stage: check for compressed data, extract
		if not extract_data(instance_inputdata):
			log.error('{}{}{}{}'.format(clr.red, 'shd__ ', clr.end, 'Error during file extraction. Please check your input data.'))

		##
		## Executing the workflow for this SHD instance
		## Ensure there are even amount of files for forward/reverse sequence pairings
		data_pairs = sequence_pairings(instance_inputdata, self.instance_rundir)
		for i in range(len(data_pairs)):
			log.info('{}{}{}{}{}{}{}'.format(clr.bold, 'shd__ ', clr.end, 'Processing sequence pair: ', str(i+1), '/', str(len(data_pairs))))
			for sequence_label, sequencepair_data in data_pairs[i].iteritems():

				##
				## For the Sequence Pair dictionary we're currently in
				## create object of the desired stage paths..
				qc_path = sequencepair_data[2]
				align_path = sequencepair_data[3]
				predict_path = sequencepair_data[4]

				##
				## Stage 1: QC and subflags
				seq_qc_flag = self.instance_params.config_dict['instance_flags']['@quality_control']
				seq_qc_dmpx = self.instance_params.config_dict['dmplex_flags']['@demultiplex_data']
				seq_qc_trim = self.instance_params.config_dict['trim_flags']['@trim_data']

				if seq_qc_flag == 'True':
					log.info('{}{}{}{}'.format(clr.yellow,'shd__ ',clr.end,'Executing sequence quality control workflow..'))
					if seq_qc.SeqQC(sequence_label, sequencepair_data, qc_path, 'valid', self.instance_params):
						if seq_qc_dmpx == 'True':
							log.info('{}{}{}{}'.format(clr.bold, 'shd__ ', clr.end, 'Initialising demultiplexing.'))
							seq_qc.SeqQC(sequence_label, sequencepair_data, qc_path, 'dmpx', self.instance_params)
							log.info('{}{}{}{}'.format(clr.green, 'shd__ ', clr.end, 'Demultiplexing complete!'))
						if seq_qc_trim == 'True':
							log.info('{}{}{}{}'.format(clr.bold, 'shd__ ', clr.end, 'Initialising trimming.'))
							seq_qc.SeqQC(sequence_label, sequencepair_data, qc_path, 'trim', self.instance_params)
							log.info('{}{}{}{}'.format(clr.green, 'shd__ ', clr.end, 'Trimming complete!'))

				##
				## Stage 2: Alignment flags
				alignment_flag = self.instance_params.config_dict['instance_flags']['@sequence_alignment']
				if alignment_flag == 'True':
					log.info('{}{}{}{}'.format(clr.yellow,'shd__ ',clr.end,'Executing alignment workflow..'))
					align.SeqAlign(sequence_label, sequencepair_data, align_path, self.instance_params)
					log.info('{}{}{}{}'.format(clr.green, 'shd__ ', clr.end, 'Sequence alignment workflow complete!'))

				##
				## Stage 3: Genotyping flags
				genotyping_flag = self.instance_params.config_dict['instance_flags']['@genotype_prediction']
				if genotyping_flag == 'True':
					log.info('{}{}{}{}'.format(clr.yellow,'shd__ ',clr.end,'Executing genotyping workflow..'))
					log.info('{}{}{}{}'.format(clr.green,'shd__ ',clr.end,'Genotyping workflow complete!'))

				##
				## Finished all desired stages for this file pair, inform user if -v
				log.info('{}{}{}{}'.format(clr.green, 'shd__ ', clr.end, 'Sequence pair workflow complete!\n'))










































def main():
	try:
		ScaleHD()
	except KeyboardInterrupt:
		log.error('{}{}{}{}'.format(clr.red,'shd__ ',clr.end,'Fatal: Keyboard Interrupt detected. Exiting.'))
		sys.exit(2)