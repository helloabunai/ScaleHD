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
		if not self.args.config: self.assembly_workflow()
		else: self.sequence_workflow()
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

	def assembly_workflow(self):

		##
		## Input path, create pairs of data in said path
		instance_inputdata = self.args.batch[0]
		assembly_pairs = sequence_pairings(instance_inputdata, self.instance_rundir, 'assembly')

		##
		## Executing the workflow for this SHD instance
		for i in range(len(assembly_pairs)):
			log.info('{}{}{}{}{}{}{}'.format(clr.bold, 'shd__ ', clr.end, 'Processing assembly pair: ', str(i+1), '/', str(len(assembly_pairs))))
			for assembly_label, assembly_data in assembly_pairs[i].iteritems():

				##
				## Required data to process
				forward_assembly = assembly_data[0]
				reverse_assembly = assembly_data[1]
				predict_path = assembly_data[2]
				instance_params = self.set_prediction_params()

				##
				## Specific paths to pass to distribution scraper
				## In instance_workflow these would've been created for alignment, but since we don't align here
				## They have to be made in this location instead
				forward_filename = forward_assembly.split('/')[-1].split('.')[0] ##absolutely_disgusting.jpg
				reverse_filename = reverse_assembly.split('/')[-1].split('.')[0] ##absolutely_disgusting.jpg
				forward_path = os.path.join(predict_path,forward_filename)
				reverse_path = os.path.join(predict_path,reverse_filename)
				if not os.path.exists(forward_path): os.makedirs(forward_path)
				if not os.path.exists(reverse_path): os.makedirs(reverse_path)

				##
				## Pre stage -- extract distributions from input sam files
				## Update assembly_data list; replacing forward/reverse assemblys with respective repeat distributions
				log.info('{}{}{}{}'.format(clr.yellow,'shd__ ',clr.end,'Extracting repeat distributions from pre-assembled data..'))
				assembly_data[0] = align.SeqAlign.extract_repeat_distributions(assembly_label,forward_path,forward_assembly)
				assembly_data[1] = align.SeqAlign.extract_repeat_distributions(assembly_label,reverse_path,reverse_assembly)
				log.info('{}{}{}{}'.format(clr.green,'shd__ ',clr.end,'Repeat distribution extraction complete!'))

				## Prediction Stage
				log.info('{}{}{}{}'.format(clr.yellow,'shd__ ',clr.end,'Executing genotyping workflow..'))
				predict.GenotypePrediction(assembly_data, predict_path, self.training_data, instance_params)
				log.info('{}{}{}{}'.format(clr.green,'shd__ ',clr.end,'Genotyping workflow complete!'))

				##
				## Finished all desired stages for this file pair, inform user if -v
				log.info('{}{}{}{}'.format(clr.green, 'shd__ ', clr.end, 'Assembly pair workflow complete!\n'))

	def sequence_workflow(self):

		##
		## Config generics
		instance_inputdata = self.instance_params.config_dict['@data_dir']


		##
		## Pre-stage: check for compressed data, extract
		if not extract_data(instance_inputdata):
			log.error('{}{}{}{}'.format(clr.red, 'shd__ ', clr.end, 'Error during file extraction. Please check your input data.'))

		##
		## If the user wants to align, it makes architectural sense to index references (which will be used in all pairs)
		## beforehand, so as to not repeat the computation for each cycle; thus the indexes are created outside of the
		## main "workflow" so to speak :: check flag and then get reference indexes
		reference_indexes = []
		if self.instance_params.config_dict['instance_flags']['@sequence_alignment']:
			log.info('{}{}{}{}'.format(clr.bold,'shd__ ',clr.end,'Indexing reference(s) before initialising sample pair cycle..'))
			index_path = os.path.join(self.instance_rundir,'Indexes')
			if not os.path.exists(index_path): os.makedirs(index_path)

			forward_reference = self.instance_params.config_dict['@forward_reference']
			reverse_reference = self.instance_params.config_dict['@reverse_reference']

			forward_index = align.ReferenceIndex(forward_reference, index_path).getIndexPath()
			reverse_index = align.ReferenceIndex(reverse_reference, index_path).getIndexPath()
			reference_indexes = [forward_index, reverse_index]

		##
		## Executing the workflow for this SHD instance
		## Ensure there are even amount of files for forward/reverse sequence pairings
		data_pairs = sequence_pairings(instance_inputdata, self.instance_rundir, 'sequence')
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
					if seq_qc.SeqQC(sequencepair_data, qc_path, 'valid', self.instance_params):
						if seq_qc_dmpx == 'True':
							log.info('{}{}{}{}'.format(clr.bold, 'shd__ ', clr.end, 'Initialising demultiplexing.'))
							seq_qc.SeqQC(sequencepair_data, qc_path, 'dmpx', self.instance_params)
							log.info('{}{}{}{}'.format(clr.green, 'shd__ ', clr.end, 'Demultiplexing complete!'))
						if seq_qc_trim == 'True':
							log.info('{}{}{}{}'.format(clr.bold, 'shd__ ', clr.end, 'Initialising trimming.'))
							seq_qc.SeqQC(sequencepair_data, qc_path, 'trim', self.instance_params)
							log.info('{}{}{}{}'.format(clr.green, 'shd__ ', clr.end, 'Trimming complete!'))

				##
				## Stage 2: Alignment flags
				alignment_flag = self.instance_params.config_dict['instance_flags']['@sequence_alignment']
				if alignment_flag == 'True':
					log.info('{}{}{}{}'.format(clr.yellow,'shd__ ',clr.end,'Executing alignment workflow..'))
					align.SeqAlign(sequence_label, sequencepair_data, align_path, reference_indexes, self.instance_params)
					log.info('{}{}{}{}'.format(clr.green, 'shd__ ', clr.end, 'Sequence alignment workflow complete!'))

				##
				## Stage 3: Genotyping flags
				genotyping_flag = self.instance_params.config_dict['instance_flags']['@genotype_prediction']
				if genotyping_flag == 'True':
					log.info('{}{}{}{}'.format(clr.yellow,'shd__ ',clr.end,'Executing genotyping workflow..'))
					predict.GenotypePrediction(sequencepair_data, predict_path, self.training_data, self.instance_params)
					log.info('{}{}{}{}'.format(clr.green,'shd__ ',clr.end,'Genotyping workflow complete!'))

				##
				## Scrape summary results from stages and write to 'master' output file



				##
				## Finished all desired stages for this file pair, inform user if -v
				log.info('{}{}{}{}'.format(clr.green, 'shd__ ', clr.end, 'Sequence pair workflow complete!\n'))


def main():
	try:
		ScaleHD()
	except KeyboardInterrupt:
		log.error('{}{}{}{}'.format(clr.red,'shd__ ',clr.end,'Fatal: Keyboard Interrupt detected. Exiting.'))
		sys.exit(2)