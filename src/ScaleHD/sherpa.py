#/usr/bin/python
import os
import sys
import argparse
import pkg_resources
import logging as log
import glob
from multiprocessing import cpu_count

##
## Backend junk
from backpack import ConfigReader
from backpack import Colour as clr
from backpack import initialise_libraries
from backpack import sanitise_inputs
from backpack import extract_data
from backpack import sanitise_outputs

##
## Package stages
from . import seq_qc
from . import align
from .align.alignment import get_repeat_distributions
from . import predict

##
## Globals
VERBOSE = False
LOGGING = True
THREADS = cpu_count()
DEF_OUT = os.path.join(os.path.expanduser('~'),'ScaleHD')

class BaseCamp:

	def __init__(self):

		##
		## Package data
		self.hdgenerics_dat = pkg_resources.resource_filename(__name__, 'train/hd_generic_model.csv')
		self.hdgenerics_rst = pkg_resources.resource_filename(__name__, 'train/hd_generic_model.rst')

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
		sanitise_inputs(self.args)
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
		initialise_libraries(self.instance_params)

		##
		## Depending on input mode, direct flow of functions
		## -i == single file, pass to class
		## -b == multiple files, loop files to class
		## -c == config, do as config parsed flags
		if self.args.input or self.args.batch:
			self.sam_workflow()
		if self.args.config:
			self.instance_workflow()

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

		##
		## Build list of sam_files
		sam_files = []
		if self.args.input:
			sam_files.append(self.args.input[0])
		if self.args.batch:
			for target_sam in glob.glob(os.path.join(self.args.batch[0], '*')):
				sam_files.append(target_sam)

		##
		## For cases where the user has specified input with -i/-b instead of -c
		## Some CLI feedback if verbose mode is on
		genotyping_flag = self.instance_params['genotype_prediction']
		if genotyping_flag == 'True':
			log.info('{}{}{}{}'.format(clr.bold, 'shd__ ', clr.end, 'Executing genotype prediction workflow..'))
			predict.SeqPredict(self.instance_rundir, self.instance_params, self.hdgenerics_dat, self.hdgenerics_rst, assembly_files=sam_files)
			log.info('{}{}{}{}'.format(clr.green, 'shd__ ', clr.end, 'Genotype prediction workflow complete!'))
			log.info('{}{}{}{}'.format(clr.bold, 'shd__ ', clr.end, 'Pipeline complete. Goodbye!'))

	def instance_workflow(self):

		##
		## Config generics
		instance_inputdata = self.instance_params.config_dict['@data_dir']
		processed_files = []
		distribution_files = []

		##
		## Pre-stage: check for compressed data, extract
		extract_data(instance_inputdata)

		##
		## Stage 1: QC and subflags
		seq_qc_flag = self.instance_params.config_dict['instance_flags']['@quality_control']
		seq_qc_dmpx = self.instance_params.config_dict['dmplex_flags']['@demultiplex_data']
		seq_qc_trim = self.instance_params.config_dict['trim_flags']['@trim_data']

		if seq_qc_flag == 'True':
			log.info('{}{}{}{}'.format(clr.bold,'shd__ ',clr.end,'Executing sequence quality control workflow..'))
			if seq_qc.SeqQC(instance_inputdata, self.instance_rundir, 'valid', self.instance_params):
				if seq_qc_dmpx == 'True':
					log.info('{}{}{}{}'.format(clr.bold, 'shd__ ', clr.end, 'Initialising demulitplexing.'))
					seq_qc.SeqQC(instance_inputdata, self.instance_rundir, 'dmpx', self.instance_params)
					log.info('{}{}{}{}'.format(clr.green, 'shd__ ', clr.end, 'Demultiplexing complete!'))
				if seq_qc_trim == 'True':
					log.info('{}{}{}{}'.format(clr.bold, 'shd__ ', clr.end, 'Initialising trimming.'))
					seq_qc.SeqQC(instance_inputdata, self.instance_rundir, 'trim', self.instance_params)
					log.info('{}{}{}{}'.format(clr.green, 'shd__ ', clr.end, 'Trimming complete!'))
				else:
					log.error('{}{}{}{}'.format(clr.red, 'shd__ ', clr.end, 'Flow: SeqQC=True, yet neither Trim/Demultiplex=True.'))
					sys.exit(2)
			processed_files = seq_qc.SeqQC.getProcessedFiles()
			log.info('{}{}{}{}'.format(clr.green,'shd__ ',clr.end,'Sequence quality control workflow complete!\n'))
		else:
			## When files are dmpx/trimmed, saved to diff dir and returned in list
			## If no processing done in SeqQC, just loop over input and add files to list anyway
			## Continue to Align stage with these
			for unprocessed_fqfile in glob.glob(os.path.join(instance_inputdata, '*')):
				processed_files.append(unprocessed_fqfile)

		##
		## Stage 2: Alignment flags
		alignment_flag = self.instance_params.config_dict['instance_flags']['@sequence_alignment']
		if alignment_flag == 'True':
			log.info('{}{}{}{}'.format(clr.bold, 'shd__ ', clr.end, 'Executing sequence alignment workflow..'))
			align.SeqAlign(self.args, processed_files, self.instance_rundir, self.instance_params)
			distribution_files = get_repeat_distributions()
			log.info('{}{}{}{}'.format(clr.green, 'shd__ ', clr.end, 'Sequence alignment workflow complete!\n'))

		##
		## Stage 3: Genotyping flags
		genotyping_flag = self.instance_params.config_dict['instance_flags']['@genotype_prediction']

		if genotyping_flag == 'True':
			log.info('{}{}{}{}'.format(clr.bold, 'shd__ ', clr.end, 'Executing genotype prediction workflow..'))
			predict.SeqPredict(self.instance_rundir, self.instance_params, self.hdgenerics_dat, self.hdgenerics_rst, distribution_files=distribution_files)
			log.info('{}{}{}{}'.format(clr.green, 'shd__ ', clr.end, 'Genotype prediction workflow complete!'))
			log.info('\n{}{}{}{}'.format(clr.bold, 'shd__ ', clr.end, 'Pipeline complete. Goodbye!'))

def main():
	try:
		BaseCamp()
	except KeyboardInterrupt:
		log.error('{}{}{}{}'.format(clr.red,'shd__ ',clr.end,'Fatal: Keyboard Interrupt detected. Exiting.'))
		sys.exit(2)