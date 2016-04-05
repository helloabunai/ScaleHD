#/usr/bin/python
import os
import sys
import argparse
import pkg_resources
import logging as log
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
		self.hdgenerics_dat = pkg_resources.resource_filename(__name__, 'train/placeholder.csv')
		self.hdgenerics_rst = pkg_resources.resource_filename(__name__, 'train/placeholder.rst')

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
		initialise_libraries()
		sanitise_inputs(self.args)
		self.instance_rundir = sanitise_outputs(self.args.output)

		##
		## Set up config dictionary of all params.
		## if -c used, read from XML. Else, use 'defaults' in set_params().
		script_path = os.path.dirname(__file__)
		if not self.args.config:
			self.instance_params = self.set_params()
		else:
			self.configfile = self.args.config[0]
			self.instance_params = ConfigReader(script_path, self.configfile)

		##
		## Depending on input mode, direct flow of functions
		## -i == single file, pass to class
		## -b == multiple files, loop files to class
		## -c == config, do as config parsed flags
		if self.args.input:
			predict.SeqPredict(self.args.input[0])
		if self.args.batch:
			for sam_file in self.args.batch[0]:
				predict.SeqPredict(sam_file)
		if self.args.config:
			self.instance_workflow()

	@staticmethod
	def set_params():

		##todo default params for gtype only

		param_dict = {'a':'',
					  'b':'',
					  'c':'',}

		return param_dict

	def instance_workflow(self):

		##
		## Config generics
		instance_inputdata = self.instance_params.config_dict['@data_dir']
		processed_files = []

		##
		## Pre-stage: check for compressed data, extract
		extract_data(instance_inputdata)

		##
		## Stage 1: QC and subflags
		seq_qc_flag = self.instance_params.config_dict['instance_flags']['@quality_control']
		seq_qc_dmpx = self.instance_params.config_dict['dmplex_flags']['@demultiplex_data']
		seq_qc_trim = self.instance_params.config_dict['trim_flags']['@trim_data']

		if seq_qc_flag == 'True':
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

		##
		## Stage 2: Alignment flags
		alignment_flag = self.instance_params.config_dict['instance_flags']['@sequence_alignment']

		if alignment_flag == 'True':
			print 'alignment true'

		##
		## Stage 3: Genotyping flags
		genotyping_flag = self.instance_params.config_dict['instance_flags']['@genotype_prediction']

		if genotyping_flag == 'True':
			print 'genotyping true'
























def main():
	BaseCamp()