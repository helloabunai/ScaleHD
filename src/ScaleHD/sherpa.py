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
from backpack import sanitise_inputs

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
			log.info('{}{}{}{}'.format(clr.bold, 'shd__ ', clr.end, 'alastair.maxwell@glasgow.ac.uk'))
		else:
			log.basicConfig(format='%(message)s', level=log.NOTSET)

		##
		## Check inputs
		sanitise_inputs(self.args)

		##
		## Set up config dictionary of all params.
		## if -c used, read from XML. Else, use 'defaults' in set_params().
		script_path = os.path.dirname(__file__)
		if not self.args.config:
			self.instance_params = self.set_params()
		else:
			self.configfile = self.args.config[0]
			self.instance_params = ConfigReader(script_path, self.configfile)


	@staticmethod
	def set_params():

		param_dict = {'a':'',
					  'b':'',
					  'c':'',}

		return param_dict




























def main():
	BaseCamp()