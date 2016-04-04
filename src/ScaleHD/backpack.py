import string
import os
import sys
import glob
import datetime
import subprocess
import logging as log
import numpy as np
from collections import defaultdict
from xml.etree import cElementTree
from lxml import etree

class Colour:

	def __init__(self):
		pass

	purple = '\033[95m'
	cyan = '\033[96m'
	darkcyan = '\033[36m'
	blue = '\033[94m'
	green = '\033[92m'
	yellow = '\033[93m'
	red = '\033[91m'
	bold = '\033[1m'
	underline = '\033[4m'
	end = '\033[0m'

class ScaleHDException:
	def __init__(self, err_str):
		pass

class ConfigReader(object):

	"""
	The configuration file reader.
	Opens a configuration file, and if valid, converts the parameters within the file to a dictionary object,
	reader to be viewed through accessing the config_dict variable.
	"""

	def __init__(self, scriptdir, config_filename=None):

		##
		## Instance variables
		self.scriptdir = scriptdir
		self.config_filename = config_filename
		self.dtd_filename = scriptdir + "/config/config.dtd"

		##
		## Check for configuration file (just incase)
		if self.config_filename is None:
			log.error("No configuration file specified!")
		else:
			self.config_file = etree.parse(self.config_filename)

		##
		## Check config vs dtd, parse info to dictionary, validate vs ruleset
		self.validate_against_dtd()
		self.set_dictionary()
		self.validate_config()

	def validate_against_dtd(self):

		"""
		Validate input config against DTD ruleset
		i.e. confirms conformation of XML structure
		"""

		##
		## Open > etree.DTD object
		dtd_file = open(self.dtd_filename, 'r')
		dtd_object = etree.DTD(dtd_file)

		##
		## If validation fails, close the object (memory) and raise an error
		if not dtd_object.validate(self.config_file):
			dtd_file.close()
			log.error("DTD validation failure {0}: {1}".format(self.config_filename, dtd_object.error_log.filter_from_errors()[0]))
			sys.exit(2)
		dtd_file.close()

	def set_dictionary(self):

		"""
		Takes the now validated XML and extracts information from the tree into
		a python dictionary {key: value}. This dictionary will be used for variables
		within the pipeline. Recursion adapted from http://stackoverflow.com/a/9286702
		"""
		def recursive_generation(t):

			d = {t.tag: {} if t.attrib else None}
			children = list(t)

			##
			## If list was populated, create dictionary, Append keys
			if children:
				dd = defaultdict(list)

				for dc in map(recursive_generation, children):
					for k, v in dc.iteritems():
						dd[k].append(v)
				d = {t.tag: {k: v[0] if len(v) == 1 else v for k, v in dd.iteritems()}}

			##
			## Values for key
			if t.attrib:
				d[t.tag].update(('@' + k, v) for k, v in t.attrib.iteritems())

			if t.text:
				text = t.text.strip()
				if children or t.attrib:
					if text:
						d[t.tag]['#text'] = text
				else:
					d[t.tag] = text
			return d

		##
		## Takes the formatted xml doc, puts through generator, returns dictionary
		string_repr = etree.tostring(self.config_file, pretty_print=True)
		element_tree = cElementTree.XML(string_repr)

		self.config_dict = recursive_generation(element_tree)
		self.config_dict = self.config_dict[self.config_dict.keys()[0]]

	def validate_config(self):

		"""
		Method which validates the configuration file's contents.
		If all pass, guarantees that the settings dictionary is full of valid settings!
		"""
		trigger = False

		##
		## Main configuration instance settings

		data_directory = self.config_dict['@data_dir']
		if not os.path.exists(data_directory):
			log.error('{}{}{}{}'.format(Colour.red, 'shd__', Colour.end, 'XML Config: Specified data directory could not be found.'))
			trigger = True
		for fqfile in glob.glob(os.path.join(data_directory, '*')):
			if not (fqfile.endswith('.fq') or fqfile.endswith('.fastq')):
				log.error('{}{}{}{}'.format(Colour.red, 'shd__ ', Colour.end, 'XML Config: Non FastQ data detected in specified input directory.'))
				trigger = True
		reference_directory = self.config_dict['@reference_file']
		if not os.path.isfile(reference_directory):
			log.error('{}{}{}{}'.format(Colour.red, 'shd__ ', Colour.end, 'XML Config: Specified reference file could not be found.'))
			trigger = True
		if not (reference_directory.endswith('.fa') or reference_directory.endswith('.fas')):
			log.error('{}{}{}{}'.format(Colour.red, 'shd__ ', Colour.end, 'XML Config: Specified reference file is not an fa/fas file.'))
			trigger = True

		##
		## Instance flag settings

		sequence_qc_flag = self.config_dict['instance_flags']['@quality_control']
		if not (sequence_qc_flag == 'True' or sequence_qc_flag == 'False'):
			log.error('{}{}{}{}'.format(Colour.red, 'shd__ ', Colour.end, 'XML Config: Sequence Quality control flag is not set to True/False.'))
			trigger = True
		alignment_flag = self.config_dict['instance_flags']['@sequence_alignment']
		if not (alignment_flag == 'True' or alignment_flag == 'False'):
			log.error('{}{}{}{}'.format(Colour.red, 'shd__ ', Colour.end, 'XML Config: Sequence Alignment flag is not set to True/False.'))
			trigger = True
		genotype_flag = self.config_dict['instance_flags']['@genotype_prediction']
		if not (genotype_flag == 'True' or genotype_flag == 'False'):
			log.error('{}{}{}{}'.format(Colour.red, 'shd__ ', Colour.end, 'XML Config: Genotype Prediction control flag is not True/False.'))
			trigger = True

		##
		## Demultiplexing flag settings
		if sequence_qc_flag:
			demultiplex_flag = self.config_dict['dmplex_flags']['@demultiplex_data']
			if not (demultiplex_flag == 'True' or demultiplex_flag == 'False'):
				log.error('{}{}{}{}'.format(Colour.red, 'shd__ ', Colour.end, 'XML Config: Demultiplexing flag is not True/False.'))
				trigger = True
			barcode_file = self.config_dict['dmplex_flags']['@barcode_file']
			if not os.path.isfile(barcode_file):
				log.error('{}{}{}{}'.format(Colour.red, 'shd__ ', Colour.end, 'XML Config: Specified barcode file could not be found.'))
				trigger = True
			demultiplex_mismatch = self.config_dict['dmplex_flags']['@max_mismatch']
			if not demultiplex_mismatch.isdigit():
				log.error('{}{}{}{}'.format(Colour.red, 'shd__ ', Colour.end, 'XML Config: Specified barcode mismatch integer is invalid.'))
				trigger = True

		##
		## Trimming flag settings
		if sequence_qc_flag:
			trimming_flag = self.config_dict['trim_flags']['@trim_data']
			if not (trimming_flag == 'True' or trimming_flag == 'False'):
				log.error('{}{}{}{}'.format(Colour.red, 'shd__ ', Colour.end, 'XML Config: Trimming flag is not True/False.'))
				trigger = True
			trimming_type = self.config_dict['trim_flags']['@trim_type']
			if not (trimming_type == 'Quality' or trimming_type	== 'Adapter'):
				log.error('{}{}{}{}'.format(Colour.red, 'shd__  ', Colour.end, 'XML Config: Trimming type is not Quality/Adapter.'))
				trigger = True
			quality_threshold = self.config_dict['trim_flags']['@quality_threshold']
			if not quality_threshold.isdigit():
				log.error('{}{}{}{}'.format(Colour.red, 'shd__ ', Colour.end, 'XML Config: Specified quality threshold integer is invalid.'))
				trigger = True
			elif not int(quality_threshold) in range(0,38):
				log.error('{}{}{}{}'.format(Colour.red, 'shd__ ', Colour.end, 'XML Config: Specified quality threshold integer out of range (0-38).'))
				trigger = True
			trim_adapters = ['-a','-g','a$','-g^','-b']
			adapter_flag = self.config_dict['trim_flags']['@adapter_flag']
			if not (adapter_flag in trim_adapters):
				log.error('{}{}{}{}'.format(Colour.red, 'shd__ ', Colour.end, 'XML Config: Specified trimming adapter not valid selection.'))
				trigger = True
			trim_adapter_base = ['A','G','C','T']
			adapter_sequence = self.config_dict['trim_flags']['@adapter']
			for charbase in adapter_sequence:
				if charbase not in trim_adapter_base:
					log.error('{}{}{}{}'.format(Colour.red, 'shd__ ', Colour.end, 'XML Config: Invalid character detected in adapter sequence.'))
					trigger = True

		##
		## Alignment flag settings
		if alignment_flag:
			extension_threshold = self.config_dict['alignment_flags']['@extension_threshold']
			if not extension_threshold.isdigit():
				log.error('{}{}{}{}'.format(Colour.red, 'shd__ ', Colour.end, 'XML Config: Specified extension threshold integer is invalid.'))
				trigger = True
			seed_size = self.config_dict['alignment_flags']['@seed_size']
			if not seed_size.isdigit():
				log.error('{}{}{}{}'.format(Colour.red, 'shd__ ', Colour.end, 'XML Config: Specified seed size integer is invalid.'))
				trigger = True
			align_mismatch = self.config_dict['alignment_flags']['@align_mismatch']
			if not align_mismatch.isdigit():
				log.error('{}{}{}{}'.format(Colour.red, 'shd__ ', Colour.end, 'XML Config: Specified align mismatch integer is invalid.'))
				trigger = True
			substring_interval_start = self.config_dict['alignment_flags']['@substr_interval_start']
			if not int(substring_interval_start) in range(0,2):
				log.error('{}{}{}{}'.format(Colour.red, 'shd__ ', Colour.end, 'XML Config: Seed Substring Interval (start) integer is out of range (0,1).'))
				trigger = True
			substring_interval_end = self.config_dict['alignment_flags']['@substr_interval_end']
			if float(substring_interval_end) in np.arange(0.5,2.55,0.05):
				log.error('{}{}{}{}'.format(Colour.red, 'shd__ ', Colour.end, 'XML Config: Seed Substring Interval (end) float is out of range (0.50, 2.50).'))
				trigger = True

		##
		## Genotype prediction flag settings
		if genotype_flag:
			sample_flag = self.config_dict['prediction_flags']['@sample_flag']
			if not (sample_flag == 'True' or sample_flag == 'False'):
				log.error('{}{}{}{}'.format(Colour.red, 'shd__  ', Colour.end, 'XML Config: Specified sample_flag is not True/False.'))
				trigger = True

		if trigger:
			log.error('{}{}{}{}'.format(Colour.red, 'shd__ ', Colour.end, 'XML Config: Failure, exiting.'))
			sys.exit(2)
		else:
			log.info('{}{}{}{}'.format(Colour.green, 'shd__ ', Colour.end, 'XML Config: Successful parsing!'))


def parse_boolean(boolean_value):

	"""
	Given a string (boolean_value), returns a boolean value representing the string contents.
	For example, a string with 'true', 't', 'y' or 'yes' will yield True.
	"""

	boolean_value = string.lower(boolean_value) in ('yes', 'y', 'true', 't', '1')
	return boolean_value

def empty_string_check(string, raise_exception=True):

	"""
	Simple check to see if the string provided by parameter string is empty. False indicates the string is NOT empty.
	Parameter raise_exception determines if a ValueError exception should be raised if the string is empty.
	If raise_exception is False and the string is empty, True is returned.
	"""

	if string != '':
		return False
	if raise_exception:
		raise ValueError("Empty string detected!")
	return True

def sanitise_inputs(parsed_arguments):

	"""
	Utilises filesystem_exists_check and check_input_files
	if either return false, path is invalid or unsupported files present
	so, quit
	"""

	if parsed_arguments.input:
		if not filesystem_exists_check(parsed_arguments.input[0]):
			log.error('{}{}{}{}'.format(Colour.red, 'shd__  ', Colour.end, 'Specified input file could not be found.'))
			sys.exit(2)
		for samfile in parsed_arguments.input:
			if not check_input_files('.sam',samfile):
				log.error('{}{}{}{}'.format(Colour.red, 'shd__  ', Colour.end, 'Specified input file is not a SAM file.'))
				sys.exit(2)

	if parsed_arguments.batch:
		if not filesystem_exists_check(parsed_arguments.batch[0]):
			log.error('{}{}{}{}'.format(Colour.red, 'shd__  ', Colour.end, 'Specified batch folder could not be found.'))
			sys.exit(2)
		for samfile in parsed_arguments.batch:
			if not check_input_files('.sam',samfile):
				log.error('{}{}{}{}'.format(Colour.red, 'shd__  ', Colour.end, 'Specified batch folder contains non SAM files.'))
				sys.exit(2)

	if parsed_arguments.config:
		if not filesystem_exists_check(parsed_arguments.config[0]):
			log.error('{}{}{}{}'.format(Colour.red, 'shd__  ', Colour.end, 'Specified config file could not be found.'))
			sys.exit(2)
		for xmlfile in parsed_arguments.config:
			if not check_input_files('.xml',xmlfile):
				log.error('{}{}{}{}'.format(Colour.red, 'shd__  ', Colour.end, 'Specified config file is not an XML file.'))
				sys.exit(2)

def filesystem_exists_check(path, raise_exception=True):

	"""
	Checks to see if the path, specified by parameter path, exists. Can be either a directory or file.
	If the path exists, True is returned. If the path does not exist, and raise_exception is set to True,
	an IOError is raised - else False is returned.
	"""

	if os.path.lexists(path):
		return True
	if raise_exception:
		log.error('{}{}{}{}'.format(Colour.red,'shd__ ',Colour.end,'Specified -i/-b/-c path could not be found.'))
	return False

def check_input_files(input_format, input_file, raise_exception=True):

	if input_file.endswith(input_format):
		return True
	if raise_exception:
		log.error('{}{}{}{}'.format(Colour.red,'shd__ ',Colour.end,'Unrecognised file format found in -i/-b/-c path.'))
	return False

def initialise_libraries():

	trigger = False

	##
	## Check for cutadapt
	which_cutadapt = subprocess.Popen(['which', 'cutadapt'],stdout=subprocess.PIPE,stderr=subprocess.PIPE)
	cutadapt_dir = which_cutadapt.communicate()[0]
	which_cutadapt.wait()
	if not 'cutadapt' in cutadapt_dir:
		log.error('{}{}{}{}'.format(Colour.red, 'shd__ ', Colour.end, 'Missing library: Cutadapt. Not installed or not on $PATH.'))
		trigger = True

	##
	## Check for sabre
	which_sabre = subprocess.Popen(['which', 'sabre'],stdout=subprocess.PIPE,stderr=subprocess.PIPE)
	sabre_dir = which_sabre.communicate()[0]
	which_sabre.wait()
	if not 'sabre' in sabre_dir:
		log.error('{}{}{}{}'.format(Colour.red, 'shd__ ', Colour.end, 'Missing library: Sabre. Not installed or not on $PATH.'))
		trigger = True

	##
	## Check for bowtie2
	which_bowtie = subprocess.Popen(['which', 'bowtie2'],stdout=subprocess.PIPE,stderr=subprocess.PIPE)
	bowtie_dir = which_bowtie.communicate()[0]
	which_bowtie.wait()
	if not 'bowtie2' in bowtie_dir:
		log.error('{}{}{}{}'.format(Colour.red, 'shd__ ', Colour.end, 'Missing library: Bowtie2. Not installed or not on $PATH.'))
		trigger = True

	##
	## Pass error
	if trigger:
		log.error('{}{}{}{}'.format(Colour.red, 'shd__ ', Colour.end, 'Cannot progress without all third party libraries. Exiting.'))
		sys.exit(2)
	else:
		log.info('{}{}{}{}'.format(Colour.green, 'shd__ ', Colour.end, 'All libraries present. Assume OK!'))

def sanitise_outputs(output_argument):

	output_root = output_argument[0]

	## Ensures root output is a real directory
	## Generates folder name based on date (for run ident)
	date = datetime.date.today().strftime('%d-%m-%Y')
	time = datetime.datetime.now().strftime('%H%M%S')
	today = date + '-' + time

	## If the user specified root doesn't exist, make it
	## Then make the run directory for datetime
	if not os.path.exists(output_root):
		log.info('{}{}{}{}'.format(Colour.bold, 'shd__ ', Colour.end, 'Creating output root... '))
		os.mkdir(output_root)
	run_dir = output_root + '/ScaleHDRun_' + today
	log.info('{}{}{}{}'.format(Colour.bold, 'shd__ ', Colour.end, 'Creating instance run directory.. '))
	os.mkdir(run_dir)

	## Inform user it's all gonna be okaaaayyyy
	log.info('{}{}{}{}'.format(Colour.green, 'shd__ ', Colour.end, 'Output directories OK!'))
	return run_dir