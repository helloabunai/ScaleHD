import string
import os
import sys
import logging as log
from multiprocessing import cpu_count
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

		invalid_flags = {}

		##
		## Main configuration settings
		rawdata_dir = self.config_dict['@rawdata_directory']
		if not os.path.exists(rawdata_dir):
			invalid_flags['RawDataDirectory'] =  'User specified raw data directory is non existent!'

		rfnc_dir = self.config_dict['@rfnc_directory']
		if not rfnc_dir:
			pass
		else:
			if not os.path.isfile(rfnc_dir):
				invalid_flags['ReferenceDirectory'] = 'User specified reference file doesn\'t exist!'

		outp_dir = self.config_dict['@outp_directory']
		if not os.path.exists(outp_dir):
			os.makedirs(outp_dir)

		##
		## System_init settings
		disease_flag = self.config_dict['system_init']['@disease_flag']
		if not (disease_flag == 'DM1' or disease_flag == 'HD'):
			invalid_flags['DiseaseFlag'] = 'User specified disease_flag invalid! (valid: DM1, HD)'

		treatment_stage = self.config_dict['system_init']['@treatment_stage']
		if not (treatment_stage == 'True' or treatment_stage == 'False'):
			invalid_flags['TreatmentStage'] = 'User specified treatment_stage flag invalid! (valid: True, False)'

		alignment_stage = self.config_dict['system_init']['@alignment_stage']
		if not (alignment_stage == 'True' or alignment_stage == 'False'):
			invalid_flags['AlignmentStage'] = 'User specified alignment_stage flag invalid! (valid: True, False)'

		genotypes_stage = self.config_dict['system_init']['@genotypes_stage']
		if not (genotypes_stage == 'True' or genotypes_stage == 'False'):
			invalid_flags['GenotypesStage'] = 'User specified genotypes_stage flag invalid! (valid: True, False)'


		##
		## dmpx_exec settings (demultiplexing)
		dmpx_data = self.config_dict['dmpx_exec']['@dmpx_data']
		if not (dmpx_data == 'True' or dmpx_data == 'False'):
			invalid_flags['DMPXData'] = 'User specified dmpx_data flag is not valid! (valid: True, False)'

		dmpx_paired = self.config_dict['dmpx_exec']['@paired_read']
		if not (dmpx_paired == 'True' or dmpx_paired == 'False'):
			invalid_flags['DMPXPairedRead'] = 'User specified (dmpx) paired_read flag is not valid! (valid: True, False)'

		dmpx_barcode = self.config_dict['dmpx_exec']['@barcode_file']
		if not os.path.isfile(dmpx_barcode):
			invalid_flags['DMPXBarcodeFile'] = 'User specified barcode_file flag is not a valid file! Check directory/filename!'

		dmpx_maxmismatch = self.config_dict['dmpx_exec']['@max_mismatch']
		if dmpx_maxmismatch is '' or None:
			pass
		else:
			if not dmpx_maxmismatch.isdigit():
				invalid_flags['DMPXMaxMismatch'] = 'User specified max_mismatch flag is not an integer!'
			if int(dmpx_maxmismatch) not in range(0,20):
				invalid_flags['DMPXMaxMismatch'] = 'User specified max_mismatch flag is out of range! (0-20)'

		##
		## Trim_exec settings (Quality/Adapter trimming)
		trim_data = self.config_dict['trim_exec']['@trim_data']
		if not (trim_data == 'True' or trim_data == 'False'):
			invalid_flags['TrimData'] = 'User specified trim_data flag is not valid! (valid: True, False)'

		valid_trim = ['Adapter', 'Quality', 'Both', '']
		trim_flag = self.config_dict['trim_exec']['@trim_type']
		if not (trim_flag in valid_trim):
			invalid_flags['TrimFlag'] = 'User specified trim_flag invalid! (valid: Adapter, Quality, Both)'

		quality_flag = self.config_dict['trim_exec']['@quality_threshold']
		if not quality_flag.isdigit():
			invalid_flags['QualityFlag'] = 'User specified quality_threshold is not an integer!'
		elif not int(quality_flag) in range(0,38):
			invalid_flags['QualityFlag'] = 'User specified quality threshold outside of range(0, 38)'

		valid_adapter = ['-a', '-g', '-a$', '-g^', '-b', '']
		adapter_flag = self.config_dict['trim_exec']['@adapter_flag']
		if not (adapter_flag in valid_adapter):
			invalid_flags['AdapterFlag'] = 'User specified adapter_flag invalid! (valid: -a, -g, -a$, -g^, -b)'

		valid_base = ['A','G','C','T', '']
		adapter_seq = self.config_dict['trim_exec']['@adapter']
		for charbase in adapter_seq:
			if charbase not in valid_base:
				invalid_flags['AdapterSequence'] = 'User specified adapter is invalid (bases not ATGC)'

		##
		## Gen_ref settings (Reference Generator)
		generate_ref = self.config_dict['genref_exec']['@generate_ref']
		if not (generate_ref == 'True' or generate_ref == 'False'):
			invalid_flags['GenerateReference'] = 'User specified generate_ref flag is not valid! (valid: True, False)'

		fiveprime = self.config_dict['genref_exec']['@fiveprime']
		for charbase in fiveprime:
			if charbase not in valid_base:
				invalid_flags['FivePrime'] = 'User specified fiveprime flag is not \'ATGC\''

		min_cag = self.config_dict['genref_exec']['@min_cag']
		if not min_cag.isdigit():
			invalid_flags['MinCAG'] = 'User specified min_cag flag is not an integer!'

		max_cag = self.config_dict['genref_exec']['@max_cag']
		if not max_cag.isdigit():
			invalid_flags['MaxCAG'] = 'User specified max_cag flag is not an integer!'

		intervening = self.config_dict['genref_exec']['@intervening']
		for charbase in intervening:
			if charbase not in valid_base:
				invalid_flags['Intervening'] = 'User specified intervening flag is not \'ATGC\''

		min_ccg = self.config_dict['genref_exec']['@min_ccg']
		if not min_ccg.isdigit():
			invalid_flags['MinCCG'] = 'User specified min_ccg flag is not an integer!'

		max_ccg = self.config_dict['genref_exec']['@max_ccg']
		if not max_ccg.isdigit():
			invalid_flags['MaxCCG'] = 'User specified max_ccg flag is not an integer!'

		threeprime = self.config_dict['genref_exec']['@threeprime']
		for charbase in threeprime:
			if charbase not in valid_base:
				invalid_flags['ThreePrime'] = 'User specified threeprime flag is not \'ATGC\''

		##
		## CLCMap_exec settings (CLC Alignment)
		cpu_threads = self.config_dict['clcmap_exec']['@cpu_threads']
		if not cpu_threads.isdigit():
			invalid_flags['CPUThreads'] = 'User specified cpu_threads flag is not an integer!'
		elif int(cpu_threads) == 0 or int(cpu_threads) > cpu_count():
			invalid_flags['CPUThreads'] = 'User specified cpu_threads flag outside of range (1-' + str(cpu_count()) + ')!'

		paired_read = self.config_dict['clcmap_exec']['@paired_read']
		if not (paired_read == 'True' or paired_read == 'False'):
			invalid_flags['PairedRead'] = 'User specified paired_read flag is not valid! (valid: True, False)'

		mismatch_cost = self.config_dict['clcmap_exec']['@mismatch_cost']
		if not mismatch_cost.isdigit():
			invalid_flags['MismatchCost'] = 'User specified mismatch_cost flag is not an integer!'
		elif int(mismatch_cost) not in range(1,4):
			invalid_flags['MismatchCost'] = 'User specified mistmatch_cost flag is out of bounds! (valid: 1-3)'

		gap_cost = self.config_dict['clcmap_exec']['@gap_cost']
		if not gap_cost.isdigit():
			invalid_flags['GapCost'] = 'User specified gap_cost flag is not an integer!'
		elif int(gap_cost) not in range(1,4):
			invalid_flags['GapCost'] = 'User specified gap_cost flag is out of bounds! (valid: 1-3)'

		delete_cost = self.config_dict['clcmap_exec']['@delete_cost']
		if not delete_cost.isdigit():
			invalid_flags['DeleteCost'] = 'User specified delete_cost flag is not an integer!'
		elif int(delete_cost) not in range(1,4):
			invalid_flags['DeleteCost'] = 'User specified delete_cost flag is out of bounds! (valid: 1-3)'

		max_align = self.config_dict['clcmap_exec']['@max_align']
		if not max_align.isdigit():
			invalid_flags['MaxAlign'] = 'User specified max_align flag is not an integer!'

		similar_score = self.config_dict['clcmap_exec']['@similar_score']
		if not float(similar_score):
			invalid_flags['SimilarityScore'] = 'User specified similiar_score flag is not a decimal integer (float)!'

		length_frac = self.config_dict['clcmap_exec']['@length_frac']
		if not float(length_frac):
			invalid_flags['LengthFraction'] = 'User specified length_frac flag is not a decimal integer (float)!'

		if len(invalid_flags) > 0:
			errors = ''
			for key in invalid_flags:
				errors += invalid_flags[key] + '\n'
			errors = errors.rstrip()
			errstr = str(errors)
			log.error(errstr)

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
			sys.exit(2)
		for samfile in parsed_arguments.input:
			if not check_input_files('.sam',samfile):
				sys.exit(2)

	if parsed_arguments.batch:
		if not filesystem_exists_check(parsed_arguments.batch[0]):
			sys.exit(2)
		for samfile in parsed_arguments.batch:
			if not check_input_files('.sam',samfile):
				sys.exit(2)

	if parsed_arguments.config:
		if not filesystem_exists_check(parsed_arguments.config[0]):
			sys.exit(2)
		for xmlfile in parsed_arguments.config:
			if not check_input_files('.xml',xmlfile):
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