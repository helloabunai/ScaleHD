#/usr/bin/python
__version__ = 0.321
__author__ = 'alastair.maxwell@glasgow.ac.uk'

##
## Imports
import string
import os
import errno
import shutil
import sys
import glob
import datetime
import subprocess
import logging as log
import numpy as np
import csv
import StringIO
import PyPDF2
from sklearn import preprocessing
from collections import defaultdict
from xml.etree import cElementTree
from lxml import etree
from reportlab.pdfgen import canvas

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
	def __init__(self):
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
			log.error('{}{}{}{}'.format(Colour.red, 'shd__ ', Colour.end, 'XML Config: Specified data directory could not be found.'))
			trigger = True
		for fqfile in glob.glob(os.path.join(data_directory, '*')):
			if not (fqfile.endswith('.fq') or fqfile.endswith('.fastq') or fqfile.endswith('.fq.gz') or fqfile.endswith('.fastq.gz')):
				log.error('{}{}{}{}'.format(Colour.red, 'shd__ ', Colour.end, 'XML Config: Non FastQ/GZ data detected in specified input directory.'))
				trigger = True
		forward_reference = self.config_dict['@forward_reference']
		if not os.path.isfile(forward_reference):
			log.error('{}{}{}{}'.format(Colour.red, 'shd__ ', Colour.end, 'XML Config: Specified forward reference file could not be found.'))
			trigger = True
		if not (forward_reference.endswith('.fa') or forward_reference.endswith('.fasta')):
			log.error('{}{}{}{}'.format(Colour.red, 'shd__ ', Colour.end, 'XML Config: Specified forward reference file is not a fa/fas file.'))
			trigger = True
		reverse_reference = self.config_dict['@reverse_reference']
		if not os.path.isfile(reverse_reference):
			log.error('{}{}{}{}'.format(Colour.red, 'shd__ ', Colour.end, 'XML Config: Specified reverse reference file could not be found.'))
			trigger = True
		if not (reverse_reference.endswith('fa') or reverse_reference.endswith('.fasta')):
			log.error('{}{}{}{}'.format(Colour.red, 'shd__ ', Colour.end, 'XML Config: Specified reverse reference file is not a fa/fas file.'))
			trigger = True
		if forward_reference.split('/')[-1] == reverse_reference.split('/')[-1]:
			log.error('{}{}{}{}'.format(Colour.red, 'shd__ ', Colour.end, 'XML Config: FW and RV references have identical filenames. Will create indexing issue.'))
			trigger = True

		##
		## Instance flag settings
		demultiplexing_flag = self.config_dict['instance_flags']['@demultiplex']
		if not (demultiplexing_flag == 'True' or demultiplexing_flag == 'False'):
			log.error('{}{}{}{}'.format(Colour.red, 'shd__ ', Colour.end, 'XML Config: Demultiplexing flag is not set to True/False.'))
			trigger = True
		sequence_qc_flag = self.config_dict['instance_flags']['@quality_control']
		if not (sequence_qc_flag == 'True' or sequence_qc_flag == 'False'):
			log.error('{}{}{}{}'.format(Colour.red, 'shd__ ', Colour.end, 'XML Config: Sequence Quality control flag is not set to True/False.'))
			trigger = True
		alignment_flag = self.config_dict['instance_flags']['@sequence_alignment']
		if not (alignment_flag == 'True' or alignment_flag == 'False'):
			log.error('{}{}{}{}'.format(Colour.red, 'shd__ ', Colour.end, 'XML Config: Sequence Alignment flag is not set to True/False.'))
			trigger = True
		atypical_flag = self.config_dict['instance_flags']['@atypical_realignment']
		if not (atypical_flag == 'True' or atypical_flag == 'False'):
			log.error('{}{}{}{}'.format(Colour.red, 'shd__ ', Colour.end, 'XML Config: Atypical Realignment flag is not True/False.'))
			trigger = True
		genotype_flag = self.config_dict['instance_flags']['@genotype_prediction']
		if not (genotype_flag == 'True' or genotype_flag == 'False'):
			log.error('{}{}{}{}'.format(Colour.red, 'shd__ ', Colour.end, 'XML Config: Genotype Prediction control flag is not True/False.'))
			trigger = True
		snpcall_flag = self.config_dict['instance_flags']['@snp_calling']
		if not (snpcall_flag == 'True' or snpcall_flag == 'False'):
			log.error('{}{}{}{}'.format(Colour.red, 'shd__ ', Colour.end, 'XML Config: SNP Calling flag is not True/False.'))
			trigger = True

		##
		## Demultiplexing flag settings
		trim_adapter_base = ['A', 'G', 'C', 'T']
		if demultiplexing_flag == 'True':
			forward_adapter = self.config_dict['demultiplex_flags']['@forward_adapter']
			for charbase in forward_adapter:
				if charbase not in trim_adapter_base:
					log.error('{}{}{}{}'.format(Colour.red, 'shd__ ', Colour.end, 'XML Config: Invalid character detected in forward_adapter demultiplexing flag.'))
					trigger = True
			forward_position = self.config_dict['demultiplex_flags']['@forward_position']
			if forward_position not in ['5P', '3P', 'AP']:
				log.error('{}{}{}{}'.format(Colour.red, 'shd__ ', Colour.end, 'XML Config: Given demultiplexing forward adapter position invalid! [5P, 3P, AP]'))
				trigger = True

			reverse_adapter = self.config_dict['demultiplex_flags']['@reverse_adapter']
			for charbase in reverse_adapter:
				if charbase not in trim_adapter_base:
					log.error('{}{}{}{}'.format(Colour.red, 'shd__ ', Colour.end, 'XML Config: Invalid character detected in reverse_adapter demultiplexing flag.'))
					trigger = True
			reverse_position = self.config_dict['demultiplex_flags']['@reverse_position']
			if reverse_position not in ['5P', '3P', 'AP']:
				log.error('{}{}{}{}'.format(Colour.red, 'shd__ ', Colour.end, 'XML Config: Given demultiplexing reverse adapter position invalid! [5P, 3P, AP]'))
				trigger = True

			error_rate = self.config_dict['demultiplex_flags']['@error_rate']
			if not error_rate.isdigit():
				log.error('{}{}{}{}'.format(Colour.red, 'shd__ ', Colour.end, 'XML Config: Specified error_rate is not a valid integer.'))
				trigger = True
			minimum_overlap = self.config_dict['demultiplex_flags']['@min_overlap']
			if not minimum_overlap.isdigit():
				log.error('{}{}{}{}'.format(Colour.red, 'shd__ ', Colour.end, 'XML Config: Specified min_overlap is not a valid integer.'))
				trigger = True
			minimum_length = self.config_dict['demultiplex_flags']['@min_length']
			if not minimum_length == '':
				if not minimum_length.isdigit():
					log.error('{}{}{}{}'.format(Colour.red, 'shd__ ', Colour.end, 'XML Config: Specified min_length is not a valid integer.'))
					trigger = True
			maximum_length = self.config_dict['demultiplex_flags']['@max_length']
			if not maximum_length == '':
				if not maximum_length.isdigit():
					log.error('{}{}{}{}'.format(Colour.red, 'shd__ ', Colour.end, 'XML Config: Specified max_length is not a valid integer.'))
					trigger = True

		##
		## Trimming flag settings
		if sequence_qc_flag == 'True':
			trimming_type = self.config_dict['trim_flags']['@trim_type']
			if not (trimming_type == 'Quality' or trimming_type	== 'Adapter' or trimming_type == 'Both'):
				log.error('{}{}{}{}'.format(Colour.red, 'shd__  ', Colour.end, 'XML Config: Trimming type is not Quality/Adapter/Both.'))
				trigger = True
			quality_threshold = self.config_dict['trim_flags']['@quality_threshold']
			if not quality_threshold.isdigit():
				log.error('{}{}{}{}'.format(Colour.red, 'shd__ ', Colour.end, 'XML Config: Specified quality threshold integer is invalid.'))
				trigger = True
			elif not int(quality_threshold) in range(0,39):
				log.error('{}{}{}{}'.format(Colour.red, 'shd__ ', Colour.end, 'XML Config: Specified quality threshold integer out of range (0-38).'))
				trigger = True
			trim_adapters = ['-a','-g','-a$','-g^','-b']
			adapter_flag = self.config_dict['trim_flags']['@adapter_flag']
			if not (adapter_flag in trim_adapters):
				log.error('{}{}{}{}'.format(Colour.red, 'shd__ ', Colour.end, 'XML Config: Specified trimming adapter not valid selection.'))
				trigger = True
			forward_adapter = self.config_dict['trim_flags']['@forward_adapter']
			for charbase in forward_adapter:
				if charbase not in trim_adapter_base:
					log.error('{}{}{}{}'.format(Colour.red, 'shd__ ', Colour.end, 'XML Config: Invalid character detected in FW adapter sequence.'))
					trigger = True
			reverse_adapter = self.config_dict['trim_flags']['@reverse_adapter']
			for charbase in reverse_adapter:
				if charbase not in trim_adapter_base:
					log.error('{}{}{}{}'.format(Colour.red, 'shd__ ', Colour.end, 'XML Config: Invalid character detected in RV adapter sequence.'))
					trigger = True
			error_tolerance = self.config_dict['trim_flags']['@error_tolerance']
			if not isinstance(float(error_tolerance), float):
				log.error('{}{}{}{}'.format(Colour.red, 'shd__ ', Colour.end, 'XML Config: Specified error tolerance is not a valid float.'))
				trigger = True
			if not float(error_tolerance) in np.arange(0,1.1,0.01):
				log.error('{}{}{}{}'.format(Colour.red, 'shd__ ', Colour.end, 'XML Config: Specified error tolerance is not 0.0 < x < 1.0.'))
				trigger = True

		##
		## Alignment flag settings
		if alignment_flag == 'True':
			min_seed_length = self.config_dict['alignment_flags']['@min_seed_length']
			if not min_seed_length.isdigit():
				log.error('{}{}{}{}'.format(Colour.red, 'shd__ ', Colour.end, 'XML Config: Specified min_seed_length integer is invalid.'))
				trigger=True

			band_width = self.config_dict['alignment_flags']['@band_width']
			if not band_width.isdigit():
				log.error('{}{}{}{}'.format(Colour.red, 'shd__ ', Colour.end, 'XML Config: Specified band_width integer is invalid.'))
				trigger=True

			seed_length_extension = self.config_dict['alignment_flags']['@seed_length_extension']
			if not isinstance(float(seed_length_extension), float):
				log.error('{}{}{}{}'.format(Colour.red, 'shd__ ', Colour.end, 'XML Config: Specified seed_length_extension float is invalid.'))
				trigger=True

			skip_seed_with_occurrence = self.config_dict['alignment_flags']['@skip_seed_with_occurrence']
			if not skip_seed_with_occurrence.isdigit():
				log.error('{}{}{}{}'.format(Colour.red, 'shd__ ', Colour.end, 'XML Config: Specified skip_seed_with_occurrence integer is invalid.'))
				trigger=True

			chain_drop = self.config_dict['alignment_flags']['@chain_drop']
			if not isinstance(float(chain_drop), float):
				log.error('{}{}{}{}'.format(Colour.red, 'shd__ ', Colour.end, 'XML Config: Specified chain_drop float is invalid.'))
				trigger=True

			seeded_chain_drop = self.config_dict['alignment_flags']['@seeded_chain_drop']
			if not seeded_chain_drop.isdigit():
				log.error('{}{}{}{}'.format(Colour.red, 'shd__ ', Colour.end, 'XML Config: Specified seeded_chain_drop integer is invalid.'))
				trigger=True

			seq_match_score = self.config_dict['alignment_flags']['@seq_match_score']
			if not seq_match_score.isdigit():
				log.error('{}{}{}{}'.format(Colour.red, 'shd__ ', Colour.end, 'XML Config: Specified seq_match_score integer is invalid.'))
				trigger=True

			mismatch_penalty = self.config_dict['alignment_flags']['@mismatch_penalty']
			if not mismatch_penalty.isdigit():
				log.error('{}{}{}{}'.format(Colour.red, 'shd__ ', Colour.end, 'XML Config: Specified mismatch_penalty integer is invalid.'))
				trigger=True

			indel_penalty_raw = self.config_dict['alignment_flags']['@indel_penalty']
			indel_penalty = indel_penalty_raw.split(',')
			for individual_indelpen in indel_penalty:
				if not individual_indelpen.isdigit():
					log.error('{}{}{}{}'.format(Colour.red, 'shd__ ', Colour.end, 'XML Config: Specified indel_penalty integer(s) is(are) invalid.'))
					trigger=True

			gap_extend_penalty_raw = self.config_dict['alignment_flags']['@gap_extend_penalty']
			gap_extend_penalty = gap_extend_penalty_raw.split(',')
			for individual_gaextend in gap_extend_penalty:
				if not individual_gaextend.isdigit():
					log.error('{}{}{}{}'.format(Colour.red, 'shd__ ', Colour.end, 'XML Config: Specified gap_extend_penalty integer(s) is(are) invalid.'))
					trigger=True

			prime_clipping_penalty_raw = self.config_dict['alignment_flags']['@prime_clipping_penalty']
			prime_clipping_penalty = prime_clipping_penalty_raw.split(',')
			for individual_prclip in prime_clipping_penalty:
				if not individual_prclip.isdigit():
					log.error('{}{}{}{}'.format(Colour.red, 'shd__ ', Colour.end, 'XML Config: Specified prime_clipping_penalty integer(s) is(are) invalid.'))
					trigger=True

			unpaired_pairing_penalty = self.config_dict['alignment_flags']['@unpaired_pairing_penalty']
			if not unpaired_pairing_penalty.isdigit():
				log.error('{}{}{}{}'.format(Colour.red, 'shd__ ', Colour.end, 'XML Config: Specified unpaired_pairing_penalty integer is invalid.'))
				trigger=True
		##
		## Genotype prediction flag settings
		if genotype_flag == 'True':
			snp_observation_pcnt = self.config_dict['prediction_flags']['@snp_observation_threshold']
			if not snp_observation_pcnt.isdigit():
				if not int(snp_observation_pcnt) in range(1,5):
					log.error('{}{}{}{}'.format(Colour.red, 'shd__ ', Colour.end, 'XML Config: SNP Observation value invalid! Please use 1-10.'))
					trigger = True

		variant_algorithm = self.config_dict['prediction_flags']['@algorithm_utilisation']
		if not variant_algorithm in ['freebayes', 'gatk']:
			log.error('{}{}{}{}'.format(Colour.red, 'shd__ ', Colour.end, 'XML Config: Specified variant_algorithm value is invalid. [freebayes/gatk]'))
			trigger = True

		quality_cutoff = self.config_dict['prediction_flags']['@quality_cutoff']
		if not quality_cutoff.isdigit():
			log.error('{}{}{}{}'.format(Colour.red, 'shd__ ', Colour.end, 'XML Config: SNP Quality Cutoff value is not an integer.'))
			trigger = True

		if trigger:
			log.error('{}{}{}{}'.format(Colour.red, 'shd__ ', Colour.end, 'XML Config: Failure, exiting.'))
			sys.exit(2)
		else:
			log.info('{}{}{}{}'.format(Colour.green, 'shd__ ', Colour.end, 'XML Config: Parsing parameters successful!'))


class DataClump(dict):
	"""Container object for datasets: dictionary-like object that
       exposes its keys as attributes."""

	def __init__(self, **kwargs):
		dict.__init__(self, kwargs)
		self.__dict__ = self


class DataLoader:

	def __init__(self, database, descriptor):

		self.database = database
		self.descriptor = descriptor

	def load_model(self):

			## Loads description file for respective data set
			modeldescr_name = self.descriptor
			with open(modeldescr_name) as f:
				descr_text = f.read()

			## Loads data set from csv, into objects in preparation for bunch()
			data_file_name = self.database
			with open(data_file_name) as f:
				data_file = csv.reader(f)
				temp = next(data_file)
				n_samples = int(temp[0])
				n_features = int(temp[1])
				data = np.empty((n_samples, n_features))
				temp = next(data_file)
				feature_names = np.array(temp)

				labels = []
				for i, d in enumerate(data_file):
					data[i] = d[:-1]
					label = d[-1]
					labels.append(label)

				le = preprocessing.LabelEncoder()
				le.fit(labels)
				hash_int_labels = le.transform(labels)

			return DataClump(DATA=data,
						 	 TARGET=hash_int_labels,
						 	 FTRNAME=feature_names[:-1],
							 DESCR=descr_text,
							 ENCDR=le)


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

	trigger = False

	##
	## Jobname prefix validity check
	if parsed_arguments.jobname:
		for character in parsed_arguments.jobname:
			if character is ' ' or character is '/':
				log.error('{}{}{}{}'.format(Colour.red,'shd__ ',Colour.end,'Specified Job Name has invalid characters: "', character, '"'))
				trigger = True

	##
	## Config mode check
	if parsed_arguments.config:
		if not filesystem_exists_check(parsed_arguments.config[0]):
			log.error('{}{}{}{}'.format(Colour.red, 'shd__ ', Colour.end, 'Specified config file could not be found.'))
			trigger = True

		for xmlfile in parsed_arguments.config:
			if not check_input_files('.xml',xmlfile):
				log.error('{}{}{}{}'.format(Colour.red, 'shd__ ', Colour.end, 'Specified config file is not an XML file.'))
				trigger = True
	return trigger

def extract_data(input_data_directory):

	target_files = glob.glob(os.path.join(input_data_directory, '*'))
	for extract_target in target_files:
		if extract_target.lower().endswith(('.fq.gz', '.fastq.gz')):
			log.info('{}{}{}{}'.format(Colour.bold, 'shd__ ', Colour.end, 'Detected compressed input data. Extracting!'))
			break

	for extract_target in target_files:
		unzipd = subprocess.Popen(['gzip', '-q', '-f', '-d', extract_target], stderr=subprocess.PIPE)
		unzipd.wait()

	return True

def sequence_pairings(data_path, instance_rundir):

	##
	## Get input files from data path
	## Sort so that ordering isn't screwy on linux
	input_files = glob.glob(os.path.join(data_path, '*'))
	sorted_input = sorted(input_files)
	sequence_pairs = []

	file_count = len(sorted_input)
	if not file_count % 2 == 0:
		log.error('{}{}{}{}'.format(Colour.red, 'shd__ ', Colour.end, 'I/O: Non-even number of input files specified. Cannot continue without pairing!'))
		sys.exit(2)

	##
	## Optimise so code isn't recycled
	for i in range(0, len(sorted_input), 2):
		file_pair = {}
		forward_data = sorted_input[i]
		reverse_data = sorted_input[i+1]

		##
		## Check forward ends with R1
		forward_data_name = sorted_input[i].split('/')[-1].split('.')[0]
		if not forward_data_name.endswith('_R1'):
			log.error('{}{}{}{}{}'.format(Colour.red,'shd__ ',Colour.end,'I/O: Forward input file does not end in _R1. ', forward_data))
			sys.exit(2)

		##
		## Check reverse ends with R2
		reverse_data_name = sorted_input[i+1].split('/')[-1].split('.')[0]
		if not reverse_data_name.endswith('_R2'):
			log.error('{}{}{}{}{}'.format(Colour.red,'shd__ ',Colour.end,'I/O: Reverse input file does not end in _R2. ', reverse_data))
			sys.exit(2)

		##
		## Make Stage outputs for use in everywhere else in pipeline
		sample_root = '_'.join(forward_data_name.split('_')[:-1])
		instance_path = os.path.join(instance_rundir)
		seq_qc_path = os.path.join(instance_rundir, sample_root, 'SeqQC')
		align_path = os.path.join(instance_rundir, sample_root, 'Align')
		predict_path = os.path.join(instance_rundir, sample_root, 'Predict')
		file_pair[sample_root] = [forward_data, reverse_data, instance_path, seq_qc_path, align_path, predict_path]
		sequence_pairs.append(file_pair)

	return sequence_pairs

def filesystem_exists_check(path, raise_exception=True):

	"""
	Checks to see if the path, specified by parameter path, exists. Can be either a directory or file.
	If the path exists, True is returned. If the path does not exist, and raise_exception is set to True,
	an IOError is raised - else False is returned.
	"""

	if os.path.lexists(path):
		return True
	if raise_exception:
		log.error('{}{}{}{}'.format(Colour.red,'shd__ ',Colour.end,'Specified input path could not be found.'))
	return False

def check_input_files(input_format, input_file):

	if input_file.endswith(input_format):
		return True
	return False

def initialise_libraries(instance_params):

	trigger = False

	##
	## Subfunction for recycling code
	## Calls UNIX type for checking binaries present
	## Changed from WHICH as apparently type functions over different shells/config files
	def type_func(binary):
		binary_result = []
		binary_string = 'type {}'.format(binary)
		binary_subprocess = subprocess.Popen([binary_string], shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
		binary_result = binary_subprocess.communicate()
		binary_subprocess.wait()

		if 'not found' in binary_result[0] or binary_result[1]:
			log.critical('{}{}{}{}{}'.format(Colour.red,'shd__ ',Colour.end,'Missing binary: ', binary, '!'))
			raise ScaleHDException

	##
	## To determine which binaries to check for
	## AttributeError in the situation where instance_params origin differs
	## try for -c style, except AttributeError for -b style
	try:
		quality_control = instance_params.config_dict['instance_flags']['@quality_control']
		alignment = instance_params.config_dict['instance_flags']['@sequence_alignment']
		genotyping = instance_params.config_dict['instance_flags']['@genotype_prediction']
		snp_calling = instance_params.config_dict['instance_flags']['@snp_calling']
	except AttributeError:
		quality_control = instance_params['quality_control']
		alignment = instance_params['sequence_alignment']
		genotyping = instance_params['genotype_prediction']
		snp_calling = instance_params['snp_calling']

	if quality_control == 'True':
		try:type_func('java')
		except ScaleHDException: trigger=True
		try:type_func('fastqc')
		except ScaleHDException: trigger=True
		try:type_func('cutadapt')
		except ScaleHDException: trigger=True
	if alignment == 'True':
		try:type_func('seqtk')
		except ScaleHDException: trigger=True
		try:type_func('bwa')
		except ScaleHDException: trigger=True
		try:type_func('samtools')
		except ScaleHDException: trigger=True
		try:type_func('generatr')
		except ScaleHDException: trigger=True
	if genotyping == 'True':
		try:type_func('samtools')
		except ScaleHDException: trigger=True
		try:type_func('generatr')
		except ScaleHDException: trigger=True
	if snp_calling == 'True':
		try: type_func('picard')
		except ScaleHDException: trigger=True
		try: type_func('freebayes')
		except ScaleHDException: trigger=True
		try: type_func('gatk')
		except ScaleHDException: trigger = True

	return trigger

def sanitise_outputs(jobname, output_argument):

	run_dir = ''
	output_root = output_argument[0]
	if jobname:
		target_output = os.path.join(output_root, jobname)
		if not os.path.exists(target_output):
			log.info('{}{}{}{}{}'.format(Colour.bold, 'shd__ ', Colour.end, 'Creating Output with prefix: ', jobname))
			run_dir = os.path.join(output_root, jobname)
			mkdir_p(run_dir)
		else:
			purge_choice = ''
			while True:
				purge_choice = raw_input('{}{}{}{}'.format(Colour.bold, 'shd__ ', Colour.end, 'Job folder already exists. Delete existing folder? Y/N: '))
				if not (purge_choice.lower() == 'y') and not (purge_choice.lower() == 'n'):
					log.info('{}{}{}{}'.format(Colour.red, 'shd__ ', Colour.end, 'Invalid input. Please input Y or N.'))
					continue
				else:
					break

			if purge_choice.lower() == 'y':
				log.info('{}{}{}{}{}'.format(Colour.bold, 'shd__ ', Colour.end, 'Clearing pre-existing Jobname Prefix: ', jobname))
				run_dir = os.path.join(output_root, jobname)
				if os.path.exists(run_dir):
					shutil.rmtree(run_dir, ignore_errors=True)
				mkdir_p(run_dir)
			else:
				raise Exception('User chose not to delete pre-existing Job folder. Cannot write output.')

	else:
		## Ensures root output is a real directory
		## Generates folder name based on date (for run ident)
		date = datetime.date.today().strftime('%d-%m-%Y')
		walltime = datetime.datetime.now().strftime('%H%M%S')
		today = date + '-' + walltime

		## If the user specified root doesn't exist, make it
		## Then make the run directory for datetime
		if not os.path.exists(output_root):
			log.info('{}{}{}{}'.format(Colour.bold, 'shd__ ', Colour.end, 'Creating output root... '))
			mkdir_p(output_root)
		run_dir = os.path.join(output_root, 'ScaleHDRun_'+today)
		log.info('{}{}{}{}'.format(Colour.bold, 'shd__ ', Colour.end, 'Creating instance run directory.. '))
		mkdir_p(run_dir)

		## Inform user it's all gonna be okaaaayyyy
		log.info('{}{}{}{}'.format(Colour.green, 'shd__ ', Colour.end, 'Output directories OK!'))

	return run_dir

def replace_fqfile(mutate_list, target_fqfile, altered_path):

	if target_fqfile in mutate_list:
		loc = mutate_list.index(target_fqfile)
		mutate_list[loc] = altered_path
	return mutate_list

def scrape_summary_data(stage, input_report_file):
	##
	## If the argument input_report_file is from trimming..
	if stage == 'trim':
		with open(input_report_file, 'r') as trpf:
			trim_lines = trpf.readlines()
			##
			## Determine buffer size to slice from above array
			scraping_buffer = 8
			if '-q' in trim_lines[1]:
				scraping_buffer += 1
			##
			## Get Anchor
			summary_start = 0
			for i in range(0, len(trim_lines)):
				if '== Summary ==' in trim_lines[i]:
					summary_start = i
			##
			## Slice and close
			summary_data = trim_lines[summary_start:summary_start + scraping_buffer]
			trpf.close()
		return summary_data[2:]

	##
	## If the argument input_report_file is from alignment..
	if stage == 'align':
		with open(input_report_file, 'r') as alnrpf:
			align_lines = alnrpf.readlines()
			alnrpf.close()
		##
		## No ranges required, only skip first line
		return align_lines[1:]

	##
	## No need to tidy up report for genotyping
	## since we already have the data from our own objects
	if stage == 'gtype':
		pass

def collate_peaks(predict_path, sample_prefix):

	##TODO docstring

	##
	## Generate single page PDF of current sample name (for ctrl+f)
	packet = StringIO.StringIO()
	can = canvas.Canvas(packet, pagesize=(350,150))
	can.drawCentredString(175, 75, sample_prefix)
	can.save()

	packet.seek(0)
	header_page = PyPDF2.PdfFileReader(packet).getPage(0)
	header_output = PyPDF2.PdfFileWriter()
	header_output.addPage(header_page)
	test = os.path.join(predict_path, 'Header.pdf')
	header_output.write(file(test, 'wb'))

	##
	## Scan prediction folder for the passed sample instance
	sample_merge = os.path.join(predict_path, 'SampleResults.pdf')
	merge_object = PyPDF2.PdfFileMerger()
	for peak_graph in glob.glob(os.path.join(predict_path, '*')):
		if 'PeakDetection' in peak_graph:
			try:
				curr_file = file(peak_graph, 'rb')
				merge_object.append(PyPDF2.PdfFileReader(curr_file))
				curr_file.close()
			except PyPDF2.utils.PdfReadError:
				pass
		if 'Header' in peak_graph:
			try:
				curr_file = file(peak_graph, 'rb')
				merge_object.merge(0, PyPDF2.PdfFileReader(curr_file))
				curr_file.close()
			except PyPDF2.utils.PdfReadError:
				pass

	merge_object.write(sample_merge)
	merge_object.close()
	os.remove(os.path.join(predict_path, 'Header.pdf'))
	return sample_merge

def generate_atypical_xml(label, allele_object, index_path, direction):

	"""
	:param allele_object:
	:param index_path:
	:return:
	"""
	##TODO docstring

	atypical_path = os.path.join(index_path, '{}{}_{}.xml'.format(direction, label, allele_object.get_reflabel()))
	fp_flank = 'GCGACCCTGGAAAAGCTGATGAAGGCCTTCGAGTCCCTCAAGTCCTTC'
	cagstart = ''; cagend = ''
	intv = allele_object.get_intervening()
	ccgstart = ''; ccgend = ''
	ccglen = allele_object.get_ccg()
	cctlen = allele_object.get_cct()
	tp_flank = 'CAGCTTCCTCAGCCGCCGCCGCAGGCACAGCCGCTGCT'

	if direction == 'fw':
		cagstart = '1'; cagend = '200'
		ccgstart = '1'; ccgend = '20'
	if direction == 'rv':
		cagstart = '100'; cagend = '100'
		ccgstart = '1'; ccgend = '20'

	##
	## Create XML
	data_root = etree.Element('data')
	loci_root = etree.Element('loci', label=allele_object.get_reflabel()); data_root.append(loci_root)

	##
	## Loci Nodes
	fp_input = etree.Element('input', type='fiveprime', flank=fp_flank)
	cag_region = etree.Element('input', type='repeat_region', order='1', unit='CAG', start=cagstart, end=cagend)
	intervening = etree.Element('input', type='intervening', sequence=intv, prior='1')
	ccg_region = etree.Element('input', type='repeat_region', order='2', unit='CCG', start=ccgstart, end=ccgend)
	cct_region = etree.Element('input', type='repeat_region', order='3', unit='CCT', start=str(cctlen), end=str(cctlen))
	tp_input = etree.Element('input', type='threeprime', flank=tp_flank)

	for node in [fp_input, cag_region, intervening, ccg_region, cct_region, tp_input]:
		loci_root.append(node)

	s = etree.tostring(data_root, pretty_print=True)
	with open(atypical_path, 'w') as xmlfi:
		xmlfi.write(s)
		xmlfi.close()

	return atypical_path

def generate_reference(input_xml, index_path, ref_indexes, direction):

	##TODO docstring

	label = input_xml.split('/')[-1].split('.')[0]
	target_output = os.path.join(index_path, label + '.fa')
	temp_output = os.path.join(index_path, label + '_concat.fa')
	gen_process = subprocess.Popen(['generatr', '-i', input_xml, '-o', target_output], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
	gen_process.wait()

	##
	## Join typical and atypical reference into one file
	if direction == 'fw':
		toutfi = open(temp_output, 'w')
		cat_process = subprocess.Popen(['cat', target_output, ref_indexes[0]], stdout=toutfi, stderr=subprocess.PIPE)
		cat_process.wait()
		toutfi.close()
		target_output = temp_output

	return target_output

def seek_target(input_list, target):

	for i in range(0, len(input_list)):
		if target in input_list[i]:
			return i

def sanitise_trimming_output(input_object, input_list):

	if type(input_object) is int:
		cleanse_target = input_list[input_object].split(':')[1].lstrip().rstrip()
		return cleanse_target
	else:
		return '*'

def sanitise_alignment_output(input_object, input_list, stage):

	if type(input_object) is int:

		if stage == 3:
			cleanse_target = input_list[input_object].lstrip().rstrip().split(' ')[0:1]
			return ''.join(cleanse_target)
		else:
			cleanse_target = input_list[input_object].lstrip().rstrip().split(' ')[0:2]
			return ' '.join(cleanse_target)
	else:
		return '*'

def mkdir_p(path):
	try: os.makedirs(path)
	except OSError as exc:
		if exc.errno == errno.EEXIST and os.path.isdir(path): pass
		else: raise
