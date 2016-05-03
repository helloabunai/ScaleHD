#/usr/bin/python
__version__ = 0.01
__author__ = 'alastair.maxwell@glasgow.ac.uk'

##
## Imports
import string
import os
import sys
import glob
import datetime
import subprocess
import logging as log
import numpy as np
import csv
from sklearn import preprocessing
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
			log.error('{}{}{}{}'.format(Colour.red, 'shd__', Colour.end, 'XML Config: Specified data directory could not be found.'))
			trigger = True
		for fqfile in glob.glob(os.path.join(data_directory, '*')):
			if not (fqfile.endswith('.fq') or fqfile.endswith('.fastq') or fqfile.endswith('.fq.gz') or fqfile.endswith('.fastq.gz')):
				log.error('{}{}{}{}'.format(Colour.red, 'shd__ ', Colour.end, 'XML Config: Non FastQ/GZ data detected in specified input directory.'))
				trigger = True
		forward_reference = self.config_dict['@forward_reference']
		if not os.path.isfile(forward_reference):
			log.error('{}{}{}{}'.format(Colour.red, 'shd__ ', Colour.end, 'XML Config: Specified forward reference file could not be found.'))
			trigger = True
		if not (forward_reference.endswith('.fa') or forward_reference.endswith('.fas')):
			log.error('{}{}{}{}'.format(Colour.red, 'shd__ ', Colour.end, 'XML Config: Specified forward reference file is not a fa/fas file.'))
			trigger = True
		reverse_reference = self.config_dict['@reverse_reference']
		if not os.path.isfile(reverse_reference):
			log.error('{}{}{}{}'.format(Colour.red, 'shd__ ', Colour.end, 'XML Config: Specified reverse reference file could not be found.'))
			trigger = True
		if not (reverse_reference.endswith('fa') or reverse_reference.endswith('.fas')):
			log.error('{}{}{}{}'.format(Colour.red, 'shd__ ', Colour.end, 'XML Config: Specified reverse reference file is not a fa/fas file.'))
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
		## SeqQC check
		demultiplex_flag = self.config_dict['dmpx_flags']['@demultiplex_data']
		trimming_flag = self.config_dict['trim_flags']['@trim_data']

		##
		## Demultiplexing flag settings
		if sequence_qc_flag == 'True':
			if not (demultiplex_flag == 'True' or demultiplex_flag == 'False'):
				log.error('{}{}{}{}'.format(Colour.red, 'shd__ ', Colour.end, 'XML Config: Demultiplexing flag is not True/False.'))
				trigger = True
			if demultiplex_flag == 'True':
				barcode_bases = ['A','G','C','T','U','N']
				barcode_sequence = self.config_dict['dmpx_flags']['@barcode']
				for charbase in barcode_sequence:
					if charbase not in barcode_bases:
						log.error('{}{}{}{}'.format(Colour.red, 'shd__ ', Colour.end, 'XML Config: Invalid character detected in barcode sequence.'))
						trigger = True
				demultiplex_mismatch = self.config_dict['dmpx_flags']['@max_mismatch']
				if not demultiplex_mismatch.isdigit():
					log.error('{}{}{}{}'.format(Colour.red, 'shd__ ', Colour.end, 'XML Config: Specified barcode mismatch integer is invalid.'))
					trigger = True
				retain_unmatched = self.config_dict['dmpx_flags']['@retain_unmatched']
				if not (retain_unmatched == 'True' or retain_unmatched == 'False'):
					log.error('{}{}{}{}'.format(Colour.red, 'shd__ ', Colour.end, 'XML Config: Specified retain unmatched reads value is not True/False.'))
					trigger = True

		##
		## Trimming flag settings
		if sequence_qc_flag == 'True':
			if not (trimming_flag == 'True' or trimming_flag == 'False'):
				log.error('{}{}{}{}'.format(Colour.red, 'shd__ ', Colour.end, 'XML Config: Trimming flag is not True/False.'))
				trigger = True
			if trimming_flag == 'True':
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
				trim_adapter_base = ['A','G','C','T']
				adapter_sequence = self.config_dict['trim_flags']['@adapter']
				for charbase in adapter_sequence:
					if charbase not in trim_adapter_base:
						log.error('{}{}{}{}'.format(Colour.red, 'shd__ ', Colour.end, 'XML Config: Invalid character detected in adapter sequence.'))
						trigger = True
				error_tolerance = self.config_dict['trim_flags']['@error_tolerance']
				if not isinstance(float(error_tolerance), float):
					log.error('{}{}{}{}'.format(Colour.red, 'shd__ ', Colour.end, 'XML Config: Specified error tolerance is not a valid float.'))
					trigger = True
				if not float(error_tolerance) in np.arange(0,1.1,0.1):
					log.error('{}{}{}{}'.format(Colour.red, 'shd__ ', Colour.end, 'XML Config: Specified error tolerance is not 0.0 < x < 1.0.'))
					trigger = True

		##
		## SeqQC stage check
		if sequence_qc_flag == 'True':
			if not (demultiplex_flag == 'True' or trimming_flag == 'True'):
				log.error('{}{}{}{}'.format(Colour.red, 'shd__ ', Colour.end, 'XML Config: Sequence Quality Control flag is true, but neither Demultiplexing or Trimming are.'))
				trigger = True

		##
		## Alignment flag settings
		if alignment_flag == 'True':
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
			if not int(align_mismatch) in range(0,2):
				log.error('{}{}{}{}'.format(Colour.red, 'shd__ ', Colour.end, 'XML Config: Align Mismatch integer is out of range (0,1).'))
				trigger = True
			substring_length = self.config_dict['alignment_flags']['@substr_length']
			if not substring_length.isdigit():
				log.error('{}{}{}{}'.format(Colour.red, 'shd__ ', Colour.end, 'XML Config: Specified substring length integer is invalid.'))
				trigger=True
			if not int(substring_length) in range(3, 33):
				log.error('{}{}{}{}'.format(Colour.red, 'shd__ ', Colour.end, 'XML Config: Substring length integer out of range (3,32).'))
				trigger=True
			substring_interval_start = self.config_dict['alignment_flags']['@substr_interval_start']
			if not int(substring_interval_start) in range(0,2):
				log.error('{}{}{}{}'.format(Colour.red, 'shd__ ', Colour.end, 'XML Config: Seed Substring Interval (start) integer is out of range (0,1).'))
				trigger = True
			substring_interval_end = self.config_dict['alignment_flags']['@substr_interval_end']
			if float(substring_interval_end) in np.arange(0.5,2.55,0.05):
				log.error('{}{}{}{}'.format(Colour.red, 'shd__ ', Colour.end, 'XML Config: Seed Substring Interval (end) float is out of range (0.50, 2.50).'))
				trigger = True
			read_gap_open = self.config_dict['alignment_flags']['@read_gap_open'] ##TODO range check
			if not read_gap_open.isdigit():
				log.error('{}{}{}{}'.format(Colour.red, 'shd__ ',Colour.end, 'XML Config: Read Gap Open value is not an integer.'))
				trigger=True
			read_gap_extend = self.config_dict['alignment_flags']['@read_gap_extend'] ## TODO range check
			if not read_gap_extend.isdigit():
				log.error('{}{}{}{}'.format(Colour.red, 'shd__ ',Colour.end, 'XML Config: Read Gap Extend value is not an integer.'))
				trigger=True
			ref_gap_open = self.config_dict['alignment_flags']['@ref_gap_open'] ## TODO range check
			if not ref_gap_open.isdigit():
				log.error('{}{}{}{}'.format(Colour.red,'shd__ ',Colour.end, 'XML Config: Reference gap open value is not an integer.'))
				trigger=True
			ref_gap_extend = self.config_dict['alignment_flags']['@ref_gap_extend'] ## TODO range check
			if not ref_gap_extend.isdigit():
				log.error('{}{}{}{}'.format(Colour.red,'shd__ ',Colour.end,'XML Config: Reference gap extend value is not an integer.'))
				trigger=True
			min_mismatch_penality = self.config_dict['alignment_flags']['@min_mismatch_pen'] ##TODO range check
			if not min_mismatch_penality.isdigit():
				log.error('{}{}{}{}'.format(Colour.red,'shd__ ',Colour.end,'XML Config: Minimum mismatch penalty value is not an integer.'))
				trigger=True
			max_mismatch_penalty = self.config_dict['alignment_flags']['@max_mismatch_pen'] ##TODO range check
			if not max_mismatch_penalty.isdigit():
				log.error('{}{}{}{}'.format(Colour.red,'shd__ ',Colour.end,'XML Config: Maximum mismatch penalty value is not an integer.'))
				trigger=True

			## TODO Ensure minimums are >> maximum for these todo ranges

		##
		## Genotype prediction flag settings
		if genotype_flag == 'True':
			probability_estimate = self.config_dict['prediction_flags']['@probablity_estimate']
			if not (probability_estimate == 'True' or probability_estimate == 'False'):
				log.error('{}{}{}{}'.format(Colour.red, 'shd__ ', Colour.end, 'XML Config: Probability estimate is not True/False.'))
				trigger = True
			max_iteration = self.config_dict['prediction_flags']['@max_iteration']
			try: int(max_iteration)
			except ValueError:
				log.error('{}{}{}{}'.format(Colour.red, 'shd__ ', Colour.end, 'XML Config: SVM Max iteration integer is invalid.'))
				trigger = True
			if not int(max_iteration) in range(-1, 100001):
				log.error('{}{}{}{}'.format(Colour.red, 'shd__ ', Colour.end, 'XML Config: SVM Max iteration integer is out of range (-1, 100000).'))
				trigger	= True

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


class DataLoader():

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

	if parsed_arguments.batch:
		if not filesystem_exists_check(parsed_arguments.batch[0]):
			log.error('{}{}{}{}'.format(Colour.red, 'shd__ ', Colour.end, 'Specified batch folder could not be found.')); print 'batch not found'
			trigger = True
		for samfile in glob.glob(os.path.join(parsed_arguments.batch[0], '*')):
			if not check_input_files('.sam',samfile):
				log.error('{}{}{}{}'.format(Colour.red, 'shd__ ', Colour.end, 'Specified batch folder contains non SAM files.')); print 'batch not sam'
				trigger = True

	if parsed_arguments.config:

		if not filesystem_exists_check(parsed_arguments.config[0]):
			log.error('{}{}{}{}'.format(Colour.red, 'shd__ ', Colour.end, 'Specified config file could not be found.')); print 'config not found'
			trigger = True

		for xmlfile in parsed_arguments.config:
			if not check_input_files('.xml',xmlfile):
				log.error('{}{}{}{}'.format(Colour.red, 'shd__ ', Colour.end, 'Specified config file is not an XML file.')); print 'config not xml'
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

def sequence_pairings(data_path, instance_rundir, workflow_type):

	input_files = glob.glob(os.path.join(data_path, '*'))
	sequence_pairs = []

	file_count = len(input_files)
	if not file_count % 2 == 0:
		log.error('{}{}{}{}'.format(Colour.red, 'shd__ ', Colour.end, 'I/O: Non-even number of input files specified. Cannot continue without pairing!'))
		sys.exit(2)

	##
	## Optimise so code isn't recycled
	for i in range(0, len(input_files), 2):
		file_pair = {}
		forward_data = input_files[i]
		reverse_data = input_files[i+1]

		##
		## Check forward ends with R1
		forward_data_name = input_files[i].split('/')[-1].split('.')[0]
		if not forward_data_name.endswith('R1'):
			log.error('{}{}{}{}{}'.format(Colour.red,'shd__ ',Colour.end,'I/O: Forward input file does not end in _R1. ', forward_data))
			sys.exit(2)

		##
		## Check reverse ends with R2
		reverse_data_name = input_files[i+1].split('/')[-1].split('.')[0]
		if not reverse_data_name.endswith('R2'):
			log.error('{}{}{}{}{}'.format(Colour.red,'shd__ ',Colour.end,'I/O: Reverse input file does not end in _R2. ', reverse_data))
			sys.exit(2)

		if workflow_type == 'sequence':

			##
			## Make Stage outputs for use in everywhere else in pipeline
			sample_root = '_'.join(forward_data_name.split('_')[:-1])
			seq_qc_path = os.path.join(instance_rundir, sample_root, 'SeqQC')
			align_path = os.path.join(instance_rundir, sample_root, 'Align')
			predict_path = os.path.join(instance_rundir, sample_root, 'Predict')

			paths = [seq_qc_path, align_path, predict_path]
			for path in paths:
				if not os.path.exists(path): os.makedirs(path)

			file_pair[sample_root] = [forward_data, reverse_data, seq_qc_path, align_path, predict_path]
			sequence_pairs.append(file_pair)

		if workflow_type == 'assembly':

			##
			## Assembly only requires a prediction folder so we do a slightly different thing here
			sample_root = '_'.join(forward_data_name.split('_')[:-1])
			predict_path = os.path.join(instance_rundir, sample_root, 'Predict')
			if not os.path.exists(predict_path): os.makedirs(predict_path)

			file_pair[sample_root] = [forward_data, reverse_data, predict_path]
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
		log.error('{}{}{}{}'.format(Colour.red,'shd__ ',Colour.end,'Specified -b/-c path could not be found.'))
	return False

def check_input_files(input_format, input_file, raise_exception=True):

	if input_file.endswith(input_format):
		return True
	if raise_exception:
		log.error('{}{}{}{}'.format(Colour.red,'shd__ ',Colour.end,'Unrecognised file format found in -b/-c path.'))
	return False

def initialise_libraries(instance_params):

	trigger = False

	##
	## Subfunction for recycling code
	## Calls UNIX which for checking binaries present
	def which(library):
		library_subprocess = subprocess.Popen(['which', library], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
		library_directory = library_subprocess.communicate()
		library_subprocess.wait()
		if not library in library_directory[0]:
			log.critical('{}{}{}{}{}'.format(Colour.red, 'shd__ ', Colour.end, 'Missing library: ', library, '. Not installed or not on $PATH'))
			raise ScaleHDException

	##
	## To determine which binaries to check for
	## AttributeError in the situation where instance_params origin differs
	## try for -c style, except AttributeError for -b style
	try:
		quality_control = instance_params.config_dict['instance_flags']['@quality_control']
		alignment = instance_params.config_dict['instance_flags']['@sequence_alignment']
		genotyping = instance_params.config_dict['instance_flags']['@genotype_prediction']
	except AttributeError:
		quality_control = instance_params['quality_control']
		alignment = instance_params['sequence_alignment']
		genotyping = instance_params['genotype_prediction']

	if quality_control:
		try:which('fastqc')
		except ScaleHDException: trigger=True
		try:which('cutadapt')
		except ScaleHDException: trigger=True
		try:which('sabre')
		except ScaleHDException: trigger=True
	if alignment:
		try:which('bowtie2')
		except ScaleHDException: trigger=True
		try:which('bowtie2-build')
		except ScaleHDException: trigger=True
		try:which('samtools')
		except ScaleHDException: trigger=True
	if genotyping:
		try:which('samtools')
		except ScaleHDException: trigger=True

	return trigger

def sanitise_outputs(output_argument):

	output_root = output_argument[0]

	## Ensures root output is a real directory
	## Generates folder name based on date (for run ident)
	date = datetime.date.today().strftime('%d-%m-%Y')
	walltime = datetime.datetime.now().strftime('%H%M%S')
	today = date + '-' + walltime

	## If the user specified root doesn't exist, make it
	## Then make the run directory for datetime
	if not os.path.exists(output_root):
		log.info('{}{}{}{}'.format(Colour.bold, 'shd__ ', Colour.end, 'Creating output root... '))
		os.mkdir(output_root)
	run_dir = output_root + 'ScaleHDRun_' + today
	log.info('{}{}{}{}'.format(Colour.bold, 'shd__ ', Colour.end, 'Creating instance run directory.. '))
	os.mkdir(run_dir)

	## Inform user it's all gonna be okaaaayyyy
	log.info('{}{}{}{}'.format(Colour.green, 'shd__ ', Colour.end, 'Output directories OK!'))
	return run_dir

def replace_fqfile(mutate_list, target_fqfile, altered_path):

	if target_fqfile in mutate_list:
		loc = mutate_list.index(target_fqfile)
		mutate_list[loc] = altered_path
	return mutate_list