#/usr/bin/python
__version__ = 0.01
__author__ = 'alastair.maxwell@glasgow.ac.uk'

##
## Generals
import os
import sys
import subprocess
import logging as log

##
## Backend junk
from ..__backend import Colour as clr
from ..__backend import replace_fqfile
from multiprocessing import cpu_count

THREADS = str(cpu_count())

class SeqQC:

	def __init__(self, sequencepair_data, target_output, stage, instance_params):
		self.sequencepair_data = sequencepair_data
		self.input_filepair = [sequencepair_data[0], sequencepair_data[1]]
		self.target_output = target_output
		self.instance_params = instance_params
		self.trimming_errors = False

		if stage.lower()=='valid':
			self.verify_input()
			self.execute_fastQC()
		if stage.lower()=='dmpx':
			self.execute_demultiplex()
		if stage.lower()=='trim':
			self.execute_trimming()

	def verify_input(self, raise_exception=True):

		for fqfile in self.input_filepair:
			if fqfile.endswith('.fq') or fqfile.endswith('.fastq') or fqfile.endswith('.fq.gz') or fqfile.endswith('.fastq.gz'):
				return True

		if raise_exception:
			log.error('{}{}{}{}'.format(clr.red,'shd__ ',clr.end,'I/O: Invalid file format detected in input. Check input data.'))
		return False

	def execute_fastQC(self):
		##
		## For the files in the current file pair, make FastQC output folder and run FastQC
		for fqfile in self.input_filepair:
			fastqc_outdir = os.path.join(self.target_output, 'FastQC')
			if not os.path.exists(fastqc_outdir): os.makedirs(fastqc_outdir)
			fastqc_process = subprocess.Popen(['fastqc','--quiet','--extract','-t',THREADS,'-o',fastqc_outdir,fqfile], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
			fastqc_process.wait()

	def execute_demultiplex(self):
		##TODO this would be first, but do later
		pass

	def execute_trimming(self):

		##
		## Generic subprocess for use in trimming
		def execute_cutadapt(arguments_split, filename_root, sample_output):
			trimming_subprocess = subprocess.Popen(['cutadapt'] + arguments_split, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
			cutadapt_raw_output = trimming_subprocess.communicate()
			cutadapt_report = cutadapt_raw_output[0]
			cutadapt_errors = cutadapt_raw_output[1]
			trimming_subprocess.wait()

			report_directory = os.path.join(sample_output, filename_root + '_TrimmingReport.txt')
			report_file = open(report_directory, 'w')
			report_file.write(cutadapt_report)
			report_file.write(cutadapt_errors)
			report_file.close()

			if cutadapt_errors is not None:
				self.trimming_errors = True

		##
		## Determine what we want to trim from parameters dictionary
		## Be paranoid, do if test for trimming just incase
		## Then go into setting up instance
		if self.instance_params.config_dict['trim_flags']['@trim_data']:
			trim_type = self.instance_params.config_dict['trim_flags']['@trim_type']

			if trim_type.lower()=='quality':
				for i in range(0,len(self.input_filepair)):
					file_root = self.input_filepair[i].split('/')[-1].split('.')[0] ##absolutely_disgusting.jpg
					trimmed_outdir = '{}/{}{}{}'.format(self.target_output,'trimmed_',file_root,'.fastq')
					quality_threshold = self.instance_params.config_dict['trim_flags']['@quality_threshold']

					argument_list = ['-q', quality_threshold, self.input_filepair[i], '-o', trimmed_outdir]
					execute_cutadapt(argument_list, file_root, self.target_output)
					self.sequencepair_data = replace_fqfile(self.sequencepair_data, self.input_filepair[i], trimmed_outdir)

			if trim_type.lower()=='adapter':
				for i in range(0,len(self.input_filepair)):
					file_root = self.input_filepair[i].split('/')[-1].split('.')[0] ##absolutely_disgusting.jpg
					trimmed_outdir = '{}/{}{}{}'.format(self.target_output,'trimmed_',file_root,'.fastq')
					adapter_anchor = self.instance_params.config_dict['trim_flags']['@adapter_flag']
					adapter_string = self.instance_params.config_dict['trim_flags']['@adapter']

					##
					## Alter string based on anchor, messy but whatever
					if adapter_anchor == '-a$':adapter_anchor = '-a';adapter_string += '$'
					if adapter_anchor == '-g^':adapter_anchor = '-g';adapter_string = '^' + adapter_string

					argument_list = [adapter_anchor, adapter_string, self.input_filepair[i], '-o', trimmed_outdir]
					execute_cutadapt(argument_list, file_root, self.target_output)
					self.sequencepair_data = replace_fqfile(self.sequencepair_data, self.input_filepair[i], trimmed_outdir)

			if trim_type.lower()=='both':
				for i in range(0,len(self.input_filepair)):
					file_root = self.input_filepair[i].split('/')[-1].split('.')[0] ##absolutely_disgusting.jpg
					trimmed_outdir = '{}/{}{}{}'.format(self.target_output,'trimmed_',file_root,'.fastq')
					quality_threshold = self.instance_params.config_dict['trim_flags']['@quality_threshold']
					adapter_anchor = self.instance_params.config_dict['trim_flags']['@adapter_flag']
					adapter_string = self.instance_params.config_dict['trim_flags']['@adapter']

					##
					## Alter string based on anchor, messy but whatever
					if adapter_anchor == '-a$':adapter_anchor = '-a';adapter_string += '$'
					if adapter_anchor == '-g^':adapter_anchor = '-g';adapter_string = '^' + adapter_string

					argument_list = ['-q', quality_threshold, adapter_anchor, adapter_string, self.input_filepair[i], '-o', trimmed_outdir]
					execute_cutadapt(argument_list, file_root, self.target_output)
					self.sequencepair_data = replace_fqfile(self.sequencepair_data, self.input_filepair[i], trimmed_outdir)

		if self.trimming_errors == 'True':
			log.error('{}{}{}{}'.format(clr.red,'shd__ ',clr.end,'Trimming errors occurred. Check logging report!'))
			sys.exit(2)