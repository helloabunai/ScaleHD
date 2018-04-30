#/usr/bin/python
__version__ = 0.311
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
from ..__backend import mkdir_p
from multiprocessing import cpu_count

THREADS = str(cpu_count())
TR_REPORT = []

class SeqQC:

	def __init__(self, sequencepair_object, instance_params, stage=None):
		self.sequencepair_data = sequencepair_object
		self.input_filepair = [sequencepair_object.get_fwreads(), sequencepair_object.get_rvreads()]
		self.target_output = sequencepair_object.get_qcpath()
		self.instance_params = instance_params
		self.trimming_errors = False
		self.trimming_report = []
		if stage.lower()=='validate': self.verify_input()
		if stage.lower()=='trim': self.execute_trimming(); self.execute_fastqc()

	def verify_input(self, raise_exception=True):

		for fqfile in self.input_filepair:
			if fqfile.endswith('.fq') or fqfile.endswith('.fastq') or fqfile.endswith('.fq.gz') or fqfile.endswith('.fastq.gz'):
				return True

		if raise_exception:
			log.error('{}{}{}{}'.format(clr.red,'shd__ ',clr.end,'I/O: Invalid file format detected in input. Check input data.'))
		return False

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
			return report_directory

		##
		## Determine what we want to trim from parameters dictionary
		## Be paranoid, do if test for trimming just incase
		## Then go into setting up instance
		if self.instance_params.config_dict['instance_flags']['@quality_control']:
			trim_type = self.instance_params.config_dict['trim_flags']['@trim_type']
			error_tolerance = self.instance_params.config_dict['trim_flags']['@error_tolerance']

			if trim_type.lower()=='quality':
				for i in range(0,len(self.input_filepair)):
					file_root = self.input_filepair[i].split('/')[-1].split('.')[0] ##absolutely_disgusting.jpg
					trimmed_outdir = '{}/{}{}{}'.format(self.target_output,'trimmed_',file_root,'.fq')
					quality_threshold = self.instance_params.config_dict['trim_flags']['@quality_threshold']

					argument_list = ['-e', error_tolerance, '-q', quality_threshold, self.input_filepair[i], '-o', trimmed_outdir]
					trim_report = execute_cutadapt(argument_list, file_root, self.target_output)
					if i == 0: self.sequencepair_data.set_fwreads(trimmed_outdir)
					if i == 1: self.sequencepair_data.set_rvreads(trimmed_outdir)
					self.trimming_report.append(trim_report)

			stepwise_counter = 0
			if trim_type.lower()=='adapter':
				for i in range(0,len(self.input_filepair)):
					file_root = self.input_filepair[i].split('/')[-1].split('.')[0] ##absolutely_disgusting.jpg
					trimmed_outdir = '{}/{}{}{}'.format(self.target_output,'trimmed_',file_root,'.fq')
					adapter_anchor = self.instance_params.config_dict['trim_flags']['@adapter_flag']
					adapter_string = ''
					if stepwise_counter == 0:
						adapter_string = self.instance_params.config_dict['trim_flags']['@forward_adapter']
					if stepwise_counter == 1:
						adapter_string = self.instance_params.config_dict['trim_flags']['@reverse_adapter']

					##
					## Alter string based on anchor, messy but whatever
					if adapter_anchor == '-a$':adapter_anchor = '-a';adapter_string += '$'
					if adapter_anchor == '-g^':adapter_anchor = '-g';adapter_string = '^' + adapter_string

					argument_list = ['-e', error_tolerance, adapter_anchor, adapter_string, self.input_filepair[i], '-o', trimmed_outdir]
					trim_report = execute_cutadapt(argument_list, file_root, self.target_output)
					if i == 0: self.sequencepair_data.set_fwreads(trimmed_outdir)
					if i == 1: self.sequencepair_data.set_rvreads(trimmed_outdir)
					self.trimming_report.append(trim_report)
					stepwise_counter += 1

			stepwise_counter = 0
			if trim_type.lower()=='both':
				for i in range(0,len(self.input_filepair)):
					file_root = self.input_filepair[i].split('/')[-1].split('.')[0] ##absolutely_disgusting.jpg
					trimmed_outdir = '{}/{}{}{}'.format(self.target_output,'trimmed_',file_root,'.fq')
					quality_threshold = self.instance_params.config_dict['trim_flags']['@quality_threshold']
					adapter_anchor = self.instance_params.config_dict['trim_flags']['@adapter_flag']
					adapter_string = ''
					if stepwise_counter == 0:
						adapter_string = self.instance_params.config_dict['trim_flags']['@forward_adapter']
					if stepwise_counter == 1:
						adapter_string = self.instance_params.config_dict['trim_flags']['@reverse_adapter']

					##
					## Alter string based on anchor, messy but whatever
					if adapter_anchor == '-a$':adapter_anchor = '-a';adapter_string += '$'
					if adapter_anchor == '-g^':adapter_anchor = '-g';adapter_string = '^' + adapter_string

					argument_list = ['-e', error_tolerance, '-q', quality_threshold, adapter_anchor, adapter_string, self.input_filepair[i], '-o', trimmed_outdir]
					trim_report = execute_cutadapt(argument_list, file_root, self.target_output)
					if i == 0: self.sequencepair_data.set_fwreads(trimmed_outdir)
					if i == 1: self.sequencepair_data.set_rvreads(trimmed_outdir)
					self.trimming_report.append(trim_report)

		if self.trimming_errors == 'True':
			log.error('{}{}{}{}'.format(clr.red,'shd__ ',clr.end,'Trimming errors occurred. Check logging report!'))
			sys.exit(2)

	def execute_fastqc(self):

		##
		## For the files in the current file pair, make FastQC output folder and run FastQC
		for fqfile in self.input_filepair[0:1]:
			fastqc_outdir = os.path.join(self.target_output, 'FastQC')
			mkdir_p(fastqc_outdir)
			fastqc_process = subprocess.Popen(['fastqc','--quiet','--extract','-t',THREADS,'-o',fastqc_outdir,fqfile], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
			fastqc_process.wait()

	def get_trimreport(self):
		return self.trimming_report