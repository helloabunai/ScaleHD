##
## Generals
import glob
import os
import sys
import subprocess
import logging as log

##
## Backend junk
from ..backpack import Colour as clr
from multiprocessing import cpu_count

THREADS = str(cpu_count())
processed_files = []

class SeqQC:

	def __init__(self, input_data, instance_rundir, stage, instance_params):
		self.input_data = input_data
		self.instance_rundir = instance_rundir
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
		for fqfile in glob.glob(os.path.join(self.input_data, '*')):
			if fqfile.endswith('.fq') or fqfile.endswith('.fastq') or fqfile.endswith('.fq.gz') or fqfile.endswith('.fastq.gz'):
				return True
		if raise_exception:
			log.error('{}{}{}{}'.format(clr.red,'shd__ ',clr.end,'I/O: Non-FQ/FastQ file found in input directory.'))
		return False

	def execute_fastQC(self):
		for fqfile in glob.glob(os.path.join(self.input_data, '*')):
			sample_root = fqfile.split('/')[-1].split('.')[0] ##absolutely_disgusting.jpg
			sample_outdir = os.path.join(self.instance_rundir, sample_root, 'SeqQC')
			if not os.path.exists(sample_outdir): os.makedirs(sample_outdir)
			fastqc_outdir = os.path.join(sample_outdir, 'FastQC')
			if not os.path.exists(fastqc_outdir): os.makedirs(fastqc_outdir)

			fastqc_process = subprocess.Popen(['fastqc','-q','--extract','-t', THREADS, '-o', fastqc_outdir, fqfile])
			fastqc_process.wait()

	def execute_demultiplex(self):
		##TODO this would be first, but do later
		pass

	def execute_trimming(self):

		##
		## Generic subprocess for use in trimming
		def execute_cutadapt(arguments_split, sample_output):
			trimming_subprocess = subprocess.Popen(['cutadapt'] + arguments_split, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
			cutadapt_raw_output = trimming_subprocess.communicate()
			cutadapt_report = cutadapt_raw_output[0]
			cutadapt_errors = cutadapt_raw_output[1]
			trimming_subprocess.wait()

			report_directory = os.path.join(sample_output,'TrimmingReport.txt')
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
				for fqfile in glob.glob(os.path.join(self.input_data, '*')):
					sample_root = fqfile.split('/')[-1].split('.')[0] ##absolutely_disgusting.jpg
					trimmed_name = sample_root + '_trimmed.fastq'
					sample_outdir = os.path.join(self.instance_rundir, sample_root, 'SeqQC')
					trimmed_outdir = os.path.join(sample_outdir,trimmed_name)
					quality_threshold = self.instance_params.config_dict['trim_flags']['@quality_threshold']

					argument_list = ['-q', quality_threshold, fqfile, '-o', trimmed_outdir]
					execute_cutadapt(argument_list, sample_outdir)
					processed_files.append(trimmed_outdir)

			if trim_type.lower()=='adapter':
				for fqfile in glob.glob(os.path.join(self.input_data, '*')):
					sample_root = fqfile.split('/')[-1].split('.')[0] ##absolutely_disgusting.jpg
					trimmed_name = sample_root + '_trimmed.fastq'
					sample_outdir = os.path.join(self.instance_rundir, sample_root, 'SeqQC')
					trimmed_outdir = os.path.join(sample_outdir,trimmed_name)
					adapter_anchor = self.instance_params.config_dict['trim_flags']['@adapter_flag']
					adapter_string = self.instance_params.config_dict['trim_flags']['@adapter']

					##
					## Alter string based on anchor, messy but whatever
					if adapter_anchor == '-a$':adapter_anchor = '-a';adapter_string += '$'
					if adapter_anchor == '-g^':adapter_anchor = '-g';adapter_string = '^' + adapter_string

					argument_list = [adapter_anchor, adapter_string, fqfile, '-o', trimmed_outdir]
					execute_cutadapt(argument_list, sample_outdir)
					processed_files.append(trimmed_outdir)

			if trim_type.lower()=='both':
				for fqfile in glob.glob(os.path.join(self.input_data, '*')):
					sample_root = fqfile.split('/')[-1].split('.')[0] ##absolutely_disgusting.jpg
					trimmed_name = sample_root + '_trimmed.fastq'
					sample_outdir = os.path.join(self.instance_rundir, sample_root, 'SeqQC')
					trimmed_outdir = os.path.join(sample_outdir,trimmed_name)

					quality_threshold = self.instance_params.config_dict['trim_flags']['@quality_threshold']
					adapter_anchor = self.instance_params.config_dict['trim_flags']['@adapter_flag']
					adapter_string = self.instance_params.config_dict['trim_flags']['@adapter']

					##
					## Alter string based on anchor, messy but whatever
					if adapter_anchor == '-a$':adapter_anchor = '-a';adapter_string += '$'
					if adapter_anchor == '-g^':adapter_anchor = '-g';adapter_string = '^' + adapter_string

					argument_list = ['-q', quality_threshold, adapter_anchor, adapter_string, fqfile, '-o', trimmed_outdir]
					execute_cutadapt(argument_list, sample_outdir)
					processed_files.append(trimmed_outdir)

		if self.trimming_errors == 'True':
			log.error('{}{}{}{}'.format(clr.red,'shd__ ',clr.end,'Trimming errors occurred. Check logging report!'))
			sys.exit(2)

	@staticmethod
	def getProcessedFiles():
		return processed_files







