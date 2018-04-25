#/usr/bin/python
__version__ = 0.3
__author__ = 'alastair.maxwell@glasgow.ac.uk'

import os
import subprocess

class DetermineMutations:
	def __init__(self, sequencepair_object, instance_params):
		self.sequencepair_object = sequencepair_object
		self.instance_params = instance_params
		self.snp_report = ''
		self.generate_variant_data()
		self.scrape_relevance()

		# samtools faidx <reference>
		# picard CreateSequenceDictioanry REFERENCE=<reference> OUTPUT=<reference.dict>
		# samtools index <assembly> <assembly.bam.bai>
		# gatk -R 4k-HD-INTER.fa -T HaplotypeCaller -I assembly_sorted.bam -o variants.vcf

	def generate_variant_data(self):

		"""
		Simple workflow function which calls various third party tools in order to call SNPs within the alignment
		which we created earlier on within the pipeline.
	
		For some damn reason while working at home on arch linux I need to call certain subprocesses
		via an interactive bash session rather than using shell == True. fuck knows why

		:return: nothing yet
		"""	
		for allele in [self.sequencepair_object.get_primaryallele(), self.sequencepair_object.get_secondaryallele()]:
			## get data
			fw_idx = allele.get_fwidx(); fw_assembly = allele.get_fwassembly()
			header = allele.get_header(); predpath = self.sequencepair_object.get_predictpath()
			faidx_subprocess = subprocess.Popen(['samtools', 'faidx', fw_idx],
			stdout=subprocess.PIPE, stderr=subprocess.PIPE)
			faidx_subprocess.wait()

			## if the allele is typical, we'll be using the standard reference
			## generically used by all typical samples so don't need to create dict repeatedly
			if allele.get_allelestatus() == "Typical":
				## check for dict, if doesn't exist, create
				dict_path = '/'.join(fw_idx.split('/')[:-1])+'/forward.dict'
				if not os.path.isfile(dict_path):
					## picard dict creation
					picard_string = 'picard {} {}'.format(fw_idx, dict_path)
					picard_subprocess = subprocess.Popen([picard_string], shell=True,
					 stdout=subprocess.PIPE, stderr=subprocess.PIPE) 
					picard_log = picard_subprocess.communicate(); picard_subprocess.wait()
					##todo error checking in picard_log

			# If allele is atypical then we generated our own reference for the INTV struct
			# need to index and dict each specific reference for these types of alleles
			if allele.get_allelestatus() == 'Atypical':
				indiv_atypical_reference_name = fw_idx.split('/')[-2:-1][0]
				dict_path = '/'.join(fw_idx.split('/')[:-1])+'/'+indiv_atypical_reference_name+'.dict'
				picard_string = 'picard {} {}'.format(fw_idx, dict_path)
				picard_subprocess = subprocess.Popen([picard_string], shell=True,
					stdout=subprocess.PIPE, stderr=subprocess.PIPE)
				picard_log = picard_subprocess.communicate(); picard_subprocess.wait()
				##todo error checking in picard_log
			
			# gatk haplotype caller
			desired_output = os.path.join(predpath, '{}_variants.vcf'.format(header))
			gatk_string = 'gatk HaplotypeCaller -R {} -I {} -O {}'.format(fw_idx, fw_assembly, desired_output)
			gatk_subprocess = subprocess.Popen([gatk_string], shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
			gatk_log = gatk_subprocess.communicate(); gatk_subprocess.wait()
			##todo error checking in gatk_log
			allele.set_variant_file(desired_output)

	def scrape_relevance(self):
		print 'hi scrape_relevance()'

		for allele in [self.sequencepair_object.get_primaryallele(), self.sequencepair_object.get_secondaryallele()]:
			print allele.get_variant_file()

	def set_report(self, input_report):
		self.snp_report = input_report
	def get_report(self):
		return self.snp_report
