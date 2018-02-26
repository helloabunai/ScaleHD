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
		self.call_variants()

		# samtools faidx <reference>
		# picard CreateSequenceDictioanry REFERENCE=<reference> OUTPUT=<reference.dict>
		# samtools index <assembly> <assembly.bam.bai>
		# gatk -R 4k-HD-INTER.fa -T HaplotypeCaller -I assembly_sorted.bam -o variants.vcf

	def call_variants(self):

		"""
		Simple workflow function which calls various third party tools in order to call SNPs within the alignment
		which we created earlier on within the pipeline.
		:return: fuck all
		"""

		for allele in [self.sequencepair_object.get_primaryallele(), self.sequencepair_object.get_secondaryallele()]:

			fw_idx = allele.get_fwidx(); fw_assembly = allele.get_fwassembly()
			header = allele.get_header(); predpath = self.sequencepair_object.get_predictpath()

			## non bwa-mem index of reference
			faidx_subprocess = subprocess.Popen(['samtools', 'faidx', fw_idx],
												stdout=subprocess.PIPE, stderr=subprocess.PIPE)
			faidx_subprocess.wait()

			## create a sequence dictionary for use in GATK
			dict_path = '/'.join(fw_idx.split('/')[:-1]) + '.dict'
			picard_subprocess = subprocess.Popen(['picard', 'CreateSequenceDictionary',
												  'REFERENCE={}'.format(fw_idx),
												  'OUTPUT={}'.format(dict_path)],
												 stdout=subprocess.PIPE, stderr=subprocess.PIPE)
			picard_log = picard_subprocess.communicate(); picard_subprocess.wait()
			## utilised HaplotypeCaller in GATK to determine variants
			desired_output = os.path.join(predpath, '{}_variants.vcf'.format(header))
			gatk_subprocess = subprocess.Popen(['gatk', '-R', fw_idx, '-T', 'HaplotypeCaller', '-I', fw_assembly
												, '-o', desired_output], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
			gatk_log = gatk_subprocess.communicate(); gatk_subprocess.wait()

	def set_report(self, input_report):
		self.snp_report = input_report
	def get_report(self):
		return self.snp_report
