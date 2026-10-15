#/usr/bin/python

import os
import vcf
import subprocess
import numpy as np
import logging as log
from ..__backend import Colour as clr

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
		# freebayes -f ref.fa -C snp_threshold assembly.bam > variants.vcf

	def generate_variant_data(self):

		"""
		Simple workflow function which calls various third party tools in order to call SNPs within the alignment
		which we created earlier on within the pipeline.
		:return: n/a
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
				indiv_typical_reference_name = fw_idx.split('/')[-2:-1][0]
				dict_path = '/'.join(fw_idx.split('/')[:-1])+'/'+indiv_typical_reference_name+'.dict'
				if not os.path.isfile(dict_path):
					## picard dict creation
					picard_string = 'picard {} {}'.format(fw_idx, dict_path)
					picard_subprocess = subprocess.Popen([picard_string], shell=True,
					 stdout=subprocess.PIPE, stderr=subprocess.PIPE)
					picard_log = picard_subprocess.communicate(); picard_subprocess.wait()
					if 'ERROR' in picard_log:
						log.error('{}{}{}{}'.format(clr.red, 'shd__ ', clr.end, 'Failure in PICARD. Check sample output log.'))
						logpath = os.path.join(predpath, 'PicardErrorLog.txt')
						with open(logpath, 'w') as logfi:
							logfi.write(picard_log[0])
							logfi.write(picard_log[1])

			# If allele is atypical then we generated our own reference for the INTV struct
			# need to index and dict each specific reference for these types of alleles
			if allele.get_allelestatus() == 'Atypical':
				indiv_atypical_reference_name = fw_idx.split('/')[-2:-1][0]
				dict_path = '/'.join(fw_idx.split('/')[:-1])+'/'+indiv_atypical_reference_name+'.dict'
				## picard dict creation
				picard_string = 'picard {} {}'.format(fw_idx, dict_path)
				picard_subprocess = subprocess.Popen([picard_string], shell=True,
					stdout=subprocess.PIPE, stderr=subprocess.PIPE)
				picard_log = picard_subprocess.communicate(); picard_subprocess.wait()
				if 'ERROR' in picard_log:
					log.error('{}{}{}{}'.format(clr.red, 'shd__ ', clr.end, 'Failure in PICARD. Check sample output log.'))
					logpath = os.path.join(predpath, 'PicardErrorLog.txt')
					with open(logpath, 'w') as logfi:
						logfi.write(picard_log[0])
						logfi.write(picard_log[1])

			## freebayes haplotype caller
			observation_threshold = self.sequencepair_object.get_snpobservationvalue()
			freebayes_output = os.path.join(predpath, '{}_FreeBayesVariantCall.vcf'.format(header))
			freebayes_string = 'freebayes -f {} -B 4000 -C {} {}'.format(
				fw_idx, observation_threshold, fw_assembly)
			freebayes_outfi = open(freebayes_output, 'w')
			freebayes_subprocess = subprocess.Popen([freebayes_string], shell=True,
													stdout=freebayes_outfi, stderr=subprocess.PIPE)
			freebayes_log = freebayes_subprocess.communicate(); freebayes_subprocess.wait()
			if 'ERROR' in freebayes_log:
				log.error('{}{}{}{}'.format(clr.red, 'shd__ ', clr.end, 'Failure in FreeBayes. Check sample output log.'))
				logpath = os.path.join(predpath, 'FreeBayesErrorLog.txt')
				with open(logpath, 'w') as logfi:
					logfi.write(freebayes_log[0])
					logfi.write(freebayes_log[1])
			allele.set_freebayes_file(freebayes_output)
			freebayes_outfi.close()

	def scrape_relevance(self):
		"""
		Given the user can prioritise which variant calling software to prefer, this function will extract
		the relevant information from associated objects here. Validity of a SNP is filtered via certain
		thresholds which can be specified by the user in the configuration XML.
		:return: nothin m8
		"""

		##
		## Todo:: generalise the code for freebayes/gatk variant calling to appear at least somewhat professional
		## For each allele we scrape the data required to call SNPs and assign to their own object variables
		variant_cutoff = int(self.instance_params.config_dict['prediction_flags']['@quality_cutoff'])
		for allele in [self.sequencepair_object.get_primaryallele(), self.sequencepair_object.get_secondaryallele()]:

			## Get variants found by freebayes
			target = ''; freebayes_call = 'N/A'; freebayes_score = 0
			freebayes_matched = []; freebayes_unmatched = []
			if allele.get_allelestatus() == 'Typical': target = allele.get_reflabel()
			if allele.get_allelestatus() == 'Atypical': target = allele.get_reflabel().split('CAG')[0]
			freebayes_reader = vcf.Reader(open(allele.get_freebayes_file(), 'r'))
			for record in freebayes_reader:
				origin = ''
				if allele.get_allelestatus() == 'Typical': origin = record.CHROM
				if allele.get_allelestatus() == 'Atypical': origin = record.CHROM.split('CAG')[0][:-1]
				if origin == target:
					freebayes_matched.append(record)
				else:
					freebayes_unmatched.append(record)

			## sort and remove records which are < user specified cutoff
			## todo again generalise this code you absolute throbber
			freebayes_sorted = sorted(freebayes_matched, key=lambda a:a.QUAL, reverse=True)
			freebayes_sorted = [x for x in freebayes_sorted if x.QUAL > variant_cutoff]

			##
			## PCR amplification results in fake mutations from the primer sequence
			## filter out anything that is outside of our region of interest i.e. smallest & largest (positions)
			positions = [x.POS for x in freebayes_sorted]
			if not len(freebayes_sorted) == 0:
				smallest = min([x for x in positions]); largest = max([x for x in positions])
				freebayes_sorted = [x for x in freebayes_sorted if not np.isclose([x.POS], [smallest], atol=5)[0]]
				freebayes_sorted = [x for x in freebayes_sorted if not np.isclose([x.POS], [largest], atol=5)[0]]

			## Determine what to set values of call/score to, then apply to allele object
			## will be written to InstanceReport.csv from whatever algo the user wanted
			if not len(freebayes_sorted) == 0:
				## we have snps!
				freebayes_str = ''; mutation_calls = ''; mutation_scores = ''
				for mutation in freebayes_sorted:
					mutation_calls += '{}->{}@{}   '.format(mutation.REF, mutation.ALT[0], mutation.POS)
					mutation_scores += '{}   '.format(mutation.QUAL)
				allele.set_variantcall(mutation_calls)
				allele.set_variantscore(freebayes_score)
			else:
				## we do not
				allele.set_variantcall(freebayes_call)
				allele.set_variantscore(freebayes_score)

			## Write unmatched to file
			target_dir = os.path.join(self.sequencepair_object.get_predictpath(), 'IrrelevantVariants.txt')
			with open(target_dir, 'w') as outfi:
				for record in freebayes_unmatched:
					record_str = 'Freebayes: {} = {} -> {} @ {}. Qual: {}'.format(record.CHROM, record.REF,
																		record.ALT, record.POS, record.QUAL)
					outfi.write(record_str+'\n')
				outfi.write('\n')

	def set_report(self, input_report):
		self.snp_report = input_report
	def get_report(self):
		return self.snp_report
