from __future__ import division
import random
import argparse
import sys
import os
from itertools import tee, izip
import matplotlib.pyplot as plt
import numpy as np
import pysam
import difflib

# class subsample:
# 	def __init__(self):
# 		pass
#
# 	parser = argparse.ArgumentParser()
# 	parser.add_argument("input", help="input FASTQ filename")
# 	parser.add_argument("output", help="output FASTQ filename")
# 	parser.add_argument("-f", "--fraction", type=float, help="fraction of reads to sample")
# 	parser.add_argument("-n", "--number", type=int, help="number of reads to sample")
# 	parser.add_argument("-s", "--sample", type=int, help="number of output files to write", default=1)
# 	args = parser.parse_args()
#
# 	if args.fraction and args.number:
# 		sys.exit("give either a fraction or a number, not both")
#
# 	if not args.fraction and not args.number:
# 		sys.exit("you must give either a fraction or a number")
#
# 	print("counting records....")
# 	with open(args.input) as input:
# 		num_lines = sum([1 for line in input])
# 	total_records = int(num_lines / 4)
#
# 	if args.fraction:
# 		args.number = int(total_records * args.fraction)
#
# 	print("sampling " + str(args.number) + " out of " + str(total_records) + " records")
#
# 	output_files = []
# 	output_sequence_sets = []
# 	for i in range(args.sample):
# 		output_files.append(open(args.output + "." + str(i), "w"))
# 		output_sequence_sets.append(set(random.sample(xrange(total_records + 1), args.number)))
#
# 	record_number = 0
# 	with open(args.input) as input:
# 			for line1 in input:
# 				line2 = input.next()
# 				line3 = input.next()
# 				line4 = input.next()
# 				for i, output in enumerate(output_files):
# 					if record_number in output_sequence_sets[i]:
# 							output.write(line1)
# 							output.write(line2)
# 							output.write(line3)
# 							output.write(line4)
# 				record_number += 1
# 				if record_number % 100000 == 0:
# 					print(str((record_number / total_records) * 100)  + " % done")
#
#
# 	for output in output_files:
# 		output.close()
# 	print("done!")

class detect_atypical:
	def __init__(self):
		"""
		Plan::

		use samtools to sort..index..
		"""

		#self.atypical_path = '/Users/alastairm/Documents/Work/ScaleHDTests/Output/IdentAtypical/atypical-701-515/Predict/atypical-701-515_R1/assembly_sorted.bam'
		#self.atypical_path = '/Users/alastairm/Documents/Work/ScaleHDTests/Output/IdentAtypical/atypical-702-505/Predict/atypical-702-505_R1/assembly_sorted.bam'
		self.atypical_path = '/Users/alastairm/Documents/Work/ScaleHDTests/Output/IdentAtypical/atypical-702-518/Predict/atypical-702-518_R1/assembly_sorted.bam'
		self.assembly_object = pysam.AlignmentFile(self.atypical_path, "rb")

		self.top1 = self.assembly_object.fetch(reference='17_1_1_7_2')
		self.top2 = self.assembly_object.fetch(reference='18_1_1_10_2')
		self.present_references = self.assembly_object.references

		##
		## Find number of reads present per reference in assembly
		## Take Top 3 for atypical allele investgation
		# for reference in self.present_references:
		# 	print 'Working on: ', reference
		# 	print 'Reads present: ', self.assembly_object.count(reference)
		# 	print '\n'

		# #
		# # Iterate over reads in the current reference contig
		# for read in self.top2:
		# 	print read.query_alignment_sequence, ' \n'

		test_sequence = 'GCGACCCTGGAAAAGCTGATGAAGGCCTTCGAGTCCCTCAAGTCCTTCCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAACAGCCGCCGCCGCCGCCGCCGCCGCCGCCGCCTCCTCAGCTTCCTCAGCCGCCGCCGCAGGCACAGCCGCTGCT' #701-515 (top2)
		#test_sequence = 'GCGACCCTGGAAAAGCTGATGAAGGCCTTCGAGTCCCTCAAGTCCTTCCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAACAGCCGCCGCCGCCGCCGCCGCCGCCGCCGCCTCCTCAGCTTCCTCAGCCGCCGCCGCAGGCACAGCCGCTGCT' #702-505 (top2)
		#test_sequence = 'GCGACCCTGGAAAAGCTGATGAAGGCCTTCGAGTCCCTCAAGTCCTTCCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAACAGCAACAGCCGCCACCGCCGCCGCCGCCGCCGCCGCCTCCTCCTCAGCTTCCTCAGCCGCCGCCGCAGGCACAGCCGCTGCT' #502-718 (top2)
		n = 3; test_split = [test_sequence[i:i + n] for i in range(0, len(test_sequence), n)]

		##
		## Get repeat regions for CAG and CCG; based on similarity to mask scores
		## Any regions returned that are > in index than the 'cutoff' flag are truncated
		## Combined list used hereafter
		cag_tract = self.get_repeat_region(test_split, 'CAG')
		ccg_tract = self.get_repeat_region(test_split, 'CCG')
		cct_tract = self.get_repeat_region(test_split, 'CCT')
		intv_range = range(cag_tract[-1]+1, ccg_tract[0]); intv_strng = ''

		for i in range(0, len(test_split)):
			if i in cag_tract: print i, test_split[i], ' <<< CAG Tract'
			elif i in ccg_tract: print i, test_split[i], ' <<< CCG Tract'
			elif i in intv_range: print i, test_split[i], ' <<< Intervening'; intv_strng += str(test_split[i])
			elif i in cct_tract: print i, test_split[i], ' <<< CCT Tract'
			else: print i, test_split[i]

		print '\nRepeat tracts: ', cag_tract + ccg_tract
		print 'CAG: ', len(cag_tract)
		print 'CCG: ', len(ccg_tract)
		print 'Intervening: ', intv_range
		print 'Intervening: ', intv_strng

		if intv_strng != 'CAACAGCCGCCA':
			print '\n!! OH SHIT SON ITS AN ATYPICAL ALLELE !!\n'

	def get_repeat_region(self, triplet_input, mask):
		"""
		hehe
		:param triplet_input:
		:param mask:
		:return:
		"""

		##
		## Calculate similarty of every triplet in the current read
		## .. given the current mask we are investigating (CAG or CCG)
		current_tract = []
		for split in triplet_input:
			curr_score = self.similar(split,mask)
			current_tract.append((split,curr_score))

		##
		## Loop over the scores, to find the 'search areas'
		## i.e. rough repeat-tract locations
		repeat_region = []
		for i in range(0, len(current_tract)):
			if current_tract[i][1] == 1.0 and current_tract[i+1][1] == 1.0: search = True
			elif current_tract[i-1][1] == 1.0 and current_tract[i][1] == 1.0: search = True
			else: search = False
			if search: repeat_region.append(i)

		##
		## Find increases for indexes where the difference is >1
		## i.e. "fake" "extra" repeat tract found in addition to legit one
		diff = 0; flagged_idx = 0
		for j in range(0, len(repeat_region)):
			try: diff = abs(repeat_region[j+1]-repeat_region[j])
			except IndexError: pass
			if diff > 1: flagged_idx = repeat_region[j]

		##
		## Iterate over repeat region index list
		## remove any indexes which are above the flagged cut-off
		for index in repeat_region:
			if flagged_idx != 0 and index > flagged_idx:
				repeat_region = [x for x in repeat_region if x != index]

		##
		## Return modified/cleaned list
		return repeat_region

	@staticmethod
	def similar(seq1, seq2):
		return difflib.SequenceMatcher(a=seq1.lower(), b=seq2.lower()).ratio()

if __name__ == '__main__':
	detect_atypical()