from __future__ import division
import pysam
import os
import difflib
import subprocess
import regex
import numpy as np
from collections import Counter
from ..__allelecontainer import IndividualAllele

class ScanAtypical:
	def __init__(self, sequencepair_object, instance_params):
		"""
		Class which utilises basic Digital Signal Processing to determine whether an aligned assembly
		contains any atypical alleles, for the most commonly aligned references contigs. Sometimes,
		variation in the Intevening sequence of the HD repeat tract can have variations which confuse
		alignment algorithms, and reads are incorrectly assigned to a reference that does not truly
		represent the length of the specific repeat tracts. This scans the sample, detects atypical alleles
		and informs the user. If --config, realignment to a specific reference is executed. If --batch,
		the user is informed in the situation of atypical detection; but since the sequence data is missing in batch
		mode, re-alignment is not possible.
		:param input_assembly_tuple: Tuple of (sample_output_folder, specific_assembly_file)
		"""

		##
		## Variables for this class/assembly data
		self.sequence_path = sequencepair_object.get_alignpath()
		self.sorted_assembly = sequencepair_object.get_fwassembly()
		self.instance_params = instance_params
		self.subsample_assembly = None
		self.subsample_index = None
		self.assembly_object = None
		self.present_references = None
		self.assembly_targets = None
		self.atypical_count = 0
		self.atypical_info = {}

		##
		## Fill objects with data
		self.process_assembly()

		##
		## Run the scanning algorithm
		## Exception for (unexpected) EOF
		try: self.scan_reference_reads()
		except StopIteration: self.assembly_object.close()

		##
		## Turn results into objects
		primary_object = IndividualAllele(); secondary_object = IndividualAllele()
		primary_data, secondary_data, atypical_count = self.organise_atypicals()
		sequencepair_object.set_atypical_count(atypical_count)
		for allele_pair in [(primary_object, primary_data, 'PRI'), (secondary_object, secondary_data, 'SEC')]:
			obj = allele_pair[0]; dat = allele_pair[1]
			obj.set_header(allele_pair[2])
			obj.set_allelestatus(dat['Status'])
			obj.set_referencelabel(dat['Reference'])
			obj.set_originalreference(dat['OriginalReference'])
			obj.set_totalreads(dat['TotalReads'])
			obj.set_typicalreads(dat['TypicalCount'])
			obj.set_typicalpcnt(dat['TypicalPcnt'])
			obj.set_atypicalreads(dat['AtypicalCount'])
			obj.set_atypicalpcnt(dat['AtypicalPcnt'])
			obj.set_fiveprime(dat['5PFlank'])
			obj.set_cagval(dat['EstimatedCAG'])
			obj.set_intervening(dat['InterveningSequence'])
			obj.set_caacagval(dat['EstimatedCAACAG'])
			obj.set_ccgccaval(dat['EstimatedCCGCCA'])
			obj.set_ccgval(dat['EstimatedCCG'])
			obj.set_cctval(dat['EstimatedCCT'])
			obj.set_threeprime(dat['3PFlank'])
		sequencepair_object.set_primary_allele(primary_object)
		sequencepair_object.set_secondary_allele(secondary_object)

		##
		## Generate an atypical report for writing
		self.atypical_report = os.path.join(self.sequence_path, 'AtypicalReport.txt')
		report_file = open(self.atypical_report, 'w')
		report_file.write('{}{}\n{}{}\n{}{}\n{}{}'.format('Primary Allele: ', primary_object.get_reflabel(),
														  'Primary Original: ', primary_object.get_originalreference(),
														  'Secondary Allele: ', secondary_object.get_reflabel(),
														  'Secondary Original: ', secondary_object.get_originalreference()))
		report_file.close()

	def process_assembly(self):
		"""
		Function to take the complete aligned assembly and subsample to 20% of the reads (for speed purposes).
		Will implement user specified subsample threshold eventually. Once subsampled, pass to processing.
		:return: None
		"""

		##
		## Determine number of reads - for subsampling float
		awk = ['awk', ' {i+=$3} END {print i}']
		count_process = subprocess.Popen(['samtools','idxstats', self.sorted_assembly], stdout=subprocess.PIPE)
		awk_process = subprocess.Popen(awk, stdin=count_process.stdout, stdout=subprocess.PIPE)
		count_process.wait(); awk_process.wait(); awk_output = int(awk_process.communicate()[0])
		if awk_output > 20000: subsample_float = 0.05
		else: subsample_float = 0.2

		##
		## Subsample reads
		## Index the subsampled assembly
		self.subsample_assembly = os.path.join(self.sequence_path,'subsample.sam')
		self.subsample_index = os.path.join(self.sequence_path,'subsample.sam.bai')
		assem_obj = open(self.subsample_assembly,'w')
		subsample_process = subprocess.Popen(['samtools','view','-s',str(subsample_float),'-b', self.sorted_assembly], stdout=assem_obj)
		subsample_process.wait(); assem_obj.close()
		index_process = subprocess.Popen(['samtools','index',self.subsample_assembly]); index_process.wait()

		##
		## Load into object, determine references
		self.assembly_object = pysam.AlignmentFile(self.subsample_assembly, 'rb')
		self.present_references = self.assembly_object.references
		assembly_refdat = []
		for reference in self.present_references:
			reference_tuple = (reference, self.assembly_object.count(reference))
			if reference_tuple[1] == 0: pass
			else: assembly_refdat.append(reference_tuple)
		self.assembly_targets = sorted(assembly_refdat, key=lambda x:x[1], reverse=True)[0:3]

	@staticmethod
	def typical_rotation(input_string):
		size1 = len('CAACAGCCGCCA')
		size2 = len(input_string)
		temp = ''

		# Check if sizes of two strings are same
		if size1 != size2: return 0

		# Create a temp string with value str1.str1
		temp = 'CAACAGCCGCCA' + 'CAACAGCCGCCA'

		# Now check if str2 is a substring of temp
		# string.count returns the number of occurences of
		# the second string in temp
		if temp.count(input_string) > 0: return 1
		else: return 0

	def scan_reference_reads(self):

		##
		## Iterate over top 3 aligned references in this assembly
		## Fetch the reads aligned to the current reference
		for investigation in self.assembly_targets:
			reference_data = self.assembly_object.fetch(reference=investigation[0])

			##
			## Counts of atypical/typical reads
			typical_count = 0; atypical_count = 0; reference_atypicals = []; fp_flanks = []; tp_flanks = []
			ref_cag = []; ref_ccg = []; ref_cct = []

			##
			## For every read in this reference, get the aligned sequence
			## Split into triplet sliding window list, remove any triplets that are < 3
			for read in reference_data:
				target_sequence = read.query_alignment_sequence
				sequence_windows = [target_sequence[i:i + 3] for i in range(0, len(target_sequence), 3)]
				sequence_windows = [x for x in sequence_windows if len(x)==3]

				##
				## Get repeat regions for CAG and CCG; based on similarity mask scores for the current window
				## Any regions returned that are (idx > end_of_region) are truncated
				## CAG and CCG repeat region index list combined
				cag_masks = ['CAG', 'AGC', 'GCA']
				ccg_masks = ['CCG', 'CGC', 'GCC']
				cct_masks = ['CCT', 'CTC', 'TCC']

				##
				## CAG/CCG Masking
				## Sort all list of tuples by length of second element (repeat tract length)
				## Select first item as 'true' tract, then calculate intervening sequence length
				cag_tracts = []; ccg_tracts = []; cct_tracts = []
				try:
					for mask in cag_masks: cag_tracts.append((mask, self.get_repeat_tract(sequence_windows, mask)))
					for mask in ccg_masks: ccg_tracts.append((mask, self.get_repeat_tract(sequence_windows, mask)))
					cag_tract = sorted(cag_tracts, key=lambda a: len(a[1]), reverse=True)[0][1]
					ccg_tract = sorted(ccg_tracts, key=lambda a: len(a[1]), reverse=True)[0][1]

					##
					## CCT Masking/Intervening calculation
					intervene_string = ''; fp_flank_string = ''; tp_flank_string = ''
					for mask in cct_masks: cct_tracts.append((mask, self.get_cct_tract(sequence_windows, mask, ccg_tract[-1])))
					cct_tract = sorted(cct_tracts, key=lambda a: len(a[1]), reverse=True)[0][1]
					intervene_range = range(cag_tract[-1]+1, ccg_tract[0])
					fp_flank_range = range(0, cag_tract[0]-1)
					tp_flank_range = range(cct_tract[-1]+1, len(sequence_windows))
				except IndexError:
					continue

				##
				## Add length to reference-run
				ref_cag.append(len(cag_tract))
				ref_ccg.append(len(ccg_tract))
				ref_cct.append(len(cct_tract))

				##
				## Count fp flank occurrences
				for i in range(0, len(sequence_windows)):
					if i in fp_flank_range:
						fp_flank_string += str(sequence_windows[i])
				fp_flanks.append(fp_flank_string)

				##
				## Count tp flank occurrences
				for i in range (0, len(sequence_windows)):
					if i in tp_flank_range:
						tp_flank_string += str(sequence_windows[i])
				tp_flanks.append(tp_flank_string)

				##
				## Atypical Detection
				for i in range(0, len(sequence_windows)):
					if i in intervene_range:
						intervene_string += str(sequence_windows[i])
				if self.rotation_check('CAACAGCCGCCA', intervene_string):
					intervene_string = 'CAACAGCCGCCA'
				if intervene_string != 'CAACAGCCGCCA':
					atypical_count += 1
					reference_atypicals.append(intervene_string)
				else:
					typical_count += 1

			##
			## Calculate the presence of each 'state' of reference
			ref_typical = format(((typical_count / investigation[1]) * 100), '.2f')
			ref_atypical = format(((atypical_count / investigation[1]) * 100), '.2f')
			est_cag = Counter(ref_cag).most_common()[0][0]
			est_ccg = Counter(ref_ccg).most_common()[0][0]
			est_cct = Counter(ref_cct).most_common()[0][0]

			##
			## Determine most frequent intervening sequence
			atypical_population = Counter(reference_atypicals).most_common()
			fp_flank_population = Counter(fp_flanks).most_common()
			tp_flank_population = Counter(tp_flanks).most_common()

			if len(atypical_population) == 0: atypical_population = [['CAACAGCCGCCA']]
			reference_dictionary = {'TotalReads':investigation[1],
									'TypicalCount': typical_count,
									'TypicalPcnt': ref_typical,
									'AtypicalCount': atypical_count,
									'AtypicalPcnt': ref_atypical,
									'Status':self.atypical_count,
									'5PFlank':fp_flank_population[0][0],
									'3PFlank':tp_flank_population[0][0],
									'EstimatedCAG': est_cag,
									'EstimatedCCG': est_ccg,
									'EstimatedCCT': est_cct,
									'InterveningSequence': atypical_population[0][0]}

			if atypical_count > typical_count:
				self.atypical_count += 1
				reference_dictionary['Status'] = 'Atypical'
			elif est_cct != 2:
				self.atypical_count += 1
				reference_dictionary['Status'] = 'Atypical'
			else:
				reference_dictionary['Status'] = 'Typical'
				reference_dictionary['InterveningSequence'] = 'CAACAGCCGCCA'

			##
			## If the intervening is longer in #2, assume poor sequencing in #1 and use #2
			try:
				if len(atypical_population[0][0]) < len(atypical_population[1][0]):
					if reference_dictionary['Status'] == 'Typical':
						reference_dictionary['InterveningSequence'] = max([atypical_population[0][0],atypical_population[1][0]], key=len)
			except IndexError:
				reference_dictionary['InterveningSequence'] = atypical_population[0][0]

			##
			## Check for mismatch just before intervening sequence
			try:
				top_hit = atypical_population[0][1]; second_hit = atypical_population[1][1]
				diff = ((top_hit-second_hit)/top_hit)*100
				if diff < 30.00:
					if len(atypical_population[0][0]) == 15:
						if np.isclose(self.similar('CAG', atypical_population[0][0][0:3]), [0.66], atol=0.1):
							reference_dictionary['InterveningSequence'] = 'CAACAGCCGCCA'
			except IndexError:
				pass

			##
			## Append results to reference label
			self.atypical_info[investigation[0]] = reference_dictionary

		os.remove(self.subsample_assembly)
		os.remove(self.subsample_index)

	def get_repeat_tract(self, triplet_input, mask):

		##
		## Score the entire read against the current mask
		current_tract = []
		for split in triplet_input:
			curr_score = self.similar(split,mask)
			current_tract.append((split,curr_score))

		##
		## Anchors
		region_start = None; region_end = None
		## Find the beginning of the CAG tract..
		## assuming streak of 3, confidence high in real start
		for i in range(0, len(current_tract)):
			try:
				if current_tract[i][1] == 1.0:
					if not region_start:
						if current_tract[i+1][1] == 1.0 and current_tract[i+2][1] == 1.0:
							region_start = i
				if current_tract[i][1] == 1.0:
					region_end = i
			except IndexError:
				pass

		##
		## If typeerror (i.e. one of the regions was None.. no start was found)
		## return empty list as there is no repeat tract for this mask
		try: first_pass_range = range(region_start, region_end+1)
		except TypeError: return []

		##
		## Loop over rough range, remove items where n-1,n+1 and n+2 are not good matches for current mask
		for j in first_pass_range:
			if not current_tract[j][1] == 1.0:
				sub_score = 0
				try:
					for sub_check in [current_tract[j-1], current_tract[j+1], current_tract[j+2]]:
						if sub_check[1] == 1.0: sub_score += 1
				except IndexError:
					pass
				if sub_score != 3: first_pass_range = [x for x in first_pass_range if x!=j]

		##
		## Some downstream matches may exist still so..
		## remove anything outside of >1 streak in pass
		diff = 0; flagged_idx = 0
		for k in range(0, len(first_pass_range)):
			try: diff = abs(first_pass_range[k+1]-first_pass_range[k])
			except IndexError: pass
			if diff > 1 and flagged_idx == 0: flagged_idx = first_pass_range[k]+1
		for index in first_pass_range:
			if flagged_idx != 0 and index > flagged_idx:
				first_pass_range = [x for x in first_pass_range if x!= index]

		##
		## Return list to call
		return first_pass_range

	def get_cct_tract(self, triplet_input, mask, anchor):

		##
		## Get all triplets after the end of the CCG tract (anchor)
		post_anchor = []
		for i in range(0, len(triplet_input)):
			if i > anchor: post_anchor.append((i, triplet_input[i]))

		##
		## If similarity matches the current mask, add that index to tract
		cct_tract = []
		for item in post_anchor:
			if self.similar(mask, item[1]) == 1.0:
				cct_tract.append(item[0])

		##
		## Remove indexes in tract list if difference between indexes > 1 (gaps dont happen in cct)
		diff = 0; flagged_idx = 0
		for i in range(0, len(cct_tract)):
			try: diff = abs(cct_tract[i+1]-cct_tract[i])
			except IndexError: pass
			if diff > 1 and flagged_idx == 0: flagged_idx = cct_tract[i]+1
		for index in cct_tract:
			if flagged_idx!=0 and index>flagged_idx:
				cct_tract = [x for x in cct_tract if x!=index]

		##
		## Return
		return cct_tract

	def organise_atypicals(self):

		##
		## Constructs
		sorted_info = sorted(self.atypical_info.iteritems(), key=lambda (x, y): y['TotalReads'], reverse=True)
		if len(sorted_info) != 3: raise IndexError('< 3 references in sorted top; alignment failure?')

		##
		## Check % dropoff in read count between #2 and #3
		beta_diff = float(abs(sorted_info[1][1]['TotalReads'] - sorted_info[2][1]['TotalReads']))
		beta_drop = float(beta_diff / sorted_info[1][1]['TotalReads'])

		##
		## TODO fix this garbage dumpster fire
		## TODO oh my god it's so ugly
		## Top1 always used
		## Check Top2 vs Top3, if drop is > 20%, ok
		## If not >= 20%, investigate further:
		##   If T2_CCG == T3_CCG, potential slippage event
		##     if T2_CAG-T3_CAG > 1, no slippage, use Top1 and Top2
		##     if T2_CAG-T3_CAG == 1, maybe slippage
		##       if T2-T3 read drop is <= 0.1, use Top1 and Top3, slippage event
		##       otherwise, use Top1 and Top2, no slippage
		##   If T2_CCG and T3_CCG differ, no slippage
		##     Check T1_CAG and T3_CAG for n-1, check diff between 2 and 3 > 10%
		primary_allele = sorted_info[0][1]; primary_allele['Reference'] = sorted_info[0][0]
		if not beta_drop >= 0.20:
			##
			## CCG match, potential slippage
			if sorted_info[1][1]['EstimatedCCG'] == sorted_info[2][1]['EstimatedCCG']:
				if abs(sorted_info[1][1]['EstimatedCAG'] - sorted_info[2][1]['EstimatedCAG']) > 1:
					secondary_allele = sorted_info[1][1]
					secondary_allele['Reference'] = sorted_info[1][0]
				elif abs(sorted_info[1][1]['EstimatedCAG'] - sorted_info[2][1]['EstimatedCAG']) == 1:
					if beta_drop <= 0.1:
						secondary_allele = sorted_info[2][1]
						secondary_allele['Reference'] = sorted_info[2][0]
					else:
						secondary_allele = sorted_info[1][1]
						secondary_allele['Reference'] = sorted_info[1][0]
				else:
					secondary_allele = sorted_info[1][1]
					secondary_allele['Reference'] = sorted_info[1][0]
			##
			## CCG don't match, no slippage
			else:
				if abs(sorted_info[0][1]['EstimatedCAG'] - sorted_info[2][1]['EstimatedCAG']) == 1 and beta_drop >= 0.10:
					secondary_allele = sorted_info[1][1]
					secondary_allele['Reference'] = sorted_info[1][0]
				else:
					secondary_allele = sorted_info[2][1]
					secondary_allele['Reference'] = sorted_info[2][0]
		##
		## Otherwise the drop between #2 and #3 is so big that #3 isn't possible
		else:
			secondary_allele = sorted_info[1][1]; secondary_allele['Reference'] = sorted_info[1][0]

		##
		## For each of the alleles we've determined..
		## Get intervening lengths, create accurate genotype string
		atypical_count = 0
		for allele in [primary_allele, secondary_allele]:
			new_genotype, caacag_count, ccgcca_count = self.create_genotype_label(allele)
			allele['OriginalReference'] = allele['Reference']
			allele['Reference'] = new_genotype
			allele['EstimatedCAACAG'] = caacag_count
			allele['EstimatedCCGCCA'] = ccgcca_count
			if allele['Status'] == 'Atypical': atypical_count += 1

		return primary_allele, secondary_allele, atypical_count

	def create_genotype_label(self, input_reference):

		intervening = input_reference['InterveningSequence']
		intervening_freq = Counter(list((intervening[0 + i:6 + i] for i in range(0, len(intervening), 6)))).items()
		caacag_count = 0; ccgcca_count = 0
		caacag_flag = False; ccgcca_flag = False

		##
		## TODO fix this dumpster garbage shit
		## TODO oh my god it's so ugly
		if len(intervening_freq) < 1: caacag_count = 1; ccgcca_count = 1
		if len(intervening_freq) < 2: caacag_count = 1; ccgcca_count = 1

		##
		## Check CAACAG
		try:
			caacag_freq = intervening_freq[0]
			if self.rotation_check('CAACAG', caacag_freq[0]):
				caacag_count = caacag_freq[1]
			else:
				if input_reference['Status'] == 'Typical':
					caacag_count = 1; ccgcca_count = 1
		except IndexError:
			caacag_flag = True

		##
		## Check CCGCCA
		try:
			ccgcca_freq = intervening_freq[1]
			if self.rotation_check('CCGCCA', ccgcca_freq[0]):
				ccgcca_count = ccgcca_freq[1]
			else:
				if input_reference['Status'] == 'Typical':
					caacag_count = 1; ccgcca_count = 1
		except IndexError:
			ccgcca_flag = True

		##
		## Parse flags in event of error
		if caacag_flag: caacag_count = 0; ccgcca_count = 0
		if ccgcca_flag:
			if not caacag_flag:
				caacag_count = 1; ccgcca_count = 0
			else:
				caacag_count = 0; ccgcca_count = 0

		##
		## Safety check
		if input_reference['Status'] == 'Typical' and (caacag_count != 1 or ccgcca_count != 1):
			caacag_count = 1; ccgcca_count = 1

		genotype_label = '{}_{}_{}_{}_{}'.format(input_reference['EstimatedCAG'], caacag_count, ccgcca_count,
												 input_reference['EstimatedCCG'], input_reference['EstimatedCCT'])
		return genotype_label, caacag_count, ccgcca_count

	def get_atypicalreport(self):
		return self.atypical_report

	@staticmethod
	def rotation_check(string1, string2):
		size1 = len(string1)
		size2 = len(string2)

		# Check if sizes of two strings are same
		if size1 != size2: return 0

		# Create a temp string with value str1.str1
		temp = string1 + string1

		# Now check if str2 is a substring of temp (with s = 1 mismatch)
		rotation_match = regex.findall(r"(?:" + string2 + "){s<=1}", temp, regex.BESTMATCH)

		if len(rotation_match) > 0:
			return 1
		else:
			return 0

	@staticmethod
	def similar(seq1, seq2):
		return difflib.SequenceMatcher(a=seq1.lower(), b=seq2.lower()).ratio()