from __future__ import division
##
##Imports
import os
import re
import regex
import pysam
import difflib
import subprocess
import numpy as np
import logging as log
import multiprocessing
from operator import itemgetter
from collections import Counter
from ..__backend import Colour as clr
from ..seq_qc.__quality_control import THREADS
from ..__allelecontainer import IndividualAllele

## required for worker thread
## (can't serialise class-bound methods)
def rotation_check(string1, string2):
	size1 = len(string1)
	size2 = len(string2)

	# Check if sizes of two strings are same
	if size1 != size2: return False

	# Create a temp string with value str1.str1
	temp = string1 + string1

	# Now check if str2 is a substring of temp (with s = 1 mismatch)
	rotation_match = regex.findall(r"(?:" + string2 + "){s<=1}", temp, regex.BESTMATCH)

	if len(rotation_match) > 0:
		return True
	else:
		return False

def get_repeat_tract(triplet_input, mask):

	##
	## Score the entire read against the current mask
	current_tract = []
	for split in triplet_input:
		curr_score = similar(split, mask)
		current_tract.append((split, curr_score))

	##
	## Anchors
	region_start = None; region_end = None
	## Find the beginning of the CAG tract..
	## assuming streak of 3, confidence high in real start
	for i in range(0, len(current_tract)):
		try:
			if current_tract[i][1] == 1.0:
				if not region_start:
					if current_tract[i + 1][1] == 1.0 and current_tract[i + 2][1] == 1.0:
						region_start = i
			if current_tract[i][1] == 1.0:
				region_end = i
		except IndexError:
			pass

	##
	## If typeerror (i.e. one of the regions was None.. no start was found)
	## return empty list as there is no repeat tract for this mask
	try:
		first_pass_range = range(region_start, region_end + 1)
	except TypeError:
		return []

	##
	## Loop over rough range, remove items where n-1,n+1 and n+2 are not good matches for current mask
	for j in first_pass_range:
		if not current_tract[j][1] == 1.0:
			sub_score = 0
			try:
				for sub_check in [current_tract[j - 1], current_tract[j + 1], current_tract[j + 2]]:
					if sub_check[1] == 1.0: sub_score += 1
			except IndexError:
				pass
			if sub_score != 3: first_pass_range = [x for x in first_pass_range if x != j]

	##
	## Some downstream matches may exist still so..
	## remove anything outside of >1 streak in pass
	diff = 0; flagged_idx = 0
	for k in range(0, len(first_pass_range)):
		try:
			diff = abs(first_pass_range[k + 1] - first_pass_range[k])
		except IndexError:
			pass
		if diff > 1 and flagged_idx == 0: flagged_idx = first_pass_range[k] + 1
	for index in first_pass_range:
		if flagged_idx != 0 and index > flagged_idx:
			first_pass_range = [x for x in first_pass_range if x != index]

	##
	## Return list to call
	return first_pass_range

def get_cct_tract(triplet_input, mask, anchor):

	##
	## Get all triplets after the end of the CCG tract (anchor)
	post_anchor = []
	for i in range(0, len(triplet_input)):
		if i > anchor: post_anchor.append((i, triplet_input[i]))

	##
	## If similarity matches the current mask, add that index to tract
	cct_tract = []
	for item in post_anchor:
		if similar(mask, item[1]) == 1.0:
			cct_tract.append(item[0])

	##
	## Remove indexes in tract list if difference between indexes > 1 (gaps dont happen in cct)
	diff = 0; flagged_idx = 0
	for i in range(0, len(cct_tract)):
		try:
			diff = abs(cct_tract[i + 1] - cct_tract[i])
		except IndexError:
			pass
		if diff > 1 and flagged_idx == 0: flagged_idx = cct_tract[i] + 1
	for index in cct_tract:
		if flagged_idx != 0 and index > flagged_idx:
			cct_tract = [x for x in cct_tract if x != index]

	##
	## Return
	return cct_tract

def similar(seq1, seq2):
	return difflib.SequenceMatcher(a=seq1.lower(), b=seq2.lower()).ratio()

## worker thread for processor assign
def scan_reference_reads(current_iterator):
	"""
	Function which determines the literal repeat regions, ignoring misalignment issues.
	We loop over every 'investigation' from this assembly <- i.e. the top 3 reference (in terms of read count)
	Each read within each reference is then scanned, to determine the structure of said read.
	:return:
	"""

	##
	## Unpack iterator values into discrete objects for manipulation
	contig_vector = current_iterator
	contig = contig_vector[0]; read_count = contig_vector[1]; assembly_path = contig_vector[2]

	##
	## object of SAM
	assembly_object = pysam.AlignmentFile(assembly_path, 'rb')

	##
	## Counts of atypical/typical reads
	typical_count = 0; atypical_count = 0; intervening_population = []
	fp_flanks = []; tp_flanks = []; ref_cag = []; ref_ccg = []; ref_cct = []

	##
	## For every read in this reference, get the aligned sequence
	## Split into triplet sliding window list, remove any triplets that are < 3
	for read in assembly_object.fetch(reference=contig):
		target_sequence = read.query_alignment_sequence
		sequence_windows = [target_sequence[i:i + 3] for i in range(0, len(target_sequence), 3)]
		sequence_windows = [x for x in sequence_windows if len(x) == 3]

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
			for mask in cag_masks: cag_tracts.append((mask, get_repeat_tract(sequence_windows, mask)))
			for mask in ccg_masks: ccg_tracts.append((mask, get_repeat_tract(sequence_windows, mask)))
			cag_tract = sorted(cag_tracts, key=lambda a: len(a[1]), reverse=True)[0][1]
			ccg_tract = sorted(ccg_tracts, key=lambda a: len(a[1]), reverse=True)[0][1]

			##
			## CCT Masking/Intervening calculation
			intervene_string = ''; fp_flank_string = ''; tp_flank_string = ''
			for mask in cct_masks: cct_tracts.append((mask, get_cct_tract(sequence_windows, mask, ccg_tract[-1])))
			cct_tract = sorted(cct_tracts, key=lambda a: len(a[1]), reverse=True)[0][1]
			intervene_range = range(cag_tract[-1] + 1, ccg_tract[0])
			fp_flank_range = range(0, cag_tract[0] - 1)
			tp_flank_range = range(cct_tract[-1] + 1, len(sequence_windows))
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
		for i in range(0, len(sequence_windows)):
			if i in tp_flank_range:
				tp_flank_string += str(sequence_windows[i])
		tp_flanks.append(tp_flank_string)

		##
		## Atypical Detection
		for i in range(0, len(sequence_windows)):
			if i in intervene_range:
				intervene_string += str(sequence_windows[i])
		if rotation_check('CAACAGCCGCCA', intervene_string): intervene_string = 'CAACAGCCGCCA'
		if intervene_string != 'CAACAGCCGCCA': atypical_count += 1
		else: typical_count += 1
		intervening_population.append(intervene_string)

	##
	## Calculate the presence of each 'state' of reference
	ref_typical = format(((typical_count / current_iterator[1]) * 100), '.2f')
	ref_atypical = format(((atypical_count / current_iterator[1]) * 100), '.2f')
	est_cag = Counter(ref_cag).most_common()[0][0]
	est_ccg = Counter(ref_ccg).most_common()[0][0]
	est_cct = Counter(ref_cct).most_common()[0][0]
	## cct fuckery
	cct_test = Counter(ref_cct).most_common()
	try:
		cct_diff = float(cct_test[0][1]) / float(cct_test[1][1])
	except IndexError:
		cct_diff = 0
	if np.isclose([cct_diff], [2.0], atol=0.3):
		est_cct = 2

	##
	## Determine most frequent intervening sequence
	common_intervening = Counter(intervening_population).most_common()
	fp_flank_population = Counter(fp_flanks).most_common()
	tp_flank_population = Counter(tp_flanks).most_common()
	single_counter = 0

	typical_tenpcnt = (typical_count / 100) * 5
	if np.isclose([typical_count], [atypical_count], atol=typical_tenpcnt):
		raise Exception('Allele(s) nearing 50/50 split atypical/typical read count.')

	if len(common_intervening) == 0: common_intervening = [['CAACAGCCGCCA']]
	reference_dictionary = {'TotalReads': current_iterator[1],
							'TypicalCount': typical_count,
							'TypicalPcnt': ref_typical,
							'AtypicalCount': atypical_count,
							'AtypicalPcnt': ref_atypical,
							'Status': '',
							'5PFlank': fp_flank_population[0][0],
							'3PFlank': tp_flank_population[0][0],
							'EstimatedCAG': est_cag,
							'EstimatedCCG': est_ccg,
							'EstimatedCCT': est_cct,
							'InterveningSequence': common_intervening[0][0]}

	if atypical_count > typical_count:
		reference_dictionary['Status'] = 'Atypical'
	elif est_cct != 2:
		reference_dictionary['Status'] = 'Atypical'
	else:
		reference_dictionary['Status'] = 'Typical'
		reference_dictionary['InterveningSequence'] = 'CAACAGCCGCCA'

	##
	## If the intervening is longer in #2, assume poor sequencing in #1 and use #2
	try:
		if len(common_intervening[0][0]) < len(common_intervening[1][0]):
			if reference_dictionary['Status'] == 'Typical':
				reference_dictionary['InterveningSequence'] = max([common_intervening[0][0], common_intervening[1][0]],
																  key=len)
	except IndexError:
		reference_dictionary['InterveningSequence'] = common_intervening[0][0]

	##
	## Check for mismatch just before intervening sequence
	try:
		top_hit = common_intervening[0][1]; second_hit = common_intervening[1][1]
		diff = ((top_hit - second_hit) / top_hit) * 100
		if diff < 30.00:
			if len(common_intervening[0][0]) == 15:
				if np.isclose(similar('CAG', common_intervening[0][0][0:3]), [0.66], atol=0.1):
					reference_dictionary['InterveningSequence'] = 'CAACAGCCGCCA'
	except IndexError:
		pass

	##
	## If all scanned items are only counted once
	for item in common_intervening:
		try:
			if item[1] == 1: single_counter += 1
		except IndexError:
			if ref_typical > ref_atypical:
				reference_dictionary['Status'] = 'Typical'
				reference_dictionary['InterveningSequence'] = 'CAACAGCCGCCA'
	if single_counter == len(common_intervening) and reference_dictionary['Status'] == 'Typical':
		reference_dictionary['InterveningSequence'] = 'CAACAGCCGCCA'

	##
	## Append results to reference label
	return [contig, reference_dictionary]

class ScanAtypical:
	def __init__(self, sequencepair_object, instance_params):
		"""
		Class which utilises Digital Signal Processing to determine the repeat tract structure of the current sample.
		General overview:
			Subsample aligned SAM file (from __alignment.py) to increase speed of function
			Loop over the references with the top 3 highest number of aligned reads (only ones relevant to genotyping)
			For each read in that reference, scan for regions using rotating masks
			Return regions, determine values
			Assign values to appropriate allele object class variables
			Return
		:param sequencepair_object: Object of the current Sequence pair being processed.. see __allelecontainer.py
		:param instance_params: Dictionary of config settings from the input XML document
		"""

		##
		## Variables for this class/assembly data
		self.sequencepair_object = sequencepair_object
		self.sequence_path = sequencepair_object.get_alignpath()
		self.sorted_assembly = sequencepair_object.get_fwassembly()
		self.instance_params = instance_params
		self.subsample_flag = sequencepair_object.get_subsampleflag()
		self.subsample_assembly = None
		self.subsample_index = None
		self.assembly_object = None
		self.present_references = None
		self.assembly_targets = None
		self.atypical_count = 0
		self.awk_output = 0
		self.atypical_info = {}
		self.alignment_warning = False

		##
		## Placeholder class object for worker pool
		## Fill objects with data
		self.processor_pool = None
		self.process_assembly()

		##
		## Run the scanning algorithm
		## Exception for (unexpected) EOF
		try:
			self.prepare_process_workers()
		except StopIteration:
			self.assembly_object.close()

		##
		## Turn results into objects
		primary_object = IndividualAllele(); secondary_object = IndividualAllele()
		primary_data, secondary_data, atypical_count = self.organise_atypicals()

		sequencepair_object.set_atypical_count(atypical_count)
		for allele_pair in [(primary_object, primary_data, 'PRI'), (secondary_object, secondary_data, 'SEC')]:
			obj = allele_pair[0]; dat = allele_pair[1]
			obj.set_header(allele_pair[2])
			obj.set_allelestatus(dat.get('Status'))
			obj.set_referencelabel(dat.get('Reference'))
			obj.set_originalreference(dat.get('OriginalReference'))
			obj.set_totalreads(dat.get('TotalReads'))
			obj.set_typicalreads(dat.get('TypicalCount'))
			obj.set_typicalpcnt(dat.get('TypicalPcnt'))
			obj.set_atypicalreads(dat.get('AtypicalCount'))
			obj.set_atypicalpcnt(dat.get('AtypicalPcnt'))
			obj.set_fiveprime(dat.get('5PFlank'))
			obj.set_cagval(dat.get('EstimatedCAG'))
			obj.set_intervening(dat.get('InterveningSequence'))
			obj.set_caacagval(dat.get('EstimatedCAACAG'))
			obj.set_ccgccaval(dat.get('EstimatedCCGCCA'))
			obj.set_ccgval(dat.get('EstimatedCCG'))
			obj.set_cctval(dat.get('EstimatedCCT'))
			obj.set_threeprime(dat.get('3PFlank'))
			obj.set_rewrittenccg(dat.get('RewrittenCCG'))
			obj.set_unrewrittenccg(dat.get('UnrewrittenCCG'))
			obj.set_differential_confusion(dat.get('DiffConfuse'))
			obj.set_neighbouring_candidate(dat.get('Neighbouring'))
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
		Function which processes the input SAM for atypical scanning.
		Determine the number of total reads present (for subsampling).
		Read file into PySAM object.. process further
		:return: None
		"""

		##
		## Determine number of reads - for subsampling float
		## Use awk to read samtools idxstats output (get total read count)
		awk = ['awk', ' {i+=$3} END {print i}']
		count_process = subprocess.Popen(['samtools', 'idxstats', self.sorted_assembly], stdout=subprocess.PIPE)
		awk_process = subprocess.Popen(awk, stdin=count_process.stdout, stdout=subprocess.PIPE)
		count_process.wait(); awk_process.wait(); awk_output = int(awk_process.communicate()[0])
		if awk_output > 20000: subsample_float = 0.350
		elif 20000 > awk_output > 15000: subsample_float = 0.450
		elif 15000 > awk_output > 10000: subsample_float = 0.550
		elif 10000 > awk_output > 5000: subsample_float = 0.650
		else: subsample_float = 1.0
		self.sequencepair_object.set_subsampled_fqcount(awk_output)

		##
		## Subsample reads
		## Index the subsampled assembly
		if not self.sequencepair_object.get_broadflag():
			self.sequencepair_object.set_subsampleflag(subsample_float)
			self.sequencepair_object.set_automatic_DSPsubsample(True)
			self.subsample_assembly = os.path.join(self.sequence_path, 'subsample.sam')
			self.subsample_index = os.path.join(self.sequence_path, 'subsample.sam.bai')
			assem_obj = open(self.subsample_assembly, 'w')
			subsample_process = subprocess.Popen(
				['samtools', 'view', '-s', str(subsample_float), '-b', self.sorted_assembly], stdout=assem_obj)
			subsample_process.wait(); assem_obj.close()
			index_process = subprocess.Popen(['samtools', 'index', self.subsample_assembly]); index_process.wait()
		else:
			self.subsample_assembly = self.sorted_assembly

		##
		## Load into object, determine references to investigate
		self.assembly_object = pysam.AlignmentFile(self.subsample_assembly, 'rb')
		self.present_references = self.assembly_object.references
		assembly_refdat = []
		for reference in self.present_references:
			reference_vector = [reference, self.assembly_object.count(reference)]
			if reference_vector[1] == 0: pass
			else: assembly_refdat.append(reference_vector)

		##
		## Assign our target references (Top 3 sorted references, sorted by read count)
		self.assembly_targets = sorted(assembly_refdat, key=itemgetter(1), reverse=True)[0:3]
		for contig in self.assembly_targets:
			temp = self.subsample_assembly
			contig.append(temp)

		##
		## Check for a mal-aligned sample
		fail_score = 0
		for target in self.assembly_targets:
			if target[1] < 20:
				fail_score += 1
		if fail_score != 0:
			log.warning('{}{}{}{}'.format(clr.red, 'shd__ ', clr.end, 'Alignment contains too few reads. Cannot guarantee precision.'))
			self.alignment_warning = True
			self.sequencepair_object.set_alignmentwarning(self.alignment_warning)

	@staticmethod
	def typical_rotation(input_string):

		"""
		Function to detect if an intervening sequence (typical structure) is located within all possible
		rotations of a given/derived intervening sequence. Easiest method is to double to target, and search.
		:param input_string:
		:return:
		"""

		##
		## Lengths of target strings
		size1 = len('CAACAGCCGCCA')
		size2 = len(input_string)
		temp = ''

		##
		## Size equality comparison
		if size1 != size2: return 0

		##
		## Duplicate string to encompass all possible rotations
		temp = 'CAACAGCCGCCA' + 'CAACAGCCGCCA'

		##
		## Now check STR2 is a substring of temp expansion
		## .count() returns number of occurrences of second string in temp
		if temp.count(input_string) > 0: return 1
		else: return 0

	def prepare_process_workers(self):

		"""
		Function which wraps scan_reference_reads inside a multi-process handler pool
		Just speeds up the process of determining HTT structure, assigning one
		investigation in our assembly targets to a discrete processor each
		:return:
		"""

		##
		## Check read count per ref
		if self.assembly_targets[0][1] < 200:
			self.sequencepair_object.set_fatalreadallele(True)
			raise Exception('<200 aligned reads in Allele #1. Data un-usable.')
		if self.assembly_targets[1][1] < 100:
			self.sequencepair_object.set_fatalreadallele(True)
			raise Exception('<100 aligned reads in Allele #2. Data un-usable.')
		if self.assembly_targets[2][1] < 50:
			if np.isclose([self.assembly_targets[2][1]],[50],atol=5):
				pass
			else:
				self.sequencepair_object.set_fatalreadallele(True)
				raise Exception('<50 aligned reads in Allele #3. Data un-usable.')

		##
		## Determine if the user's system has enough processors for the work pool
		# if int(THREADS) >= 3: self.processor_pool = multiprocessing.Pool(3)
		# else: self.processor_pool = multiprocessing.Pool(1)
		if int(THREADS) >= 3: self.processor_pool = multiprocessing.Pool(3)
		else: self.processor_pool = multiprocessing.Pool(1)

		##
		## for each worker we have, give it one 'contig' in our target_assembly to work over
		## make iterator for our processor pool to share
		## run that iterator on our pool of worker processes
		allele_iterator = iter(self.assembly_targets)
		pool_output = self.processor_pool.map(scan_reference_reads, allele_iterator)
		## assign the output of each contig result to the respective dictionary
		for contig_results in pool_output:
			self.atypical_info[contig_results[0]] = contig_results[1]
		self.processor_pool.close()

		## If broadflag = false, we subsampled
		## remove the subsampled SAM file, retaining only the original
		if not self.sequencepair_object.get_broadflag():
			try:
				os.remove(self.subsample_assembly)
				os.remove(self.subsample_index)
			except TypeError:
				pass

	def organise_atypicals(self):

		##
		## Constructs
		sorted_info = sorted(self.atypical_info.iteritems(), key=lambda (x, y): y['TotalReads'], reverse=True)
		if len(sorted_info) != 3: raise IndexError('< 3 references in sorted top; alignment failure?')

		## Top1 always used
		## Secondary reassigned after filtering with heuristics
		primary_allele = sorted_info[0][1]; primary_allele['Reference'] = sorted_info[0][0]
		secondary_allele = primary_allele;secondary_was_set = False

		## Heuristics
		## Estimated CCG
		alpha_estCCG = int(sorted_info[0][1]['EstimatedCCG'])
		beta_estCCG = int(sorted_info[1][1]['EstimatedCCG'])
		theta_estCCG = int(sorted_info[2][1]['EstimatedCCG'])

		## Estimated CAG
		alpha_estCAG = int(sorted_info[0][1]['EstimatedCAG'])
		beta_estCAG = int(sorted_info[1][1]['EstimatedCAG'])
		theta_estCAG = int(sorted_info[2][1]['EstimatedCAG'])

		## CAG differences
		alpha_beta_CAGDiff = abs(alpha_estCAG - beta_estCAG)
		beta_theta_CAGDiff = abs(beta_estCAG - theta_estCAG)
		alpha_theta_CAGDiff = abs(alpha_estCAG - theta_estCAG)

		## Read Count
		alpha_readCount = int(sorted_info[0][1]['TotalReads'])
		beta_readCount = int(sorted_info[1][1]['TotalReads'])
		theta_readCount = int(sorted_info[2][1]['TotalReads'])

		## Literal read drop values
		alpha_beta_ReadDelta = abs(alpha_readCount - beta_readCount)
		beta_theta_ReadDelta = abs(beta_readCount - theta_readCount)

		## Percentage read drops
		alpha_beta_ReadPcnt = alpha_beta_ReadDelta / alpha_readCount
		beta_theta_ReadPcnt = beta_theta_ReadDelta / beta_readCount

		## Begin filtering...
		uniform_ccg = 0
		for allele in sorted_info[1:]:
			if allele[1]['EstimatedCCG'] == primary_allele['EstimatedCCG']:
				uniform_ccg += 1

		## CCG in Top3 all equal?
		if uniform_ccg == 2:
			##top1-top2 CAG difference
			if alpha_beta_CAGDiff == 1 and beta_theta_CAGDiff != 1:
				if alpha_beta_CAGDiff == 1 and alpha_theta_CAGDiff == 1:
					if np.isclose([alpha_beta_ReadPcnt], [0.3], atol=0.1):
						## B.T.ReadPcnt > threshold.. beta == neighbouring
						if np.isclose([beta_theta_ReadPcnt],[0.75],atol=0.1):
							secondary_allele = sorted_info[1][1]
							secondary_allele['Reference'] = sorted_info[1][0]
							secondary_allele['Neighbouring'] = True
							secondary_was_set = True
						## dropoff too small, homozygous haplotype
						if 0 < beta_theta_ReadPcnt <= 0.64:
							secondary_allele = primary_allele.copy()
							secondary_was_set = True
					else:
						## B.T.ReadPcnt > threshold.. beta == neighbouring
						if np.isclose([beta_theta_ReadPcnt],[0.75],atol=0.1):
							secondary_allele = sorted_info[1][1]
							secondary_allele['Reference'] = sorted_info[1][0]
							secondary_allele['Neighbouring'] = True
							secondary_was_set = True
						## dropoff too small, homozygous haplotype
						if 0 < beta_theta_ReadPcnt <= 0.64:
							secondary_allele = primary_allele.copy()
							secondary_allele['DiffConfuse'] = True
							secondary_was_set = True
				if alpha_beta_CAGDiff == 1 and alpha_theta_CAGDiff != 1:
					secondary_allele = sorted_info[2][1]
					secondary_allele['Reference'] = sorted_info[2][0]
					secondary_was_set = True

			##top2-top3 CAG difference
			if beta_theta_CAGDiff == 1 and alpha_beta_CAGDiff != 1:
				if 0.55 <= beta_theta_ReadPcnt < 1:
					secondary_allele = sorted_info[1][1]
					secondary_allele['Reference'] = sorted_info[1][0]
					secondary_was_set = True
				if np.isclose([beta_theta_ReadPcnt],[0.45],atol=0.1):
					secondary_allele = sorted_info[1][1]
					secondary_allele['Reference'] = sorted_info[1][0]
					secondary_was_set = True
				elif beta_theta_ReadPcnt <= 0.35:
					secondary_allele = sorted_info[1][1]
					secondary_allele['Reference'] = sorted_info[1][0]
					secondary_allele['DiffConfuse'] = True
					secondary_was_set = True

			##top1-top2-top3 all within 1
			if beta_theta_CAGDiff == 1 and alpha_beta_CAGDiff == 1:
				if 0.55 <= alpha_beta_ReadPcnt < 1:
					secondary_allele = primary_allele.copy()
					secondary_was_set = True
				elif np.isclose([alpha_beta_ReadPcnt],[0.45],atol=0.1):
					secondary_allele = sorted_info[1][1]
					secondary_allele['Reference'] = sorted_info[1][0]
					secondary_was_set = True
				elif alpha_beta_ReadPcnt <= 0.35:
					secondary_allele = sorted_info[1][1]
					secondary_allele['Reference'] = sorted_info[1][0]
					secondary_allele['DiffConfuse'] = True
					secondary_was_set = True

			## Same CCG but out of step with order
			if alpha_theta_CAGDiff == 1:
				secondary_allele = sorted_info[1][1]
				secondary_allele['Reference'] = sorted_info[1][0]
				secondary_allele['DiffConfuse'] = True
				secondary_was_set = True

			## cell line DNA, skewed expansion with broad peak
			## theta is not neighbouring of beta, but beta is legitimate
			if alpha_beta_CAGDiff > 1 and np.isclose([beta_estCAG], [theta_estCAG], atol=5):
				secondary_allele = sorted_info[1][1]
				secondary_allele['Reference'] = sorted_info[1][0]
				secondary_allele['DiffConfuse'] = True
				secondary_was_set = True

		## CCG in Top3 not equal?
		if uniform_ccg < 2:
			## the CCG in #2,#3 match
			if beta_estCCG == theta_estCCG:
				## determine CAG distance
				## if "neighbours"
				if beta_theta_CAGDiff == 1:
					## read dropoff percentage of #2 determines whether #2 or #3 is legitimate
					if np.isclose([beta_theta_ReadPcnt], [0.15], atol=0.1):
						secondary_allele = sorted_info[1][1]
						secondary_allele['Reference'] = sorted_info[1][0]
						secondary_allele['DiffConfuse'] = True
						secondary_was_set = True
					else:
						secondary_allele = sorted_info[1][1]
						secondary_allele['Reference'] = sorted_info[1][0]
						secondary_was_set = True
				## theta is close to beta, but not neighbouring.. cell line DNA/broad expansion
				if np.isclose([beta_theta_CAGDiff], [1], atol=5):
					secondary_allele = sorted_info[1][1]
					secondary_allele['Reference'] = sorted_info[1][0]
					secondary_allele['DiffConfuse'] = True
					secondary_was_set = True

			## the CCG in #2,#3 don't match
			if not beta_estCCG == theta_estCCG:
				## beta CCG match and alpha.CAG-beta.CAG == 1:
				if alpha_estCCG == beta_estCCG:
					if alpha_beta_CAGDiff == 1:
						## large enough drop between alpha and beta rules out beta/theta
						if np.isclose([alpha_beta_ReadPcnt],[0.75],atol=0.25):
							if abs(beta_estCCG - theta_estCCG) >= 1:
								if np.isclose([beta_theta_ReadPcnt], [0.15], atol=0.15):
									secondary_allele = sorted_info[2][1]
									secondary_allele['Reference'] = sorted_info[2][0]
									secondary_was_set = True
								else:
									secondary_allele = primary_allele.copy()
									secondary_allele['DiffConfuse'] = True
									secondary_was_set = True
							else:
								secondary_allele = primary_allele.copy()
								secondary_was_set = True
						## drop is small and alpha.CAG-beta.CAG == 1, thus slippage (theta real)
						elif 0.25 <= alpha_beta_ReadPcnt <= 0.49:
							secondary_allele = sorted_info[2][1]
							secondary_allele['Reference'] = sorted_info[2][0]
							secondary_was_set = True
						else:
							secondary_allele = sorted_info[1][1]
							secondary_allele['Reference'] = sorted_info[1][0]
							secondary_allele['DiffConfuse'] = True
							secondary_was_set = True

				## theta is slippage of alpha, beta is legitimate
				if alpha_estCCG != beta_estCCG:
					if alpha_estCCG == theta_estCCG and alpha_theta_CAGDiff == 1:
						secondary_allele = sorted_info[1][1]
						secondary_allele['Reference'] = sorted_info[1][0]
						secondary_was_set = True
				## all alleles different CCG
				if alpha_estCCG != beta_estCCG and beta_estCCG != theta_estCCG:
					## however, cag values on alpha/theta are close
					## and CCG within 1, so probably misread by sequencing machine
					## thus beta is legitimate
					if np.isclose([theta_estCAG], [alpha_estCAG], atol=5):
						if np.isclose([alpha_estCCG],[theta_estCCG], atol=2):
							secondary_allele = sorted_info[1][1]
							secondary_allele['Reference'] = sorted_info[1][0]
							secondary_was_set = True

		##
		## For each of the alleles we've determined..
		## Get intervening lengths, create accurate genotype string
		atypical_count = 0
		for allele in [primary_allele, secondary_allele]:
			if allele['Status'] == 'Atypical': atypical_count += 1
			else: allele['InterveningSequence'] = 'CAACAGCCGCCA'
			new_genotype, caacag_count, ccgcca_count = self.create_genotype_label(allele)
			allele['OriginalReference'] = allele['Reference']
			allele['Reference'] = new_genotype
			allele['EstimatedCAACAG'] = caacag_count
			allele['EstimatedCCGCCA'] = ccgcca_count

		##
		## Check for atypical allele rewriting CCG Het to CCG Hom
		temp_zyg = []; temp_curr = []
		for allele in [primary_allele, secondary_allele]:
			## diff confusion check
			try:
				if allele['DiffConfuse']: self.sequencepair_object.set_differential_confusion(True)
			except KeyError:
				allele['DiffConfuse'] = False
			orig_ccg = allele['OriginalReference'].split('_')[3]
			curr_ccg = allele['EstimatedCCG']
			## if original ref sect isn't string
			if not type(orig_ccg) == int: orig_ccg = int(filter(str.isdigit, orig_ccg))
			temp_zyg.append(orig_ccg); temp_curr.append(curr_ccg)
			if allele['Status'] == 'Atypical':
				if int(orig_ccg) != int(curr_ccg):
					self.sequencepair_object.set_atypical_ccgrewrite(True)
					allele['RewrittenCCG'] = orig_ccg
				if int(orig_ccg) == int(curr_ccg):
					allele['UnrewrittenCCG'] = orig_ccg
					self.sequencepair_object.set_atypical_zygrewrite(True)

		if not temp_zyg[0] == temp_zyg[1]:
			if temp_curr[0] == temp_curr[1]:
				if self.sequencepair_object.get_atypical_ccgrewrite():
					self.sequencepair_object.set_atypical_zygrewrite(True)
		if temp_zyg[0] == temp_zyg[1]:
			if temp_curr[0] != temp_zyg[0] or temp_curr[1] != temp_zyg[1]:
				self.sequencepair_object.set_atypical_zygrewrite(True)

		## potential homozygous haplotype
		if primary_allele['Reference'] == secondary_allele['Reference']:
			self.sequencepair_object.set_homozygoushaplotype(True)

		self.sequencepair_object.set_heuristicfilter(secondary_was_set)
		return primary_allele, secondary_allele, atypical_count

	def create_genotype_label(self, input_reference):

		##
		## Check before typical intervening for 'new' atypicals
		## Set up data structure
		## Check for rotations of known structures...
		intervening = input_reference['InterveningSequence']
		for real in ['CAACAG','CCGCCA','CAACAGCAACAGCCGCCA','CAACAGCCGCCACCGCCA']:
			if rotation_check(real, intervening):
				intervening = real
				input_reference['InterveningSequence'] = real

		##
		## Dictionary constructs to be filled by the following manual checks
		int_one = {'Mask': 'CAACAG', 'Count': 0, 'StartIDX': 0, 'EndIDX': 0, 'Label': '', 'Suffix': ''}
		int_two = {'Mask': 'CCGCCA', 'Count': 0, 'StartIDX': 0, 'EndIDX': 0, 'Label': '', 'Suffix': ''}

		##
		## Scrape data into structure
		for mask_dict in [int_one, int_two]:
			self.scraper(mask_dict, intervening)

		intervening_flag = True; atypical_flag = True; int_one_offset_flag = False
		int_two_offset_flag = False; int_one_investigate = False; int_two_investigate = False
		caacag_count = int_one['Count']; ccgcca_count = int_two['Count']
		int_one_offset = 0; int_one_simscore = 0; int_two_simscore = 0

		##########################
		##CAACAG (int one) check##
		##########################
		if int_one['Count'] > 0:
			if not int_one['StartIDX'] == 0:
				offset_str = intervening.split(int_one['Mask'])[0]
				int_one_offset = len(offset_str)
				int_one_offset_flag = True
				int_one_investigate = True
		else:
			remainder = len(intervening) % 6
			if not remainder == 0:
				offset_mutated = intervening.split(intervening[remainder:remainder + 6])[0]
				int_one_offset = len(offset_mutated)
				potential_mask = intervening[remainder:remainder + 6]
				int_one_offset_simscore = similar('CAACAG', potential_mask)
				int_one_offset_flag = True
				if int_one_offset_simscore >= 0.5:
					int_one['Mask'] = potential_mask
					self.scraper(int_one, intervening)
					int_one_investigate = True
			else:
				if len(intervening) > 6 and not int_one_offset_flag:
					int_one_simscore = similar('CAACAG', intervening[0:6])
					if int_one_simscore >= 0.5:
						int_one['Mask'] = intervening[0:6]
						self.scraper(int_one, intervening)
						int_one_investigate = True
				elif len(intervening) == 6:
					for sub_mask in ['CAA','CAG']:
						remainder = ''
						if sub_mask in intervening:
							remainder = intervening.replace(sub_mask, '')
						if rotation_check(remainder, sub_mask):
							int_one['Count'] = 1
							int_one['EndIDX'] = 6
							caacag_count = 1

		##########################
		##CCGCCA (int two) check##
		##########################
		if int_two['Count'] > 0:
			offset = (int_one['Count'] * 6) + int_one_offset
			if not int_two['StartIDX'] == offset:
				int_two_offset_flag = True
				offset_str = intervening.split(int_two['Mask'])[0].split(int_one['Mask'])[1]
				int_two_investigate = True
		else:
			remainder = len(intervening) % 6
			if not remainder == 0:
				lhinge = remainder + 6; rhinge = lhinge + 6
				offset_mutated = intervening.split(intervening[lhinge:rhinge])[0].split(int_one['Mask'])[1]
				offset = (int_one['Count'] * 6) + len(offset_mutated)
				if not int_two['StartIDX'] == offset:
					offset_str = intervening[lhinge:rhinge]
					int_two_offset_simscore = similar('CCGCCA', offset_str)
					int_two_offset_flag = True
					if int_two_offset_simscore >= 0.5:
						int_two['Mask'] = offset_str
						self.scraper(int_two, intervening)
						int_two_investigate = True
			if len(intervening) > 6 and not int_two_offset_flag:
				int_two_simscore = similar('CCGCCA', intervening[6:12])
				if int_two_simscore >= 0.5:
					int_two['Mask'] = intervening[6:12]
					self.scraper(int_two, intervening)
					int_two_investigate = True
			elif len(intervening) == 6:
				for sub_mask in ['CCG','CCA']:
					remainder = ''
					if sub_mask in intervening:
						remainder = intervening.replace(sub_mask, '')
					if rotation_check(remainder, sub_mask):
						int_two['Count'] = 1
						int_two['StartIDX'] = int_one['EndIDX'] + 6
						int_two['EndIDX'] = int_two['StartIDX'] + 6
						ccgcca_count = 1

		##################
		## Header check ##
		##################
		if caacag_count == 0 and ccgcca_count == 0: intervening_flag = False
		if caacag_count != 1 and ccgcca_count != 1: atypical_flag = True

		#################################
		##Anything longer than typical?##
		#################################
		if intervening_flag and len(intervening) > 12:
			returned_suffix = intervening[int_two['EndIDX']:]
			if returned_suffix:
				int_two_investigate = True

		###########################
		##If not present at all..##
		###########################
		if not intervening_flag:
			int_one['Count'] = 0; int_one_investigate = True
			int_two_investigate = True; int_two['Count'] = 0

		###############
		##Easy checks##
		###############
		for count in [int_one['Count'], int_two['Count']]:
			if not int(count) in [0,1,2]:
				self.sequencepair_object.set_novel_atypical_structure(True)
		if not int(input_reference['EstimatedCCT']) in [0,1,2,3]:
			self.sequencepair_object.set_novel_atypical_structure(True)
		if input_reference['Status'] == 'Typical' and (caacag_count != 1 or ccgcca_count != 1):
			caacag_count = 1; ccgcca_count = 1

		########################
		##Build genotype label##
		########################
		if int_one_investigate:
			int_one['Suffix'] = '*'
			self.sequencepair_object.set_novel_atypical_structure(True)
		if int_two_investigate:
			int_two['Suffix'] = '*'
			self.sequencepair_object.set_novel_atypical_structure(True)


		int_one['Label'] = '{}'.format(str(int_one['Count']) + int_one['Suffix'])
		int_two['Label'] = '{}'.format(str(int_two['Count']) + int_two['Suffix'])
		genotype_label = '{}_{}_{}_{}_{}'.format(input_reference['EstimatedCAG'], int_one['Label'], int_two['Label'],
												 input_reference['EstimatedCCG'], input_reference['EstimatedCCT'])

		return genotype_label, caacag_count, ccgcca_count

	def get_atypicalreport(self):
		return self.atypical_report

	@staticmethod
	def scraper(intv_dict, intervening_str):
		for dna_module in re.finditer(intv_dict['Mask'], intervening_str):
			if intv_dict['Count'] == 0:
				intv_dict['StartIDX'] = dna_module.start(); intv_dict['EndIDX'] = dna_module.end()
			elif intv_dict['Count'] != 0:
				intv_dict['EndIDX'] = dna_module.end()
			intv_dict['Count'] += 1
		return intv_dict
