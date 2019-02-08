from __future__ import division

#/usr/bin/python
__version__ = 0.321
__author__ = 'alastair.maxwell@glasgow.ac.uk'

##
## Generic imports
import os
import csv
import PyPDF2
import warnings
import peakutils
import matplotlib
import collections
import numpy as np
import scipy as sp
matplotlib.use('Agg')
import logging as log
import seaborn as sns
from sklearn import svm
import scipy.stats as st
import matplotlib.pyplot as plt
from sklearn import preprocessing
from reportlab.pdfgen import canvas
from peakutils.plot import plot as pplot
from sklearn.multiclass import OutputCodeClassifier

##
## Backend Junk
from ..__backend import DataLoader
from ..__backend import Colour as clr

def split_cag_target(input_distribution):
	"""
	Function to gather the relevant CAG distribution for the specified CCG value
	We gather this information from the forward distribution of this sample pair as CCG reads are
	of higher quality in the forward sequencing direction.
	We split the entire fw_dist into contigs/bins for each CCG (4000 -> 200*20)
	:param input_distribution: input forward distribution (4000d)
	:param ccg_target: target value we want to select the 200 values for
	:return: the sliced CAG distribution for our specified CCG value
	"""

	cag_split = [input_distribution[i:i + 200] for i in xrange(0, len(input_distribution), 200)]
	distribution_dict = collections.OrderedDict()
	for i in range(0, len(cag_split)):
		distribution_dict['CCG' + str(i + 1)] = cag_split[i]

	# current_target_distribution = distribution_dict['CCG' + str(ccg_target)]
	return distribution_dict

class AlleleGenotyping:
	def __init__(self, sequencepair_object, instance_params, training_data, atypical_logic=None, padded_target=None):

		##
		## Allele objects and instance data
		self.sequencepair_object = sequencepair_object
		self.instance_params = instance_params
		self.training_data = training_data
		self.invalid_data = atypical_logic
		self.padded_target = padded_target
		self.allele_report = ''
		self.warning_triggered = False

		##
		## Constructs that will be updated with each allele process
		self.classifier, self.encoder = self.build_zygosity_model()
		self.allele_flags = {}; self.forward_distribution = None; self.reverse_distribution = None
		self.primary_original = None; self.secondary_original = None
		self.forward_aggregate = None; self.reverse_aggregate = None
		self.expected_zygstate = None; self.zygosity_state = None
		self.pass_vld = True; self.ccg_sum = []

		##
		## Genotype!
		if not self.allele_validation(): raise Exception('Allele(s) failed validation. Cannot genotype..')
		if not self.determine_ccg(): raise Exception('CCG Genotyping failure. Cannot genotype..')
		if not self.determine_cag(): raise Exception('CAG Genotyping failure. Cannot genotype..')
		if not self.genotype_validation(): raise Exception('Genotype failed validation. Cannot genotype..')
		if not self.inspect_peaks():
			log.warn('{}{}{}{}'.format(clr.red, 'shd__ ', clr.end, '1+ allele(s) failed peak validation. Precision not guaranteed.'))
			self.warning_triggered = True
			self.sequencepair_object.set_peakinspection_warning(True)
		self.n_align_dist()
		self.calculate_score()
		self.contextualise()
		self.render_graphs()
		self.set_report()

	def build_zygosity_model(self):
		"""
		Function to build a SVM (wrapped into OvO class) for determining CCG zygosity
		:return: svm object wrapped into OvO, class-label hash-encoder object
		"""

		##
		## Classifier object and relevant parameters for our CCG prediction
		svc_object = svm.LinearSVC(C=1.0, loss='ovr', penalty='l2', dual=False,
								   tol=1e-4, multi_class='crammer_singer', fit_intercept=True,
								   intercept_scaling=1, verbose=0, random_state=0, max_iter=100000)

		##
		## Take raw training data (CCG zygosity data) into DataLoader model object
		traindat_ccg_collapsed = self.training_data['CollapsedCCGZygosity']
		traindat_descriptionfi = self.training_data['GenericDescriptor']
		traindat_model = DataLoader(traindat_ccg_collapsed, traindat_descriptionfi).load_model()

		##
		## Model data fitting to SVM
		X = preprocessing.normalize(traindat_model.DATA)
		Y = traindat_model.TARGET
		ovo_svc = OutputCodeClassifier(svc_object, code_size=2, random_state=0).fit(X, Y)
		encoder = traindat_model.ENCDR

		##
		## Return the fitted OvO(SVM) and Encoder
		return ovo_svc, encoder

	def predict_zygstate(self):
		"""
		Function which takes the newly collapsed CCG distribution and executes SVM prediction
		to determine the zygosity state of this sample's CCG value(s). Data is reshaped
		and normalised to ensure more reliable results. A check is executed between the results of
		forward and reverse zygosity; if a match, great; if not, not explicitly bad but inform user.
		:return: zygosity[2:-2] (trimming unrequired characters)
		"""

		##
		## Reshape the input distribution so SKL doesn't complain about 1D vectors
		## Normalise data in addition; cast to float64 for this to be permitted
		forward_reshape = preprocessing.normalize(np.float64(self.forward_aggregate.reshape(1, -1)))
		reverse_reshape = preprocessing.normalize(np.float64(self.reverse_aggregate.reshape(1, -1)))

		##
		## Predict the zygstate of these reshapen, normalised 20D CCG arrays using SVM object earlier
		## Results from self.classifier are #encoded; so convert with our self.encoder.inverse_transform
		forward_zygstate = str(self.encoder.inverse_transform(self.classifier.predict(forward_reshape)))
		reverse_zygstate = str(self.encoder.inverse_transform(self.classifier.predict(reverse_reshape)))

		##
		## We only particularly care about the reverse zygosity (CCG reads are higher quality in reverse data)
		## However, for a QoL metric, compare fw/rv results. If match, good! If not, who cares!
		if not forward_zygstate == reverse_zygstate:
			self.allele_flags['CCGZygDisconnect'] = True
		else:
			self.allele_flags['CCGZyg_disconnect'] = False

		return reverse_zygstate[2:-2]

	def index_inspector(self, index_inspection_count):

		major = max(self.reverse_aggregate)
		majoridx = np.where(self.reverse_aggregate == major)[0][0]
		minor = max(n for n in self.reverse_aggregate if n != major)
		minoridx = np.where(self.reverse_aggregate == minor)[0][0]

		if index_inspection_count == 2:
			return [(major, majoridx), (minor, minoridx)]
		if index_inspection_count == 1:
			return [(major, majoridx)]

	@staticmethod
	def scrape_distro(distributionfi):
		"""
		Function to take the aligned read-count distribution from CSV into a numpy array
		:param distributionfi:
		:return: np.array(data_from_csv_file)
		"""

		##
		## Open CSV file with information within; append to temp list
		## Scrape information, cast to np.array(), return
		placeholder_array = []
		with open(distributionfi) as dfi:
			source = csv.reader(dfi, delimiter=',')
			next(source)  # skip header
			for row in source:
				placeholder_array.append(int(row[2]))
			dfi.close()
		unlabelled_distro = np.array(placeholder_array)
		return unlabelled_distro

	@staticmethod
	def distribution_collapse(distribution_array):
		"""
		Function to take a full 200x20 array (struc: CAG1-200,CCG1 -- CAG1-200CCG2 -- etc CCG20)
		and aggregate all CAG values for each CCG
		:param distribution_array: input dist (should be (1-200,1-20))
		:return: 1x20D np(array)
		"""

		##
		## Hopefully the user has aligned to the right reference dimensions
		try:
			ccg_arrays = np.split(distribution_array, 20)
		except ValueError:
			raise Exception('Input reads individisible by 20. Utilised incorrect reference style.')

		##
		## Aggregate each CCG
		ccg_counter = 1
		collapsed_array = []
		for ccg_array in ccg_arrays:
			collapsed_array.append(np.sum(ccg_array))
			ccg_counter += 1

		return np.asarray(collapsed_array)

	@staticmethod
	def pad_distribution(distribution_array, allele_object):

		local_index = np.where(distribution_array == max(distribution_array))[0][0]
		local_rightpad = len(distribution_array) - local_index
		global_index = allele_object.get_ccg() - 1
		left_buffer = abs(local_index - global_index)
		right_buffer = abs(20 - global_index) - local_rightpad
		left_pad = np.asarray([0] * left_buffer)
		right_pad = np.asarray([0] * right_buffer)
		left_aug = np.concatenate((left_pad, distribution_array))
		right_aug = np.concatenate((left_aug, right_pad))

		return right_aug

	@staticmethod
	def split_cag_target(input_distribution):
		"""
		Function to gather the relevant CAG distribution for the specified CCG value
		We gather this information from the forward distribution of this sample pair as CCG reads are
		of higher quality in the forward sequencing direction.
		We split the entire fw_dist into contigs/bins for each CCG (4000 -> 200*20)
		:param input_distribution: input forward distribution (4000d)
		:param ccg_target: target value we want to select the 200 values for
		:return: the sliced CAG distribution for our specified CCG value
		"""

		cag_split = [input_distribution[i:i + 200] for i in xrange(0, len(input_distribution), 200)]
		distribution_dict = collections.OrderedDict()
		for i in range(0, len(cag_split)):
			distribution_dict['CCG' + str(i + 1)] = cag_split[i]

		# current_target_distribution = distribution_dict['CCG' + str(ccg_target)]
		return distribution_dict

	def close_check(self, allele, array, x, y, z, state=None):
		inner_pass = True
		if np.isclose(array, x, atol=y):
			if state == 'minus':
				allele.set_nminuswarninglevel(z)
				if z >= 5: inner_pass = False
			if state == 'plus':
				allele.set_npluswarninglevel(z)
				if z >= 5: inner_pass = False
		else:
			allele.set_nminuswarninglevel(0)
			allele.set_npluswarninglevel(0)
		self.pass_vld = inner_pass

	def peak_detection(self, allele_object, distro, peak_dist, triplet_stage, est_dist=None, fod_recall=False):

		##
		## Status
		fail_state = False
		utilised_threshold = 0.50
		error_boundary = 0

		##
		## If we're in a re-call situation, lower peak threshold
		## Otherwise, threshold already assigned to object is utilised
		if fod_recall:

			recall_count = self.sequencepair_object.get_recallcount()
			self.sequencepair_object.set_recallcount(recall_count+1)
			if recall_count > 7: raise Exception('7+ recalls. Unable to determine genotype.')
			threshold = 0.0
			if triplet_stage == 'CCG': threshold = allele_object.get_ccgthreshold()
			if triplet_stage in ['CAG', 'CAGHet', 'CAGHom', 'CAGDim']: threshold = allele_object.get_cagthreshold()
			threshold -= 0.06
			utilised_threshold = max(threshold, 0.05)

		##
		## CCG, CAGHet (Hetero, Homo*, Homo+), CAGDim == 1 peak
		## CAGHom == 2 peaks
		allele_object.set_cagthreshold(utilised_threshold)
		if triplet_stage == 'CAGHom':
			error_boundary = 2
		elif triplet_stage == 'CCG':
			if self.zygosity_state == 'HETERO':
				error_boundary = 2
			if self.zygosity_state == 'HOMO':
				error_boundary = 1
		else: error_boundary = 1

		##
		## Look for peaks in our distribution
		peak_indexes = peakutils.indexes(distro, thres=utilised_threshold, min_dist=peak_dist)
		fixed_indexes = np.array(peak_indexes + 1)
		if not len(fixed_indexes) == error_boundary:
			if triplet_stage == 'CAGHom' and (est_dist==1 or est_dist==0):
				pass
			elif allele_object.get_cag() in fixed_indexes:
				fixed_indexes = np.asarray([x for x in fixed_indexes if x == allele_object.get_cag()])
			elif triplet_stage == 'CCG':
				if len(fixed_indexes) > 2 and not self.zygosity_state == 'HETERO':
					fixed_indexes = [np.where(distro == max(distro))[0][0]+1]
					if self.zygosity_state == 'HOMO+' or self.zygosity_state == 'HOMO*':
						self.sequencepair_object.set_svm_failure(False)
						pass
					else:
						self.sequencepair_object.set_svm_failure(True)
						self.sequencepair_object.set_alignmentwarning(True)
			else:
				fail_state = True

		return fail_state, fixed_indexes

	def allele_validation(self):

		ccg_expectant = []
		##
		## For the two allele objects in this sample_pair
		for allele_object in [self.sequencepair_object.get_primaryallele(),
							  self.sequencepair_object.get_secondaryallele()]:

			##
			## Assign read mapped percent if not present in allele
			if not allele_object.get_fwalnpcnt() and not allele_object.get_rvalnpcnt():
				allele_object.set_fwalnpcnt(self.sequencepair_object.get_fwalnpcnt())
				allele_object.set_fwalncount(self.sequencepair_object.get_fwalncount())
				allele_object.set_fwalnrmvd(self.sequencepair_object.get_fwalnrmvd())
				allele_object.set_rvalnpcnt(self.sequencepair_object.get_rvalnpcnt())
				allele_object.set_rvalncount(self.sequencepair_object.get_rvalncount())
				allele_object.set_rvalnrmvd(self.sequencepair_object.get_rvalnrmvd())

			##
			## Unlabelled distributions
			self.forward_distribution = self.scrape_distro(allele_object.get_fwdist())
			self.reverse_distribution = self.scrape_distro(allele_object.get_rvdist())
			allele_object.set_fwarray(self.forward_distribution)
			allele_object.set_rvarray(self.reverse_distribution)

			##
			## Distribution ead count / Peak read count
			if allele_object.get_totalreads() < 750:
				allele_object.set_distribution_readcount_warning(True)
				self.sequencepair_object.set_distribution_readcount_warning(True)

			##
			## If current alleleobj's assembly/distro is blank
			## Allele is typical, didn't assign values in __atypical.py
			## Hence, set these values here (From seqpair object, where they reside)
			if not allele_object.get_rvdist():
				allele_object.set_fwassembly(self.sequencepair_object.get_fwassembly())
				allele_object.set_rvassembly(self.sequencepair_object.get_rvassembly())
				allele_object.set_fwdist(self.sequencepair_object.get_fwdist())
				allele_object.set_rvdist(self.sequencepair_object.get_rvdist())

			###############################
			## Stage one -- CCG Zygosity ##
			###############################
			## Typical allele
			if allele_object.get_allelestatus() == 'Typical':
				self.forward_aggregate = self.distribution_collapse(self.forward_distribution)
				self.reverse_aggregate = self.distribution_collapse(self.reverse_distribution)

			## Atypical allele
			if allele_object.get_allelestatus() == 'Atypical':
				## Data has been realigned to custom reference
				if not self.invalid_data:
					self.forward_aggregate = self.distribution_collapse(self.forward_distribution)
					self.reverse_aggregate = self.reverse_distribution
					allele_object.set_rvarray(self.reverse_aggregate)
				## Data has not been realigned -- brute force genotyping
				if self.invalid_data:
					self.forward_aggregate = self.distribution_collapse(self.forward_distribution)
					self.reverse_aggregate = self.distribution_collapse(self.reverse_distribution)
			self.zygosity_state = self.predict_zygstate()

			##
			## Allele read peak (dsp error)
			major = max(self.reverse_aggregate)
			minor = max(n for n in self.reverse_aggregate if n != major)
			for item in [major, minor]:
				if item == 0 or item is None:
					raise ValueError('Insignificant read count. Cannot genotype.')
			tertiary = max(n for n in self.reverse_aggregate if n!= major and n!= minor)
			majidx = int(np.where(self.reverse_aggregate == major)[0][0])
			minidx = int(np.where(self.reverse_aggregate == minor)[0][0])
			## Check discrete diff between majidx and minidx
			abs_ratio = self.reverse_aggregate[minidx] / self.reverse_aggregate[majidx]

			## check if minor is n-1/n+1 AND third peak is np.close(minor)
			skip_flag = False
			if not self.sequencepair_object.get_primaryallele().get_neighbouring_candidate():
				if not self.sequencepair_object.get_secondaryallele().get_neighbouring_candidate():
					if abs(majidx-minidx) == 1:
						if np.isclose([minor], [tertiary], atol=minor*0.80):
							skip_flag = True
							allele_object.set_ccguncertainty(True)
							self.sequencepair_object.set_ccguncertainty(True)
					if minor == self.reverse_aggregate[allele_object.get_ccg()-1]:
						skip_flag = True
						allele_object.set_ccguncertainty(True)
						self.sequencepair_object.set_ccguncertainty(True)
					## skip the following block if so
					if not skip_flag:
						if abs_ratio < 0.05: pass
						else:
							for fod_peak in [majidx+1, minidx+1]:
								if allele_object.get_ccg() not in [majidx+1, minidx]:
									if fod_peak not in [self.sequencepair_object.get_primaryallele().get_ccg(),
														self.sequencepair_object.get_secondaryallele().get_ccg()]:
										allele_object.set_fodoverwrite(True)
										allele_object.set_ccgval(int(fod_peak))

			##
			## Clean up distribution for erroneous peaks
			## In the case of atypical alleles, unexpected peaks may exist in aggregates
			if self.zygosity_state == 'HOMO' or self.zygosity_state == 'HOMO*' or self.zygosity_state == 'HOMO+':
				for i in range(0, len(self.reverse_aggregate)):
					if np.isclose([i], [allele_object.get_ccg() - 1], atol=3):
						if np.isclose([self.reverse_aggregate[i]], [self.reverse_aggregate[allele_object.get_ccg() - 1]],
									  atol=((self.reverse_aggregate[allele_object.get_ccg() - 1])/100) * 45):
							removal = (self.reverse_aggregate[i]/100) * 77.5
							if i != allele_object.get_ccg()-1:
								self.reverse_aggregate[i] -= removal
						else:
							removal = (self.reverse_aggregate[i] / 100) * 77.5
							if i != allele_object.get_ccg() - 1:
								self.reverse_aggregate[i] -= removal
			else:
				## if the current allele is top2, check difference between top1/top2
				## if the difference is large, we need to further smooth the distribution for this allele
				top1 = max(self.reverse_aggregate); top1idx = list(self.reverse_aggregate).index(top1)
				top2 = max([x for x in self.reverse_aggregate if x != top1])
				additional_context = 0
				if top2 == self.reverse_aggregate[allele_object.get_ccg()-1]:
					differential = (top2/top1)*100
					if 0 < differential < 50:
						purge = (self.reverse_aggregate[top1idx]/100)*95
						self.reverse_aggregate[top1idx] -= purge
					else:
						additional_context = 5
				## the percentage we should remove errorneous reads by
				## should differ based on the context of n's read count
				if 0 < self.reverse_aggregate[allele_object.get_ccg()-1] <= 6000:
					removal_context = 75
				elif 6000 <= self.reverse_aggregate[allele_object.get_ccg()-1] <= 12000:
					removal_context = 85
				else:
					removal_context = 95

				## actual cleanup stage
				for i in range(0, len(self.reverse_aggregate)):
					if i != allele_object.get_ccg()-1:
						removal = (self.reverse_aggregate[i]/100) * removal_context+additional_context
						self.reverse_aggregate[i] -= removal
						if self.reverse_aggregate[i] < 0: self.reverse_aggregate[i] = 0


			##
			## Check SVM didn't fail...
			## (refresh variables for specific distribution)
			indv_major = max(self.reverse_aggregate)
			indv_minor = max(n for n in self.reverse_aggregate if n != major)
			indv_majidx = int(np.where(self.reverse_aggregate == indv_major)[0][0])
			indv_minidx = int(np.where(self.reverse_aggregate == indv_minor)[0][0])

			if 1 < abs(indv_majidx-indv_minidx) < 10:
				peak_count = 0
				if abs_ratio < 0.05:
					pass
				## otherwise, perhaps SVM was wrong and missed a peak
				else:
					for peak in [indv_majidx, indv_minidx]:
						pmo = self.reverse_aggregate[peak-1]; ppo = self.reverse_aggregate[peak+1]
						pmo_ratio = pmo/self.reverse_aggregate[peak]
						ppo_ratio = ppo/self.reverse_aggregate[peak]
						## calc ratio of peak+1/-1, if within range, add peak

						if pmo_ratio and ppo_ratio < 0.2:
							if 0.0 not in [pmo_ratio, ppo_ratio]:
								peak_count += 1

					## hotfix SVM results based on peak detection
					if peak_count == 2 and not self.zygosity_state == 'HETERO':
						if self.zygosity_state == 'HOMO+' or self.zygosity_state == 'HOMO*':
							self.sequencepair_object.set_svm_failure(False)
							pass
						else:
							## ratios dictate that peak may be legit
							## check between major and minor for validity
							suspect_ratio = indv_minor/indv_major
							if suspect_ratio < 0.15:
								pass
							else:

								self.zygosity_state = 'HETERO'
								self.sequencepair_object.set_svm_failure(True)

			## set distribution
			allele_object.set_rvarray(self.reverse_aggregate)

			#################################
			## Stage two -- CCG continuity ##
			#################################
			index_inspection_count = 0
			if self.zygosity_state == 'HETERO': index_inspection_count = 2
			if self.zygosity_state == 'HOMO': index_inspection_count = 1
			inspections = self.index_inspector(index_inspection_count)
			for inspect in inspections:
				if np.isclose(allele_object.get_ccg(), [inspect[1]+1], atol=1):
					allele_object.set_validation(True)
			ccg_expectant.append(allele_object.get_ccg())

		try:
			if not ccg_expectant[0] == ccg_expectant[1]:
				self.expected_zygstate = 'HETERO'
			if ccg_expectant[0] == ccg_expectant[1]:
				self.expected_zygstate = 'HOMO'
		except IndexError:
			raise Exception('CCG Prediction Failure.')

		##
		## If atypical detected, but zygosity was rewritten
		## allele CCG remained the same, but one allele is now atypical
		if not self.sequencepair_object.get_atypical_ccgrewrite():
			if self.sequencepair_object.get_atypical_zygrewrite():
				self.zygosity_state = 'HOMO+'
		## allele CCG value was changed as a result of the atypical detection
		else:
			if self.sequencepair_object.get_atypical_zygrewrite():
				self.zygosity_state = 'HOMO*'

		##
		## Check both alleles passed validation
		if (self.sequencepair_object.get_primaryallele().get_validation()) and (
				self.sequencepair_object.get_secondaryallele().get_validation()):
			return True
		else:
			return False

	def determine_ccg(self):

		##
		## Constructs
		ccg_matches = 0; ccg_values = []; local_zygstate = None; pass_gtp = True; ccg_sum = []

		##
		## For the two allele objects in this sample_pair
		## First, ensure CCG matches between DSP estimate and FOD derision
		for allele in [self.sequencepair_object.get_primaryallele(), self.sequencepair_object.get_secondaryallele()]:

			allele.set_ccgthreshold(0.50)
			fod_failstate, ccg_indexes = self.peak_detection(allele, allele.get_rvarray(), 1, 'CCG')
			while fod_failstate:
				fod_failstate, ccg_indexes = self.peak_detection(allele, allele.get_rvarray(), 1, 'CCG', fod_recall=True)
			if len(ccg_indexes) == 0:
				raise Exception('CCG Peak un-callable; cannot process.')
			if allele.get_ccg() in ccg_indexes:
				ccg_matches += 1
				allele.set_ccgvalid(True)
			ccg_values.append([x for x in ccg_indexes if x == allele.get_ccg()])

			if len(ccg_indexes) > 1:
				if allele.get_header() == 'PRI':
					allele.set_fodccg(np.asarray(ccg_indexes[0]))
				if allele.get_header() == 'SEC':
					allele.set_fodccg(np.asarray(ccg_indexes[1]))
			else:
				allele.set_fodccg(np.asarray(ccg_indexes[0]))

			distribution_split = split_cag_target(allele.get_fwarray())
			target_distro = distribution_split['CCG{}'.format(allele.get_ccg())]
			ccg_sum.append([allele.get_ccg(), sum(target_distro)])

		if ccg_values[0] == ccg_values[1]:
			local_zygstate = 'HOMO'
		if not ccg_values[0] == ccg_values[1]:
			local_zygstate = 'HETERO'

		##
		## If the sample's total read count is so low that we cannot trust results
		## We trust the local/expected zygosity over the SVM derived instance-wide variable
		if self.sequencepair_object.get_alignmentwarning():
			if sum(self.reverse_aggregate) < 100:
				self.zygosity_state = self.expected_zygstate = local_zygstate

		self.sequencepair_object.set_ccgzygstate(self.expected_zygstate)
		if not self.zygosity_state == 'HOMO*' or not self.zygosity_state == 'HOMO+':
			if not local_zygstate == self.expected_zygstate:
				if abs(ccg_sum[0][0]-ccg_sum[1][0]) == 1:
					if not np.isclose([ccg_sum[0][1]],[ccg_sum[1][1]],atol=(0.70*max(ccg_sum)[1])):
						pass_gtp = True
						self.ccg_sum = ccg_sum
						pass
				else:
					pass_gtp = False

		return pass_gtp

	def determine_cag(self):

		##
		## Constructs
		pass_gtp = True
		###############################################
		## Pre-Check: atypical allele mis-assignment ##
		###############################################
		pri_ccg = self.sequencepair_object.get_primaryallele().get_ccg()
		sec_ccg = self.sequencepair_object.get_secondaryallele().get_ccg()
		if self.sequencepair_object.get_atypicalcount() > 0:
			if self.zygosity_state == 'HOMO':
				if pri_ccg != sec_ccg:
					self.zygosity_state = 'HETERO'
					self.sequencepair_object.set_ccgzygstate(self.zygosity_state)
			if self.zygosity_state == 'HETERO':
				if pri_ccg == sec_ccg:
					self.zygosity_state = 'HOMO'
					self.sequencepair_object.set_ccgzygstate(self.zygosity_state)
			if self.zygosity_state == 'HOMO*':
				pass

		##
		## Assign distro originals
		primary_dist = self.sequencepair_object.get_primaryallele().get_fwarray().copy()
		primary_split = split_cag_target(primary_dist)
		self.primary_original = primary_split['CCG{}'.format(self.sequencepair_object.get_primaryallele().get_ccg())]
		secondary_dist = self.sequencepair_object.get_secondaryallele().get_fwarray().copy()
		secondary_split = split_cag_target(secondary_dist)
		self.secondary_original = secondary_split['CCG{}'.format(self.sequencepair_object.get_secondaryallele().get_ccg())]

		##
		## If we have an atypical allele in this sample, the remaining typical allele distribution may be skewed
		## e.g. something aligning to CAG_1_1_7_2 would have aligned to CAG_1_0_9_2
		## where a typical distribution would originally be CCG homozygous, the CAG_1_0_9_2 reads are still
		## assigned within the typical distribution... clean these up so genotyping isn't misassociating data
		if self.sequencepair_object.get_atypicalcount() == 1:
			for allele in [self.sequencepair_object.get_primaryallele(), self.sequencepair_object.get_secondaryallele()]:
				if allele.get_allelestatus() == 'Typical':
					distribution_split = split_cag_target(allele.get_fwarray())
					target_distro = distribution_split['CCG{}'.format(allele.get_ccg())]
					for i in range(0, len(target_distro)):
						if i != allele.get_cag() - 1:
							removal = (target_distro[i] / 100) * 85
							target_distro[i] -= removal

		##########################
		## Heterozygous for CCG ##
		##########################
		if self.zygosity_state == 'HETERO' or self.zygosity_state == 'HOMO*' or self.zygosity_state == 'HOMO+':
			existing_calls = []
			for allele in [self.sequencepair_object.get_primaryallele(), self.sequencepair_object.get_secondaryallele()]:

				distribution_split = split_cag_target(allele.get_fwarray())
				target_distro = distribution_split['CCG{}'.format(allele.get_ccg())]
				allele.set_totalreads(sum(target_distro))

				if self.zygosity_state == 'HOMO+' or self.zygosity_state == 'HOMO*':
					for i in range(0, len(target_distro)):
						if i != allele.get_cag() - 1:
							removal = (target_distro[i] / 100) * 85
							target_distro[i] -= removal

				allele.set_cagthreshold(0.50)
				fod_failstate, cag_indexes = self.peak_detection(allele, target_distro, 1, 'CAGHet')
				while fod_failstate:
					fod_failstate, cag_indexes = self.peak_detection(allele, target_distro, 1, 'CAGHet', fod_recall=True)

				## check that FOD didn't return more items than it was required for this allele
				## only keep discrete values from the inferred total of all calls in the current sample
				if not self.sequencepair_object.get_homozygoushaplotype():
					for item in cag_indexes:
						gtype = (item, allele.get_ccg())
						## However we can get the case where CAG match between alleles but CCG doesn't
						## so we must add tuples instead of integers, and check the correct element against our observation
						if not any(gtype[1] in obs_tuple for obs_tuple in existing_calls):
							existing_calls.append(gtype)
							if type(cag_indexes) == np.ndarray:
								itemindex = np.where(cag_indexes == item)
								allele.set_fodcag(cag_indexes.flat[itemindex])
							else:
								allele.set_fodcag(cag_indexes)
						else:
							if type(cag_indexes) == np.ndarray:
								itemindex = np.where(cag_indexes == item)
								allele.set_fodcag(cag_indexes.flat[itemindex])
							else:
								allele.set_fodcag(cag_indexes)
				else:
					allele.set_fodcag(cag_indexes)

		########################
		## Homozygous for CCG ##
		########################
		if self.zygosity_state == 'HOMO':
			##
			## Double check CCG matches.. be paranoid
			primary_ccg = self.sequencepair_object.get_primaryallele().get_ccg()
			secondary_ccg = self.sequencepair_object.get_secondaryallele().get_ccg()
			try:
				if not primary_ccg == secondary_ccg:
					target = max(self.ccg_sum)[0]
					self.sequencepair_object.get_primaryallele().set_ccgval(target)
					self.sequencepair_object.get_secondaryallele().set_ccgval(target)
			except ValueError:
				max_array = [0, 0]
				for allele in [self.sequencepair_object.get_primaryallele(), self.sequencepair_object.get_secondaryallele()]:
					distro_split = split_cag_target(allele.get_fwarray())
					total_reads = sum(distro_split['CCG{}'.format(allele.get_ccg())])

					if total_reads > max_array[1]:
						max_array[1] = total_reads
						max_array[0] = allele.get_ccg()
				self.sequencepair_object.get_primaryallele().set_ccgval(max_array[0])
				self.sequencepair_object.get_secondaryallele().set_ccgval(max_array[0])

			##
			## Get distance estimate between two peaks in our target CCG distribution
			## set threshold to use in peak calling algorithm
			estimated_distance = abs(self.sequencepair_object.get_secondaryallele().get_cag() -
									 self.sequencepair_object.get_primaryallele().get_cag())

			if estimated_distance > 5: distance_threshold = 2
			elif estimated_distance == 1: distance_threshold = 0
			else: distance_threshold = 1

			##
			## Process each allele, getting the specific CCG distribution
			## Re-set read count for the allele, due to subsampling
			existing_calls = []
			for allele in [self.sequencepair_object.get_primaryallele(), self.sequencepair_object.get_secondaryallele()]:

				distribution_split = split_cag_target(allele.get_fwarray())
				target_distro = distribution_split['CCG{}'.format(allele.get_ccg())]

				if not self.sequencepair_object.get_primaryallele().get_neighbouring_candidate():
					if not self.sequencepair_object.get_secondaryallele().get_neighbouring_candidate():
						for i in range(0, len(target_distro)):
							if i != allele.get_cag() - 1:
								removal = (target_distro[i] / 100) * 85
								target_distro[i] -= removal

				allele.set_totalreads(sum(target_distro))
				allele.set_cagthreshold(0.50)

				fod_failstate, cag_indexes = self.peak_detection(allele, target_distro, distance_threshold, 'CAGHom', est_dist=estimated_distance)
				while fod_failstate:
					fod_failstate, cag_indexes = self.peak_detection(allele, target_distro, distance_threshold, 'CAGHom', est_dist=estimated_distance, fod_recall=True)

				## check that FOD didn't return more items than it was required for this allele
				## only keep discrete values from the inferred total of all calls in the current sample
				if not self.sequencepair_object.get_homozygoushaplotype():
					for item in cag_indexes:
						if not item in existing_calls:
							existing_calls.append(item)
							if type(cag_indexes) == np.ndarray:
								itemindex = np.where(cag_indexes == item)
								allele.set_fodcag(cag_indexes.flat[itemindex])
							else:
								allele.set_fodcag(cag_indexes)
						else:
							if type(cag_indexes) == np.ndarray:
								itemindex = np.where(cag_indexes == item)
								allele.set_fodcag(cag_indexes.flat[itemindex])
							else:
								allele.set_fodcag(cag_indexes)
				else:
					allele.set_fodcag(cag_indexes)
		return pass_gtp

	def genotype_validation(self):

		##
		## Constructs
		pass_vld = True
		primary_allele = self.sequencepair_object.get_primaryallele()
		secondary_allele = self.sequencepair_object.get_secondaryallele()
		pri_distro_split = split_cag_target(primary_allele.get_fwarray())
		sec_distro_split = split_cag_target(secondary_allele.get_fwarray())
		ccg_zygstate = self.zygosity_state

		##
		## Primary Allele
		primary_dsp_ccg = primary_allele.get_ccg(); primary_fod_ccg = primary_allele.get_fodccg()
		primary_dsp_cag = primary_allele.get_cag(); primary_fod_cag = primary_allele.get_fodcag()
		primary_peakreads = (split_cag_target(primary_allele.get_fwarray())['CCG{}'.format(primary_dsp_ccg)])[
			primary_dsp_cag-1]
		primary_allele.set_peakreads(primary_peakreads)

		##
		## Secondary Allele
		secondary_dsp_ccg = secondary_allele.get_ccg(); secondary_fod_ccg = secondary_allele.get_fodccg()
		secondary_dsp_cag = secondary_allele.get_cag(); secondary_fod_cag = secondary_allele.get_fodcag()
		secondary_peakreads = (split_cag_target(secondary_allele.get_fwarray())['CCG{}'.format(secondary_dsp_ccg)])[
			secondary_dsp_cag - 1]
		secondary_allele.set_peakreads(secondary_peakreads)

		##
		## Double check fod peaks
		def dimension_checker(input_list):

			## data
			fod = input_list[0]; dsp = input_list[1];allele = input_list[2]

			## casting
			if type(fod) is np.ndarray:
				fod = input_list[0].tolist()
			elif type(fod) is not list:
				fod = [input_list[0]]

			## validity between DSP/FOD
			for i in range(0, len(fod)):
				if np.isclose([fod[i]], [dsp], atol=1.0):
					allele.set_fodcag(fod[i])

		for item in [[primary_fod_cag, primary_dsp_cag, primary_allele],
					 [secondary_fod_cag, secondary_dsp_cag, secondary_allele]]:
			dimension_checker(item)
			primary_fod_cag = [primary_allele.get_fodcag()]; secondary_fod_cag = [secondary_allele.get_fodcag()]

		##
		## Subfunctions
		def read_comparison(val1, val2):
			if np.isclose(val1, val2, atol=1):
				return val2
			else:
				return val1

		def ensure_integrity():

			##
			## Ensure integrity
			inner_pass = True
			try:
				if not primary_dsp_ccg == int(primary_fod_ccg):
					if read_comparison(primary_dsp_ccg, int(primary_fod_ccg)) == primary_fod_ccg:
						self.sequencepair_object.get_primaryallele().set_fodccg(primary_dsp_ccg)
						inner_pass = True
					else:
						inner_pass = False
			except TypeError:
				self.sequencepair_object.get_primaryallele().set_fodccg(primary_dsp_ccg)
				inner_pass = True

			try:
				if not primary_dsp_cag == int(primary_fod_cag):
					if read_comparison(primary_dsp_cag, int(primary_fod_cag)) == primary_fod_cag:
						self.sequencepair_object.get_primaryallele().set_fodcag(primary_dsp_cag)
						inner_pass = True
					else:
						inner_pass = False
			except TypeError:
				self.sequencepair_object.get_primaryallele().set_fodcag(primary_dsp_cag)
				inner_pass = True

			try:
				if not secondary_dsp_ccg == int(secondary_fod_ccg):
					if read_comparison(secondary_dsp_ccg, int(secondary_fod_ccg)) == secondary_fod_ccg:
						self.sequencepair_object.get_secondaryallele().set_fodccg(secondary_dsp_ccg)
						inner_pass = True
					else:
						inner_pass = False
			except TypeError:
				self.sequencepair_object.get_secondaryallele().set_fodccg(secondary_dsp_ccg)
				inner_pass = True

			try:
				if not secondary_dsp_cag == int(secondary_fod_cag):
					if read_comparison(secondary_dsp_cag, int(secondary_fod_cag)) == secondary_fod_cag:
						self.sequencepair_object.get_secondaryallele().set_fodcag(secondary_dsp_cag)
						inner_pass = True
					else:
						inner_pass = False
			except TypeError:
				self.sequencepair_object.get_secondaryallele().set_fodcag(secondary_dsp_cag)

			return inner_pass

		##
		## Brute force zygosity
		if not (primary_fod_ccg == secondary_fod_ccg) and (ccg_zygstate == 'HOMO' or ccg_zygstate == 'HOMO*' or ccg_zygstate == 'HOMO+'):
			self.zygosity_state = 'HETERO'; ccg_zygstate = 'HETERO'
			self.sequencepair_object.set_svm_failure(True)
		if (primary_fod_ccg == secondary_fod_ccg) and ccg_zygstate == 'HETERO':
			self.zygosity_state = 'HOMO'; ccg_zygstate = 'HOMO'
			self.sequencepair_object.set_svm_failure(True)

		##
		## Check for potential homozygous haplotype/neighbouring peak
		if ccg_zygstate == 'HOMO' and np.isclose(primary_dsp_cag, secondary_dsp_cag, atol=1):
			primary_target = pri_distro_split['CCG{}'.format(primary_allele.get_ccg())]
			secondary_target = sec_distro_split['CCG{}'.format(secondary_allele.get_ccg())]

			primary_reads = primary_target[primary_allele.get_cag()-1]
			secondary_reads = secondary_target[secondary_allele.get_cag()-1]
			diff = abs(primary_reads-secondary_reads)
			pcnt = (diff/max([primary_reads, secondary_reads]))

			## If read count is so close (and distance is atol=1)
			## Neighbouring peak...
			if 0.0 < pcnt < 0.20:
				self.sequencepair_object.set_neighbouringpeaks(True)
				pass_vld = ensure_integrity()
				return pass_vld
			if np.isclose([pcnt], [0.25], atol=0.05):
				if 0.0 < pcnt <= 0.30:
					self.sequencepair_object.set_neighbouringpeaks(True)
					pass_vld = ensure_integrity()
					return pass_vld
				else:
					self.sequencepair_object.set_homozygoushaplotype(True)
					self.sequencepair_object.set_secondary_allele(self.sequencepair_object.get_primaryallele())
					##no need to call ensure_integrity as secondary allele is a copy of primary object
					return True

			## fucky bug with homozygotes not filtering properly
			## leaving unassigned value from FOD
			try:
				primary_fod_cag.all(); secondary_fod_cag.all()
			except AttributeError:
				secondary_fod_cag = primary_fod_cag

			if primary_fod_cag == secondary_fod_cag:

				self.sequencepair_object.set_homozygoushaplotype(True)
				self.sequencepair_object.set_secondary_allele(self.sequencepair_object.get_primaryallele())
				##no need to call ensure_integrity as secondary allele is a copy of primary object
				for allele in [self.sequencepair_object.get_primaryallele(), self.sequencepair_object.get_secondaryallele()]:
					if allele.get_peakreads() < 250:
						allele.set_fatalalignmentwarning(True)
						self.sequencepair_object.set_fatalreadallele(False)
					else:
						allele.set_fatalalignmentwarning(False)
						self.sequencepair_object.set_fatalreadallele(False)
				## Check if homozygous haplotype & that FW/RV distributions agree...
				if self.sequencepair_object.get_homozygoushaplotype():
					inferred_fwarray = []
					## infer CCG from fw dist (sum x200)
					for contig, distribution in pri_distro_split.iteritems():
						inferred_fwarray.append(sum(distribution))
					## get values for 'peaks'
					inferred_fwarray = np.asarray(inferred_fwarray)
					fwidx = int(np.where(inferred_fwarray == max(inferred_fwarray))[0][0])
					rvidx = int(np.where(primary_allele.get_rvarray() == max(primary_allele.get_rvarray()))[0][0])
					## compare against FOD
					fwpeak = peakutils.indexes(inferred_fwarray, thres=0.15, min_dist=1)
					rvpeak = peakutils.indexes(primary_allele.get_rvarray(), thres=0.15, min_dist=1)
					## suspected peak match, inspect for further noise
					if fwidx == rvidx:
						if len(rvpeak) > len(fwpeak):
							self.sequencepair_object.set_ccguncertainty(True)
						if len(fwpeak) > len(rvpeak):
							self.sequencepair_object.set_ccguncertainty(True)
				return pass_vld

		##
		## Check for diminished peaks (incase DSP failure / read count is low)
		## Primary read info
		primary_dist = split_cag_target(primary_allele.get_fwarray())
		primary_target = primary_dist['CCG{}'.format(primary_allele.get_ccg())]
		primary_reads = primary_target[primary_allele.get_cag() - 1]
		primary_total = sum(primary_target)
		## Secondary read info
		secondary_dist = split_cag_target(secondary_allele.get_fwarray())
		secondary_target = secondary_dist['CCG{}'.format(secondary_allele.get_ccg())]
		secondary_reads = secondary_target[secondary_allele.get_cag() - 1]
		secondary_total = sum(secondary_target)

		## Set specifics for zygstate
		peak_total = sum([primary_reads, secondary_reads]); dist_total = 0
		if ccg_zygstate == 'HOMO':
			dist_total = sum([primary_total])
		if ccg_zygstate == 'HOMO*' or ccg_zygstate == 'HOMO+' or ccg_zygstate == 'HETERO':
			dist_total = sum([primary_total, secondary_total])

		##
		## In the case where the peak isn't 65% of all current_ccg distribution reads
		## check that there's no diminished peak (i.e. very small expanded peak)
		## Except when atypical; due to re-alignment occurring to a specific ref
		## diminished peaks won't occur
		if not peak_total/dist_total >= 0.65:
			if np.isclose([peak_total/dist_total], [0.65], atol=0.175):
				pass
			elif primary_fod_ccg == secondary_fod_ccg and primary_dsp_cag != secondary_dsp_cag:
				primary_target = pri_distro_split['CCG{}'.format(primary_allele.get_ccg())]
				split_target = primary_target[primary_allele.get_cag()+5:-1]
				difference_buffer = len(primary_target)-len(split_target)
				fod_failstate, cag_diminished = self.peak_detection(primary_allele, split_target, 1, 'CAGDim')
				while fod_failstate:
					fod_failstate, cag_diminished = self.peak_detection(primary_allele, split_target, 1, 'CAGDim', fod_recall=True)
				if split_target[cag_diminished] > 100:
					if not primary_allele.get_allelestatus()=='Atypical' and not secondary_allele.get_allelestatus()=='Atypical':
						## bypass integrity checks
						secondary_allele.set_cagval(int(cag_diminished+difference_buffer-1))
						secondary_allele.set_fodcag(int(cag_diminished+difference_buffer-1))
						secondary_allele.set_fodoverwrite(True)
						for peak in [primary_reads, secondary_reads]:
							if peak < 750:
								self.sequencepair_object.set_diminishedpeaks(True)
						return pass_vld

		return pass_vld

	def inspect_peaks(self):

		primary_allele = self.sequencepair_object.get_primaryallele()
		secondary_allele = self.sequencepair_object.get_secondaryallele()
		for allele in [primary_allele, secondary_allele]:

			distribution_split = split_cag_target(allele.get_fwarray())
			target = distribution_split['CCG{}'.format(allele.get_ccg())]
			linspace = np.linspace(0,199,200)

			##
			## Here we set the peak reads for this allele (i.e. number of reads aligned to N)
			## In the case of WOEFUL atypical re-alignments, sometimes more than one 'peak' is found
			## Hence we detect for that, raise an exception as this allele is un-genotype-able
			try:
				if allele.get_peakreads() < 250:
					allele.set_fatalalignmentwarning(True)
					self.sequencepair_object.set_fatalreadallele(True)
			except ValueError:
				if not len(allele.get_peakreads()) == 1:
					self.sequencepair_object.set_atypical_alignmentwarning(True)
					raise Exception('Atypical re-alignment inaccuracy. Cannot genotype.')
			##
			## fucking weird interp bug filtering
			## Interp a gaussian to suspected peak
			warnings.filterwarnings('error')
			try:
				warnings.warn(Warning())
				peaks_interp = peakutils.interpolate(linspace, target, ind=[allele.get_fodcag() - 1])
				if np.isclose([peaks_interp], [allele.get_fodcag() - 1], rtol=0.5):
					interp_distance = abs(peaks_interp - float(allele.get_fodcag()) - 1)
					allele.set_interpdistance(interp_distance[0])
				else:
					allele.raise_interpolation_warning(True)
			except Warning:
				allele.raise_interpolation_warning(True)
				pass

			##
			## Calculate % of reads located near peak
			spread_reads = sum(target[allele.get_cag()-6:allele.get_cag()+5])
			spread_pcnt = (spread_reads/sum(target))
			allele.set_vicinityreads(spread_pcnt)

			##
			## Calculate peak dropoff
			nminus = target[allele.get_cag()-2]; n = target[allele.get_cag()-1]; nplus = target[allele.get_cag()]
			nminus_overn = nminus/n; nplus_overn = nplus/n
			dropoff_list = [nminus_overn, nplus_overn]
			allele.set_immediate_dropoff(dropoff_list)

			##
			## Sometimes, alignment parameters can result in invalid genotyping (i.e. 2 peaks when expecting 1)
			## Test for this, inform user..
			if self.zygosity_state == 'HETERO':
				major = max(target)
				majoridx = np.where(target == major)[0][0]
				minor = max(n for n in target if n != major)
				minoridx = np.where(target == minor)[0][0]
				thresh = (major/100)*55

				if abs(majoridx-minoridx) > 2:
					if np.isclose([major],[minor], atol=thresh):
						allele.set_unexpectedpeaks(True)
						self.pass_vld = False

			##
			## Slippage
			## Gather from N-2:N-1, sum and ratio:N
			nmt = allele.get_cag() - 3; nmo = allele.get_cag() - 1
			slip_ratio = (sum(target[nmt:nmo])) / target[allele.get_cag() - 1]
			allele.set_backwardsslippage(slip_ratio)

			rv_ratio = (target[allele.get_fodcag()-2]/target[allele.get_fodcag()-1])
			fw_ratio = (target[allele.get_fodcag()]/target[allele.get_fodcag()-1])
			if not self.sequencepair_object.get_homozygoushaplotype() and not self.sequencepair_object.get_neighbouringpeaks():
				if np.isclose([fw_ratio], [0.85], atol=0.075):
					if rv_ratio > 0.65:
						allele.set_fodcag(allele.get_fodcag()+1)
						allele.set_slippageoverwrite(True)
				if np.isclose([rv_ratio], [0.80], atol=0.150):
					if fw_ratio > 0.65:
						allele.set_fodcag(allele.get_fodcag()-1)
						allele.set_slippageoverwrite(True)
			##
			## If we're not homozygous or neighbouring, 'normal' peaks..
			## Check dropoffs are legitimate and 'clean'
			if not self.sequencepair_object.get_homozygoushaplotype() and not self.sequencepair_object.get_neighbouringpeaks():
				self.close_check(allele, nminus_overn, [0.25], 0.02, 1, state='minus') ## inform user
				self.close_check(allele, nplus_overn, [0.05], 0.02, 1, state='plus')   ## inform user
				self.close_check(allele, nminus_overn, [0.35], 0.04, 2, state='minus') ## warn user
				self.close_check(allele, nplus_overn, [0.15], 0.03, 2, state='plus')   ## warn user
				self.close_check(allele, nminus_overn, [0.45], 0.05, 3, state='minus') ## severe warning
				self.close_check(allele, nplus_overn, [0.27], 0.03, 3, state='plus')   ## severe warning
				self.close_check(allele, nminus_overn, [0.60], 0.05, 4, state='minus') ## extreme warning
				self.close_check(allele, nplus_overn, [0.37], 0.05, 4, state='plus')   ## extreme warning
				self.close_check(allele, nminus_overn, [0.75], 0.05, 5, state='minus') ## failure
				self.close_check(allele, nplus_overn, [0.65], 0.05, 5, state='plus')   ## failure
				if nminus_overn > 0.75: allele.set_nminuswarninglevel(6); self.pass_vld = False ## failure
				if nplus_overn > 0.65: allele.set_npluswarninglevel(6); self.pass_vld = False   ## failure
			else:
				allele.set_nminuswarninglevel(2)
				allele.set_npluswarninglevel(2)

			##
			## Somatic mosaicism
			## Gather from N+1:N+10, sum and ratio:N
			npo = allele.get_cag(); npt = allele.get_cag()+10
			if allele.get_header() == 'PRI': dist = self.primary_original
			if allele.get_header() == 'SEC': dist = self.secondary_original
			somatic_ratio = (sum(dist[npo:npt]))/dist[allele.get_cag()-1]
			allele.set_somaticmosaicism(somatic_ratio)

			##
			## If we get here; alleles are valid
			allele.set_ccgvalid(True)
			allele.set_cagvalid(True)
			allele.set_genotypestatus(True)

			novel_caacag = allele.get_reflabel().split('_')[1]; novel_ccgcca = allele.get_reflabel().split('_')[2]
			allele.set_allelegenotype('{}_{}_{}_{}_{}'.format(allele.get_fodcag(), novel_caacag,
															  novel_ccgcca, allele.get_fodccg(),
															  allele.get_cct()))

			##
			## Check DSP generated allele label vs FOD results
			if int(allele.get_reflabel().split('_')[3]) != int(allele.get_fodccg()):
				allele.set_referencelabel('{}_{}_{}_{}_{}'.format(allele.get_fodcag(), novel_caacag,
															  novel_ccgcca, allele.get_fodccg(),
															  allele.get_cct()))
				allele.set_fodoverwrite(True)
			if int(allele.get_reflabel().split('_')[0]) != int(allele.get_fodcag()):
				allele.set_referencelabel('{}_{}_{}_{}_{}'.format(allele.get_fodcag(), novel_caacag,
															  novel_ccgcca, allele.get_fodccg(),
															  allele.get_cct()))
				allele.set_fodoverwrite(True)

			##
			## Homozygous from SVM, genotypes match, allele DSP differ, hzygous flag False
			## aka this sample has an large expanded allele with almost no reads
			## and the pipeline missed it
			if self.zygosity_state != 'HETERO':
				if not self.sequencepair_object.get_homozygoushaplotype():
					if allele.get_fodoverwrite():
						self.sequencepair_object.set_missed_expansion(True)
						self.sequencepair_object.set_diminishedpeaks(True)

			##
			## If failed, write intermediate data to report
			if not self.pass_vld:
				inspection_logfi = os.path.join(self.sequencepair_object.get_predictpath(),
												'{}{}'.format(allele.get_header(), 'PeakInspectionLog.txt'))
				inspection_str = '{}  {}\n{}: {}\n{}: {}\n' \
								 '{}: {}\n{}: {}\n{}: {}\n' \
								 '{}: {}\n{}: {}\n{}: {}\n' \
								 '{}: {}\n{}: {}\n'.format(
								 '>> Peak Inspection Failure','Intermediate results log',
								 'Investigating CCG', allele.get_ccg(),
								 'Interpolation warning', allele.get_interpolation_warning(),
								 'Interpolation distance', allele.get_interpdistance(),
								 'Reads (%) surrounding peak', allele.get_vicinityreads(),
								 'Peak dropoff', dropoff_list,
								 'NMinus ratio', nminus_overn,
								 'NMinus warning', allele.get_nminuswarninglevel(),
								 'NPlus ratio', nplus_overn,
								 'NPlus warning', allele.get_npluswarninglevel(),
								 'Unexpected Peaks', allele.get_unexpectedpeaks())
				with open(inspection_logfi,'w') as logfi:
					logfi.write(inspection_str)
					logfi.close()

		return self.pass_vld

	def n_align_dist(self):

		"""
		Function to align alleles of the current sample to the same n-point as any previous alleles processed
		in this run. Append the padded distributions to the same CSV file for in-depth somatic mosaicism
		studies..
		:return: fuck all
		"""

		##
		## For each allele in the sample get the target CCG distribution
		## Pad it so that N (the determined genotype) is at the same index in the output file
		## output the padded distribution and close the file
		for allele in [self.sequencepair_object.get_primaryallele(), self.sequencepair_object.get_secondaryallele()]:
			if allele.get_header() == 'PRI': target = self.primary_original
			if allele.get_header() == 'SEC': target = self.secondary_original
			fix_target = ','.join(['%.5f' % num for num in target])

			anchor = 203
			anchor_port = anchor - allele.get_cag()
			anchor_starboard = anchor_port + 200
			left_buffer = '-,'*anchor_port
			right_buffer = '-,'*(403-anchor_starboard)
			padded_dist = left_buffer+fix_target+right_buffer[:-1]

			###
			### distribution fixed but csv writing incorrect list still
			sample_output = '{},{},{},{}\n'.format(self.sequencepair_object.get_label(), allele.get_header(),
													  allele.get_allelegenotype(), padded_dist)

			with open(self.padded_target, 'a') as distfi: distfi.write(sample_output)
			distfi.close()

	def render_graphs(self):

		##
		## Data for graph rendering (prevents frequent calls/messy code [[lol irony]])
		pri_fodccg = self.sequencepair_object.get_primaryallele().get_fodccg()-1
		sec_fodccg = self.sequencepair_object.get_secondaryallele().get_fodccg()-1
		pri_fodcag = self.sequencepair_object.get_primaryallele().get_fodcag()-1
		sec_fodcag = self.sequencepair_object.get_secondaryallele().get_fodcag()-1
		pri_rvarray = self.sequencepair_object.get_primaryallele().get_rvarray()
		sec_rvarray = self.sequencepair_object.get_secondaryallele().get_rvarray()
		pri_fwarray = self.sequencepair_object.get_primaryallele().get_fwarray()
		sec_fwarray = self.sequencepair_object.get_secondaryallele().get_fwarray()

		predpath = self.sequencepair_object.get_predictpath()

		def graph_subfunction(x, y, axis_labels, xticks, peak_index, predict_path, file_handle, prefix='', graph_type=None, neg_anchor=False):

			#seaborn palette
			sns.set(style='darkgrid')

			## force fonts because matplotlib spam on some people's systems?
			plt.rcParams['pdf.fonttype']=42
			plt.rcParams['ps.fonttype']=42
			mpll = log.getLogger('matplotlib'); mpll.setLevel(log.WARNING)

			x = np.linspace(x[0],x[1],x[2])
			fig, ax = plt.subplots(figsize=(10, 6)); plt.title(prefix+self.sequencepair_object.get_label())
			plt.xlabel(axis_labels[0]); plt.ylabel(axis_labels[1])
			if graph_type == 'bar':
				## ticker forced to integers (instead of floats)
				from matplotlib.ticker import FuncFormatter
				## format xtick labels correctly (CCG vs CAG distributions)
				if neg_anchor: xtickslabel = xticks[2]
				else: xtickslabel = [i-1 for i in xticks[2]]
				## plot bar data, remove x labels
				sns.barplot(x,y); plt.xticks([])
				## re-add x labels above respective bar plots
				for p, dat in zip(ax.patches, xtickslabel):
					ax.text(p.get_x() + p.get_width() / 2., p.get_height()+25, dat, ha="center", fontsize=9)

				## ticker forced to integers (instead of floats)
				plt.gca().xaxis.set_major_formatter(FuncFormatter(lambda x, _: int(x)))
			else:
				plt.xticks(np.arange(xticks[0][0], xticks[0][1], xticks[0][2]), xticks[2])
				plt.xlim(xticks[1][0], xticks[1][1])
				pplot(x,y,peak_index)
			peak_index = [i+1 for i in peak_index]
			plt.legend(['Genotype: {}'.format(peak_index)])
			plt.savefig(os.path.join(predict_path, file_handle), format='pdf')
			plt.close()

		def pagemerge_subfunction(graph_list, prediction_path, ccg_val, header=None, hplus=False):

			##
			## Readers and pages
			line_reader = PyPDF2.PdfFileReader(open(graph_list[0], 'rb')); line_page = line_reader.getPage(0)
			bar_reader = PyPDF2.PdfFileReader(open(graph_list[1], 'rb')); bar_page = bar_reader.getPage(0)

			##
			## Create new page (double width), append bar and line pages side-by-side
			translated_page = PyPDF2.pdf.PageObject.createBlankPage(None, bar_page.mediaBox.getWidth()*2, bar_page.mediaBox.getHeight())
			translated_page.mergeScaledTranslatedPage(bar_page, 1, 720, 0)
			translated_page.mergePage(line_page)

			##
			## Write to one PDF
			if hplus: suffix = 'AtypicalHomozyg'
			else: suffix = ''
			if not header: output_path = os.path.join(prediction_path, 'CCG{}CAGDetection_{}.pdf'.format(ccg_val, suffix))
			else: output_path = os.path.join(prediction_path, 'IntroCCG.pdf')
			writer = PyPDF2.PdfFileWriter()
			writer.addPage(translated_page)
			with open(output_path, 'wb') as f:
				writer.write(f)

			##
			## Return CAG plot path
			return output_path

		##########################################
		## SAMPLE CARD FOR GENOTYPE INFORMATION ##
		##########################################
		sample_pdf_path = os.path.join(predpath, '{}{}'.format(self.sequencepair_object.get_label(),'.pdf'))
		c = canvas.Canvas(sample_pdf_path, pagesize=(720,432))
		header_string = '{}{}'.format('Sample header: ', self.sequencepair_object.get_label())
		primary_string = '{}(CAG{}, CCG{}) ({}; {}; Confidence {}%)'.format('Primary: ', self.sequencepair_object.get_primaryallele().get_fodcag(),
												 self.sequencepair_object.get_primaryallele().get_fodccg(),
												 self.sequencepair_object.get_primaryallele().get_allelestatus(),
												 self.sequencepair_object.get_primaryallele().get_allelegenotype(),
												 int(self.sequencepair_object.get_primaryallele().get_alleleconfidence()))
		secondary_string = '{}(CAG{}, CCG{}) ({}; {}; Confidence {}%)'.format('Secondary: ', self.sequencepair_object.get_secondaryallele().get_fodcag(),
												   self.sequencepair_object.get_secondaryallele().get_fodccg(),
												   self.sequencepair_object.get_secondaryallele().get_allelestatus(),
												   self.sequencepair_object.get_secondaryallele().get_allelegenotype(),
												   int(self.sequencepair_object.get_secondaryallele().get_alleleconfidence()))

		##########################################################
		## Create canvas for sample 'intro card'				##
		## Set font colour depending on subsample/invalid/valid ##
		## invalid == atypical allele, no realignment			##
		## valid == atypical allele, realigned					##
		##########################################################
		if self.sequencepair_object.get_subsampleflag():
			if self.sequencepair_object.get_subsampleflag() == '0.05**':
				pass
			elif self.sequencepair_object.get_automatic_DSPsubsample():
				pass
			elif float(self.sequencepair_object.get_subsampleflag()) >= 0.5:
				pass
			else:
				c.setFillColorRGB(75, 0, 130)
				c.drawCentredString(360, 186, '!! Genotype derived from significantly subsampled data !!')
		if self.invalid_data:
			c.setFillColorRGB(255, 0, 0)
			c.drawCentredString(360, 196, '!! Atypical alleles without re-alignment !!')
		if not self.invalid_data:
			c.setFillColorRGB(0, 0, 0)
		c.drawCentredString(360, 256, header_string)
		c.drawCentredString(360, 236, primary_string)
		c.drawCentredString(360, 216, secondary_string)
		c.save()

		###############################################
		## CCG heterozygous example					 ##
		## i.e. CCG two peaks, one CAG dist per peak ##
		###############################################
		if self.zygosity_state == 'HETERO' or self.zygosity_state == 'HOMO*' or self.zygosity_state == 'HOMO+':

			##
			## Render CCG graph, append path to allele path list
			## Merge intro_ccg card with sample CCG graph
			## Append merged intro_ccg to heterozygous list
			hetero_graphs = []; ccg_peaks = [int(pri_fodccg),int(sec_fodccg)]
			concat = np.asarray([a + b for a, b in zip(pri_rvarray,sec_rvarray)])
			graph_subfunction([0, 21, 20], concat, ['CCG Value', 'Read Count'], ([1, 20, 1], [1, 20], range(1,21)),
							  ccg_peaks, predpath, 'CCGDetection.pdf', graph_type='bar', neg_anchor=True)
			intro_card = pagemerge_subfunction([sample_pdf_path, os.path.join(predpath, 'CCGDetection.pdf')],
													predpath, ccg_val=0, header=True)
			hetero_graphs.append(intro_card)
			plt.close()

			##
			## For each CCG allele in this heterozygous sample
			for allele in [self.sequencepair_object.get_primaryallele(),self.sequencepair_object.get_secondaryallele()]:

				##
				## Data for this allele (peak detection graph)
				temp_graphs = []

				target_distro = []
				if allele.get_header() == 'PRI':
					target_distro = self.primary_original
				if allele.get_header() == 'SEC':
					target_distro = self.secondary_original

				if self.zygosity_state == 'HOMO+':
					for i in range(0, len(target_distro)):
						if i != allele.get_cag() - 1:
							removal = (target_distro[i] / 100) * 75
							target_distro[i] -= removal

				if allele.get_rewrittenccg() is not None:
					peak_filename = 'CCG{}-CAGDetection_atypical_ccgdiff.pdf'.format(allele.get_fodccg())
					peak_prefix = '(CCG{}**) '.format(allele.get_fodccg())
				elif allele.get_unrewrittenccg() is not None:
					peak_filename = 'CCG{}-CAGDetection_atypical_ccgsame.pdf'.format(allele.get_fodccg())
					peak_prefix = '(CCG{}++) '.format(allele.get_fodccg())
				else:
					peak_filename = 'CCG{}-CAGDetection.pdf'.format(allele.get_fodccg())
					peak_prefix = '(CCG{}) '.format(allele.get_fodccg())
				peak_graph_path = os.path.join(predpath, peak_filename)
				## Render the graph, append to list, close plot
				graph_subfunction([0, 199, 200], target_distro, ['CAG Value', 'Read Count'],
								  ([1, 200, 50], [1, 200], [0,50,100,150,200]), [np.int64(allele.get_fodcag() - 1)],
								  predpath, peak_filename, prefix=peak_prefix)
				temp_graphs.append(peak_graph_path); plt.close()

				##
				## Inspect the peak (subslice)
				slice_range = range(allele.get_fodcag()-4, allele.get_fodcag()+7)
				if allele.get_rewrittenccg():
					slice_filename = 'CCG{}-Peak_atypical_ccgdiff.pdf'.format(allele.get_fodccg())
					slice_prefix = '(CCG{}**) '.format(allele.get_ccg())
				elif allele.get_unrewrittenccg():
					slice_filename = 'CCG{}-Peak_atypical_ccgsame.pdf'.format(allele.get_fodccg())
					slice_prefix = '(CCG{}++) '.format(allele.get_ccg())
				else:
					slice_filename = 'CCG{}-Peak.pdf'.format(allele.get_fodccg())
					slice_prefix = '(CCG{}) ' .format(allele.get_ccg())
				sub = target_distro[np.int64(allele.get_fodcag()-6):np.int64(allele.get_fodcag()+5)]
				## Render the graph, append to list, close plot
				graph_subfunction([0,10,11], sub, ['CAG Value', 'Read Count'], ([1,11,1], [1,11], slice_range),
								  [np.int64(allele.get_fodcag()-1)], predpath,slice_filename, prefix=slice_prefix, graph_type='bar')
				temp_graphs.append(os.path.join(predpath,slice_filename)); plt.close()

				##
				## Merge 'allele sample' into one page
				ccg_val = allele.get_fodccg()
				if (allele.get_unrewrittenccg() is not None) or (allele.get_rewrittenccg() is not None): hplus = True
				else: hplus = False
				merged_graph = pagemerge_subfunction(temp_graphs, predpath, ccg_val, hplus=hplus)
				hetero_graphs.append(merged_graph)

			self.sequencepair_object.get_primaryallele().set_allelegraphs(hetero_graphs)
			self.sequencepair_object.get_secondaryallele().set_allelegraphs(hetero_graphs)

		##############################################
		## CCG homozygous example					##
		## i.e. CCG one peak, one CAG dist per peak ##
		##############################################
		if self.zygosity_state == 'HOMO':

			##
			##Data for homozygous graph(s)
			homo_graphs = []
			page_graphs = []
			target_ccg = 'CCG{}'.format(self.sequencepair_object.get_primaryallele().get_ccg())
			## Peak data
			peak_filename = 'CCG{}-CAGDetection.pdf'.format(self.sequencepair_object.get_primaryallele().get_fodccg())
			peak_prefix = '(CCG{}) '.format(self.sequencepair_object.get_primaryallele().get_ccg())
			altpeak_filename = 'CCG{}-Peak.pdf'.format(self.sequencepair_object.get_primaryallele().get_fodccg())
			ccg_peaks = [int(pri_fodccg),int(sec_fodccg)]; cag_peaks = [int(pri_fodcag),int(sec_fodcag)]
			distribution_split = split_cag_target(pri_fwarray); target_distro = self.primary_original
			## Subslice data
			pri_cag = self.sequencepair_object.get_primaryallele().get_cag()
			sec_cag = self.sequencepair_object.get_secondaryallele().get_cag()
			upper = max([pri_cag, sec_cag])
			if self.sequencepair_object.get_homozygoushaplotype(): lower = upper
			else: lower = max(n for n in [pri_cag, sec_cag] if n != upper)
			sub = target_distro[lower-6:upper+5]
			slice_range = range(lower-4,upper+7)

			##
			## Render the graph, append to list, close plot
			## Merge intro_ccg card with sample CCG graph
			## Append merged intro_ccg to homozygous list, append line/bar peak to page list
			graph_subfunction([0, 21, 20], pri_rvarray, ['CCG Value', 'Read Count'], ([1, 20, 1], [1, 20], range(1,21)),
							  ccg_peaks, predpath, 'CCGDetection.pdf', graph_type='bar', neg_anchor=True); plt.close()
			graph_subfunction([0, 199, 200], target_distro, ['CAG Value', 'Read Count'],
							  ([1, 200, 50], [1, 200], [0,50,100,150,200]), cag_peaks, predpath,
							  peak_filename, prefix=peak_prefix); plt.close()
			graph_subfunction([0, len(sub)-1, len(sub)], sub, ['CAG Value', 'Read Count'],
							  ([1, len(sub), 1], [1, len(sub)], slice_range), cag_peaks, predpath, altpeak_filename,
							  prefix=peak_prefix, graph_type='bar'); plt.close()
			intro_card = pagemerge_subfunction([sample_pdf_path, os.path.join(predpath, 'CCGDetection.pdf')],
													predpath, ccg_val=0, header=True)
			homo_graphs.append(intro_card)
			page_graphs.append(os.path.join(predpath, peak_filename))
			page_graphs.append(os.path.join(predpath, altpeak_filename))

			##
			## Merge 'allele sample' into one page
			ccg_val = self.sequencepair_object.get_primaryallele().get_fodccg()
			merged_graph = pagemerge_subfunction(page_graphs, predpath, ccg_val)

			## Combine CCG and CAG graphs
			homo_graphs.append(merged_graph)
			self.sequencepair_object.get_primaryallele().set_allelegraphs(homo_graphs)
			self.sequencepair_object.get_secondaryallele().set_allelegraphs(homo_graphs)

		############################################
		## Merge graphs into a single summary PDF ##
		############################################

		##
		## Allele graphs
		## Ensure uniqueness of entries in primary/secondary (i.e. no duplicating CCG graph)
		primary_graphs = self.sequencepair_object.get_primaryallele().get_allelegraphs()[0]
		secondary_graphs = self.sequencepair_object.get_primaryallele().get_allelegraphs()[0]
		sample_graphs = primary_graphs + secondary_graphs
		def set_orderpreserve(seq, idfun=None):
			if idfun is None:
				def idfun(x): return x
			seen = {}
			result = []
			for item in seq:
				marker = idfun(item)
				if marker in seen: continue
				seen[marker] = 1
				result.append(item)
			return result
		uniques = set_orderpreserve(sample_graphs)

		##
		## Merge alleles together
		merger = PyPDF2.PdfFileMerger()
		for pdf in uniques:
			merger.append(pdf)
		merger.write(sample_pdf_path)
		merger.close()

		##
		## Remove individual plots
		clean_target = []
		for target_file in os.listdir(predpath):
			if target_file.endswith(".pdf"):
				clean_target.append(os.path.join(predpath, target_file))
		for rmpdf in clean_target:
			if not '{}{}'.format(self.sequencepair_object.get_label(),'.pdf') in rmpdf:
				os.remove(rmpdf)

	def calculate_score(self):

		##
		## For both alleles
		for allele in [self.sequencepair_object.get_primaryallele(), self.sequencepair_object.get_secondaryallele()]:
			allele_log_fi = os.path.join(self.sequencepair_object.get_predictpath(), '{}{}'.format(allele.get_header(), '_PenaltiesLog.txt'))
			with open(allele_log_fi, 'a') as penfi:
				penfi.write('{}, {}\n'.format('Flag/Warning','Score Penalty'))

				##
				## Start score high, deduct for questionable calls..
				allele_confidence = 100

				##
				## Sample based genotyping flags
				if self.sequencepair_object.get_recallcount() == 7: allele_confidence -= 25; penfi.write('{}, {}\n'.format('Recall Count','-25'))
				if 7 > self.sequencepair_object.get_recallcount() > 4: allele_confidence -= 15; penfi.write('{}, {}\n'.format('Recall Count','-15'))
				if 4 > self.sequencepair_object.get_recallcount() > 0: allele_confidence -= 5; penfi.write('{}, {}\n'.format('Recall Count', '-5'))
				else: allele_confidence += 0; penfi.write('{}, {}\n'.format('Recall Count', '+0'))

				if self.sequencepair_object.get_homozygoushaplotype():
					allele_confidence -= 25; penfi.write('{}, {}\n'.format('Homozygous Haplotype','-25'))
				elif self.sequencepair_object.get_neighbouringpeaks():
					allele_confidence -= 15; penfi.write('{}, {}\n'.format('Neighbouring Peaks', '-15'))
				else: allele_confidence += 0; penfi.write('{}, {}\n'.format('Normal Data','+0'))

				if self.sequencepair_object.get_diminishedpeaks():
					allele_confidence -= 10; penfi.write('{}, {}\n'.format('Diminished Peaks','-10'))
				if allele.get_fodoverwrite():
					allele_confidence -= 20; penfi.write('{}, {}\n'.format('Differential Overwrite','-20'))
					if self.sequencepair_object.get_missed_expansion():
						allele_confidence -= 50; penfi.write('{}, {}\n'.format('Missed Expansion', '-50'))

				##
				## Allele based genotyping flags
				## Allele typical/atypical structure
				if allele.get_allelestatus() == 'Atypical':
					allele_confidence -= 5; penfi.write('{}, {}\n'.format('Atypical Allele','-5'))
					if np.isclose([float(allele.get_atypicalpcnt())],[50.00],atol=5.00):
						allele_confidence -= 30; penfi.write('{}, {}\n'.format('Atypical reads (50%)','-30'))
					if np.isclose([float(allele.get_atypicalpcnt())],[80.00],atol=20.00):
						allele_confidence += 15; penfi.write('{}, {}\n'.format('Atypical reads (80%>)','+15'))
				if allele.get_allelestatus() == 'Typical':
					allele_confidence += 1; penfi.write('{}, {}\n'.format('Typical Allele', '+1'))
					if np.isclose([float(allele.get_typicalpcnt())],[50.00],atol=5.00):
						allele_confidence -= 30; penfi.write('{}, {}\n'.format('Typical reads (50%)','-30'))
					if np.isclose([float(allele.get_typicalpcnt())],[80.00],atol=15.00):
						allele_confidence += 2; penfi.write('{}, {}\n'.format('Typical reads (80%>)','+2'))

				##
				## Total reads in sample..
				if allele.get_totalreads() > 10000:	allele_confidence += 5; penfi.write('{}, {}\n'.format('High total read count', '+5'))
				elif allele.get_totalreads() < 1000: allele_confidence -= 15; penfi.write('{}, {}\n'.format('Low total read count', '-15'))
				else: allele_confidence += 1; penfi.write('{}, {}\n'.format('Normal total read count','+1'))

				##
				## Variance of distribution utilised
				if allele.get_vicinityreads()*100 > 85.00: allele_confidence += 1; penfi.write('{}, {}\n'.format('Reads near peak','+1'))
				elif 84.99 > allele.get_vicinityreads()*100 > 65.00: allele_confidence -= 10; penfi.write('{}, {}\n'.format('Reads near peak','-10'))
				elif 64.99 > allele.get_vicinityreads()*100 > 45.00: allele_confidence -= 15; penfi.write('{}, {}\n'.format('Reads near peak','-15'))
				elif 44.99 > allele.get_vicinityreads()*100 > 00.00: allele_confidence -= 20; penfi.write('{}, {}\n'.format('Reads near peak','-20'))

				##
				## Backwards slippage ratio ([N-2:N-1]/N]
				if 0.00 < allele.get_backwardsslippage() < 0.10: allele_confidence += 5; penfi.write('{}, {}\n'.format('Backwards slippage','+5'))
				elif 0.10 < allele.get_backwardsslippage() < 0.25: allele_confidence += 1; penfi.write('{}, {}\n'.format('Backwards slippage','+1'))
				elif allele.get_backwardsslippage() > 0.25: allele_confidence -= 1; penfi.write('{}, {}\n'.format('Backwards slippage','-1'))
				elif allele.get_backwardsslippage() > 0.45: allele_confidence -= 5; penfi.write('{}, {}\n'.format('Backwards slippage','-5'))
				elif allele.get_backwardsslippage() > 0.65: allele_confidence -= 10; penfi.write('{}, {}\n'.format('Backwards slippage','-10'))
				elif 0.65 < allele.get_backwardsslippage() < 1.00: allele_confidence -= 25; penfi.write('{}, {}\n'.format('Backwards slippage','-25'))
				if allele.get_slippageoverwrite(): allele_confidence -= 25; penfi.write('{}, {}\n'.format('Slippage overwrite','-25'))

				##
				## Somatic mosiacisim ratio ([N+1:N+10]/N]
				if 0.000 < allele.get_somaticmosaicism() < 0.010: allele_confidence += 10; penfi.write('{}, {}\n'.format('Somatic mosaicism','+10'))
				elif 0.010 < allele.get_somaticmosaicism() < 0.015: allele_confidence += 5; penfi.write('{}, {}\n'.format('Somatic mosaicism','+5'))
				elif allele.get_somaticmosaicism() > 0.015: allele_confidence -= 1; penfi.write('{}, {}\n'.format('Somatic mosaicism','-1'))
				elif allele.get_somaticmosaicism() > 0.025: allele_confidence -= 10; penfi.write('{}, {}\n'.format('Somatic mosaicism','-10'))
				elif allele.get_somaticmosaicism() > 0.035: allele_confidence -= 15; penfi.write('{}, {}\n'.format('Somatic mosaicism','-15'))
				elif 0.035 < allele.get_somaticmosaicism() < 0.100: allele_confidence -= 20; penfi.write('{}, {}\n'.format('Somatic mosaicism','-20'))
				elif allele.get_somaticmosaicism() > 0.100: allele_confidence -= 30; penfi.write('{}, {}\n'.format('Somatic mosaicism','-30'))

				##
				## Peak calling thresholds
				for contig in [allele.get_ccgthreshold(), allele.get_cagthreshold()]:
					if contig != 0.5:
						if 0.5 > contig > 0.3: allele_confidence -= 5; penfi.write('{}, {}\n'.format('Peak calling threshold','-5'))
						if 0.3 > contig > 0.0: allele_confidence -= 10; penfi.write('{}, {}\n'.format('Peak calling threshold','-10'))
					else: allele_confidence += 1; penfi.write('{}, {}\n'.format('Peak calling threshold','+1'))

				##
				## Peak dropoff warnings
				for peak_position_error in [allele.get_nminuswarninglevel(), allele.get_npluswarninglevel()]:
					if peak_position_error == 0: allele_confidence += 5; penfi.write('{}, {}\n'.format('Surrounding read ratio','+5'))
					elif peak_position_error == 1: allele_confidence -= 5; penfi.write('{}, {}\n'.format('Surrounding read ratio', '-5'))
					elif 2 >= peak_position_error > 1: allele_confidence -= 7; penfi.write('{}, {}\n'.format('Surrounding read ratio','-7'))
					elif peak_position_error >= 5: allele_confidence -= 10; penfi.write('{}, {}\n'.format('Surrounding read ratio','-10'))
					else: allele_confidence -= 15; penfi.write('{}, {}\n'.format('Surrounding read ratio','-15'))

				##
				## Multiply score by a factor if reads were subsampled
				if self.sequencepair_object.get_subsampleflag() and not self.sequencepair_object.get_subsampleflag() == '0.05**':
					subsample_penalty = []; utilised_subsample_penalty = 0.0; context_penalty = 0.0
					if 0 <= self.sequencepair_object.get_totalseqreads() <= 2000: subsample_penalty = [0.35,0.40,0.45]
					if 2000 <= self.sequencepair_object.get_totalseqreads() <= 5000: subsample_penalty = [0.55,0.65,0.70]
					if 5000 <= self.sequencepair_object.get_totalseqreads() <= 10000: subsample_penalty = [0.75,0.85,0.95]
					if self.sequencepair_object.get_totalseqreads() > 10000: subsample_penalty = [1.0,1.0,1.0]

					if 0.1 <= self.sequencepair_object.get_subsampleflag() <= 0.3:
						utilised_subsample_penalty = subsample_penalty[0]
					elif 0.3 <= self.sequencepair_object.get_subsampleflag() <= 0.6:
						utilised_subsample_penalty = subsample_penalty[1]
					elif 0.6 <= self.sequencepair_object.get_subsampleflag() <= 0.9:
						utilised_subsample_penalty = subsample_penalty[2]
					else:
						utilised_subsample_penalty = 1

					allele_read_ratio = allele.get_totalreads() / self.sequencepair_object.get_totalseqreads()
					if np.isclose([allele_read_ratio],[0.05],atol=0.05): context_penalty = 15
					if np.isclose([allele_read_ratio],[0.15],atol=0.05): context_penalty = 10
					if np.isclose([allele_read_ratio],[0.25],atol=0.05): context_penalty = 7
					if np.isclose([allele_read_ratio],[0.35],atol=0.05): context_penalty = 5
					if np.isclose([allele_read_ratio],[0.45],atol=0.05): context_penalty = 5
					if np.isclose([allele_read_ratio],[0.55],atol=0.05): context_penalty = 1

					##
					## Due to the nature of neighbour/homozygous peaks, don't cound surrounding reads...
					if self.sequencepair_object.get_neighbouringpeaks() or self.sequencepair_object.get_homozygoushaplotype():
						pass
					else:
						allele_confidence = allele_confidence * utilised_subsample_penalty
						allele_confidence -= context_penalty

					penfi.write('{}, *{}\n'.format('Subsample demultiplier', utilised_subsample_penalty))
					penfi.write('{}, -{}\n'.format('Read Ratio Context', context_penalty))

				##
				## Mapping percentage
				for map_pcnt in [allele.get_fwalnpcnt(), allele.get_rvalnpcnt()]:
					if map_pcnt > 90: allele_confidence += 2; penfi.write('{}, {}\n'.format('Mapping percentage', '+2'))
					elif 85 < map_pcnt < 90: allele_confidence += 1; penfi.write('{}, {}\n'.format('Mapping percentage', '+1'))
					else: allele_confidence -= 2; penfi.write('{}, {}\n'.format('Mapping percentage', '-2'))

				##
				## Warning penalty.. if triggered, no confidence
				if self.warning_triggered: allele_confidence -= 10; penfi.write('{}, {}\n'.format('Peak Inspection warning triggered','-10'))
				if allele.get_ccguncertainty(): allele_confidence -= 15; penfi.write('{}, {}\n'.format('CCG Uncertainty','-15'))
				if self.sequencepair_object.get_alignmentwarning(): allele_confidence -= 15; penfi.write('{}, {}\n'.format('Low read count alignment warning','-15'))
				if self.sequencepair_object.get_atypical_alignmentwarning(): allele_confidence -= 50; penfi.write('{}, {}\n'.format('Atypical re-alignment inaccurate','-50'))
				if allele.get_fatalalignmentwarning(): allele_confidence -= 25; penfi.write('{}, {}\n'.format('Fatal low read count alignment warning','-25'))

				##
				## Differential Confusion
				## i.e two peaks nearby, large difference between suspected but unsure whether homo or neighbour
				if allele.get_differential_confusion():
					allele_confidence -= 30; penfi.write('{}, {}\n'.format('Differential Confusion', '-30'))

				## Heuristic filtering of DSP results state check
				if not self.sequencepair_object.get_heuristicfilter():
					allele_confidence -= 20; penfi.write('{}, {}\n'.format('Heuristic filtering of alleles did not assign a secondary allele','-20'))

				##
				## If reflabel CAG and FOD CAG differ.. no confidence
				label_split = allele.get_reflabel().split('_')[0]
				if allele.get_allelestatus() == 'Atypical':
					if not np.isclose([int(allele.get_fodcag())],[int(label_split)],atol=1):
						allele_confidence = 0; penfi.write('{}, {}\n'.format('Atypical DSP:FOD inconsistency','-100'))

				##
				## Determine score (max out at 100), return genotype
				capped_confidence = sorted([0, allele_confidence, 100])[1]
				allele.set_alleleconfidence(capped_confidence)
				penfi.write('{}, {}\n\n'.format('Final score', capped_confidence))
				penfi.close()

	def contextualise(self):

		##
		## For both alleles
		for allele in [self.sequencepair_object.get_primaryallele(), self.sequencepair_object.get_secondaryallele()]:

			## get list of CAG sizes and number of reads aligned to each size1
			cag_sizes = [i for i in range(1,201)]
			aligned_distribution = allele.get_fwarray()
			raw_repeat_distribution = allele.get_fwarray().copy()
			split_distribution = split_cag_target(raw_repeat_distribution)
			allele_distribution = split_distribution['CCG{}'.format(allele.get_ccg())]
			index = allele.get_cag()-1

			## test clean up of distribution for specific alleles
			## todo implement a mechanism by which to remove individual alleles data
			## rather than just anything other than the index of the allele
			allele_clearance_range = range(allele.get_cag()-11, allele.get_cag()+11)
			for i in range(0, len(allele_distribution)):
				if i not in allele_clearance_range:
					removal = (allele_distribution[i] / 100) * 85
					allele_distribution[i] -= removal

			## make 'expanded distribution' of literal count
			expanded_distribution = []
			for index, aln_count in zip(cag_sizes, allele_distribution):
				to_add = [index] * aln_count
				expanded_distribution += to_add

			## calculate 95% confidence interval given this distribution
			a = 1.0 * np.array(expanded_distribution)
			n = len(a)
			m, se = np.mean(a), sp.stats.sem(a)
			h = se * sp.stats.t.ppf((1 + 0.95) / 2., n-1)
			lower_ci = int(round(m-h)); upper_ci = int(round(m+h))
			if lower_ci == upper_ci:
				if abs(upper_ci - int(allele.get_cag())) > 1:
					allele.set_alleleconfinterval('{}-{}'.format(lower_ci, allele.get_cag()))
				elif abs(upper_ci - int(allele.get_cag())) == 1:
					allele.set_alleleconfinterval('{}-{}'.format(allele.get_cag(), allele.get_cag()))
				else:
					allele.set_alleleconfinterval('{}-{}'.format(lower_ci, upper_ci))
			else:
				if allele.get_cag() > upper_ci:
					allele.set_alleleconfinterval('{}-{}'.format(lower_ci, allele.get_cag()))
				else:
					allele.set_alleleconfinterval('{}-{}'.format(lower_ci, upper_ci))

			## haha bad lazy programming wow
			if allele.get_alleleconfidence() < 50:
				garbage_code = allele.get_alleleconfinterval().split('-')
				lower = garbage_code[0]; upper = garbage_code[1]
				lower = int(lower)-(int(np.random.randint(1,3,size=1)[0]))
				upper = int(upper)+(int(np.random.randint(1,2,size=1)[0]))
				allele.set_alleleconfinterval('{}-{}'.format(lower,upper))


	def set_report(self):

		for allele in [self.sequencepair_object.get_primaryallele(), self.sequencepair_object.get_secondaryallele()]:

			##
			## Report path for this allele
			allele_filestring = '{}{}{}'.format(allele.get_header(),allele.get_allelestatus(), '_AlleleReport.txt')
			report_path = os.path.join(self.sequencepair_object.get_predictpath(), allele_filestring)
			allele.set_allelereport(report_path)
			report_string = '{}{}\n\n{}\n{}{}\n{}{}\n{}{}\n{}{}\n{}{}\n{}{}\n{}{}' \
							'\n\n{}\n{}{}\n{}{}\n{}{}\n{}{}\n{}{}\n{}{}\n{}{}\n{}{}\n{}{}\n{}{}\n{}{}\n{}{}\n\n' \
							'{}\n{}{}\n{}{}\n{}{}\n{}{}\n{}{}\n{}{}'.format(
							'Allele Report>> ', self.sequencepair_object.get_label(),
							'Summary Information>>',
							'Genotype: ', allele.get_allelegenotype(),
							'Confidence: ', allele.get_alleleconfidence(),
							'CCG Uncertain: ', self.sequencepair_object.get_ccguncertainty(),
							'Structure Status: ', allele.get_allelestatus(),
							'Typical Pcnt: ', allele.get_typicalpcnt(),
							'Atypical Pcnt: ', allele.get_atypicalpcnt(),
							'Total Reads: ', allele.get_totalreads(),
							'Flags>>',
							'Recall Count: ', self.sequencepair_object.get_recallcount(),
							'Homozygous Haplotype: ', self.sequencepair_object.get_homozygoushaplotype(),
							'Neighbouring Peaks: ', self.sequencepair_object.get_neighbouringpeaks(),
							'Diminished Peaks: ', self.sequencepair_object.get_diminishedpeaks(),
							'Backwards Slippage: ', allele.get_backwardsslippage(),
							'Somatic Mosaicism: ', allele.get_somaticmosaicism(),
							'Slippage Overwritten: ', allele.get_slippageoverwrite(),
							'Peak Interpolation Warning: ', allele.get_interpolation_warning(),
							'Peak Interpolation Distance: ', allele.get_interpdistance(),
							'Peak DSP Overwritten: ', allele.get_fodoverwrite(),
							'Low read-count alignment: ', self.sequencepair_object.get_alignmentwarning(),
							'Fatal low read-count: ', allele.get_fatalalignmentwarning(),
							'Data Quality>>',
							'Reads (%) surrounding peak: ', allele.get_vicinityreads(),
							'Immediate Dropoffs: ', allele.get_immediate_dropoff(),
							'N-1 Warning Level: ', allele.get_nminuswarninglevel(),
							'N+1 Warning Level: ', allele.get_npluswarninglevel(),
							'CCG Threshold: ', allele.get_ccgthreshold(),
							'CAG Threshold: ', allele.get_cagthreshold()
							)
			##
			## Write to file
			with open(report_path, 'w') as outfi:
				outfi.write(report_string)
				outfi.close()

	def get_report(self):

		self.allele_report = [self.sequencepair_object.get_primaryallele().get_allelereport(),
							  self.sequencepair_object.get_secondaryallele().get_allelereport()]
		return self.allele_report
