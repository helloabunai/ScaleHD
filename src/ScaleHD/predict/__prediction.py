from __future__ import division

#/usr/bin/python
__version__ = 0.01
__author__ = 'alastair.maxwell@glasgow.ac.uk'

##
## Generic imports
import os
import csv
import PyPDF2
import warnings
import peakutils
import matplotlib
import numpy as np
matplotlib.use('Agg')
from sklearn import svm
import matplotlib.pyplot as plt
from sklearn import preprocessing
from reportlab.pdfgen import canvas
from peakutils.plot import plot as pplot
from sklearn.multiclass import OutputCodeClassifier

##
## Backend Junk
from ..__backend import DataLoader

class AlleleGenotyping:
	def __init__(self, sequencepair_object, instance_params, training_data, atypical_logic=None):

		##
		## Allele objects and instance data
		self.sequencepair_object = sequencepair_object
		self.instance_params = instance_params
		self.training_data = training_data
		self.invalid_data = atypical_logic
		self.allele_report = ''

		##
		## Constructs that will be updated with each allele process
		self.classifier, self.encoder = self.build_zygosity_model()
		self.allele_flags = {}; self.forward_distribution = None; self.reverse_distribution = None
		self.forward_aggregate = None; self.reverse_aggregate = None
		self.expected_zygstate = None; self.zygosity_state = None

		##
		## Genotype!
		if not self.allele_validation(): raise Exception('Allele(s) failed validation. Cannot genotype..')
		if not self.determine_ccg(): raise Exception('CCG Genotyping failure. Cannot genotype..')
		if not self.determine_cag(): raise Exception('CAG Genotyping failure. Cannot genotype..')
		if not self.genotype_validation(): raise Exception('Genotype failed validation. Check output log..')
		if not self.inspect_peaks(): raise Exception('Peak Inspection failure. Check output log..')
		self.render_graphs()
		self.calculate_score()
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
								   intercept_scaling=1, verbose=0, random_state=0, max_iter=-1)

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
		## Predict the zygstate of these reshapen, noramlised 20D CCG arrays using SVM object earlier
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
		## Object for CCG split
		ccg_arrays = None

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
		distribution_dict = {}
		for i in range(0, len(cag_split)):
			distribution_dict['CCG' + str(i + 1)] = cag_split[i]

		# current_target_distribution = distribution_dict['CCG' + str(ccg_target)]
		return distribution_dict

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
			if recall_count > 7: raise Exception('7+ recalls. Unable to determine genotype.\n')
			threshold = 0.0
			if triplet_stage == 'CCG': threshold = allele_object.get_ccgthreshold()
			if triplet_stage == 'CAG': threshold = allele_object.get_cagthreshold()
			threshold -= 0.06
			utilised_threshold = max(threshold, 0.05)

		if triplet_stage == 'CCG':
			allele_object.set_ccgthreshold(utilised_threshold)
			error_boundary = 1
		if triplet_stage == 'CAGHet':
			allele_object.set_cagthreshold(utilised_threshold)
			error_boundary = 1
		if triplet_stage == 'CAGHom':
			allele_object.set_cagthreshold(utilised_threshold)
			error_boundary = 2
		##
		## Look for peaks in our distribution
		peak_indexes = peakutils.indexes(distro, thres=utilised_threshold, min_dist=peak_dist)
		fixed_indexes = np.array(peak_indexes + 1)
		if not len(fixed_indexes) == error_boundary:
			if triplet_stage == 'CAGHom' and est_dist == 1:
				self.sequencepair_object.set_homozygoushaplotype(True)
				fixed_indexes = np.asarray([fixed_indexes[0], fixed_indexes[0]])
			elif allele_object.get_cag() in fixed_indexes:
				fixed_indexes = np.asarray([x for x in fixed_indexes if x == allele_object.get_cag()])
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
			## Unlabelled distributions
			self.forward_distribution = self.scrape_distro(allele_object.get_fwdist())
			self.reverse_distribution = self.scrape_distro(allele_object.get_rvdist())
			allele_object.set_fwarray(self.forward_distribution)
			allele_object.set_rvarray(self.reverse_distribution)

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
					self.reverse_aggregate = self.pad_distribution(self.reverse_distribution, allele_object)
					allele_object.set_rvarray(self.reverse_aggregate)
				## Data has not been realigned -- brute force genotyping
				if self.invalid_data:
					self.forward_aggregate = self.distribution_collapse(self.forward_distribution)
					self.reverse_aggregate = self.distribution_collapse(self.reverse_distribution)
			self.zygosity_state = self.predict_zygstate()

			##
			## Clean up distribution for erroneous peaks
			for i in range(0, len(self.reverse_aggregate)):
				if i > allele_object.get_ccg()-1:
					removal = (i/100) * 90
					self.reverse_aggregate[i] = i-removal
				if i < allele_object.get_ccg()-1:
					removal = (i/100) * 90
					self.reverse_aggregate[i] = i-removal
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
		## Check both alleles passed validation
		if (self.sequencepair_object.get_primaryallele().get_validation()) and (
				self.sequencepair_object.get_secondaryallele().get_validation()):
			return True
		else:
			return False

	def determine_ccg(self):

		##
		## Constructs
		ccg_matches = 0; ccg_values = []; local_zygstate = None; pass_gtp = True

		##
		## For the two allele objects in this sample_pair
		## First, ensure CCG matches between DSP estimate and FOD derision
		for allele in [self.sequencepair_object.get_primaryallele(), self.sequencepair_object.get_secondaryallele()]:

			allele.set_ccgthreshold(0.50)
			fod_failstate, ccg_indexes = self.peak_detection(allele, allele.get_rvarray(), 1, 'CCG')
			while fod_failstate:
				fod_failstate, ccg_indexes = self.peak_detection(allele, allele.get_rvarray(), 1, 'CCG', fod_recall=True)

			if ccg_indexes[0] == allele.get_ccg():
				ccg_matches += 1
				allele.set_ccgvalid(True)
			ccg_values.append(ccg_indexes[0])
			allele.set_fodccg(np.asarray(ccg_indexes[0]))

		if ccg_values[0] == ccg_values[1]:
			local_zygstate = 'HOMO'
		if not ccg_values[0] == ccg_values[1]:
			local_zygstate = 'HETERO'

		self.zygosity_state = local_zygstate
		self.sequencepair_object.set_ccgzygstate(local_zygstate)
		if not local_zygstate == self.expected_zygstate:
			pass_gtp = False

		return pass_gtp

	def determine_cag(self):

		##
		## Constructs
		pass_gtp = True

		##########################
		## Heterozygous for CCG ##
		##########################
		if self.zygosity_state == 'HETERO':
			for allele in [self.sequencepair_object.get_primaryallele(), self.sequencepair_object.get_secondaryallele()]:
				distribution_split = self.split_cag_target(allele.get_fwarray())
				target_distro = distribution_split['CCG{}'.format(allele.get_ccg())]
				allele.set_cagthreshold(0.50)
				fod_failstate, cag_indexes = self.peak_detection(allele, target_distro, 1, 'CAGHet')
				while fod_failstate:
					fod_failstate, cag_indexes = self.peak_detection(allele, target_distro, 1, 'CAGHet', fod_recall=True)
				allele.set_fodcag(cag_indexes)

		########################
		## Homozygous for CCG ##
		########################
		if self.zygosity_state == 'HOMO':

			##
			## Double check CCG matches.. be paranoid
			primary_ccg = self.sequencepair_object.get_primaryallele().get_ccg()
			secondary_ccg = self.sequencepair_object.get_secondaryallele().get_ccg()
			if not primary_ccg == secondary_ccg:
				raise Exception('CCG Mismatch, in CAG calling. Fatal error.')

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
			for allele in [self.sequencepair_object.get_primaryallele(), self.sequencepair_object.get_secondaryallele()]:
				distribution_split = self.split_cag_target(allele.get_fwarray())
				target_distro = distribution_split['CCG{}'.format(allele.get_ccg())]
				allele.set_cagthreshold(0.50)

				fod_failstate, cag_indexes = self.peak_detection(allele, target_distro, distance_threshold, 'CAGHom', est_dist=estimated_distance)
				while fod_failstate:
					fod_failstate, cag_indexes = self.peak_detection(allele, target_distro, distance_threshold, 'CAGHom', est_dist=estimated_distance, fod_recall=True)
				allele.set_fodcag(cag_indexes)

		return pass_gtp

	def genotype_validation(self):

		##
		## Constructs
		pass_vld = True
		primary_allele = self.sequencepair_object.get_primaryallele()
		secondary_allele = self.sequencepair_object.get_secondaryallele()
		distribution_split = self.split_cag_target(primary_allele.get_fwarray())
		ccg_zygstate = self.sequencepair_object.get_ccgzygstate()

		##
		## Primary Allele
		primary_dsp_ccg = primary_allele.get_ccg(); primary_fod_ccg = primary_allele.get_fodccg()
		primary_dsp_cag = primary_allele.get_cag(); primary_fod_cag = primary_allele.get_fodcag()

		##
		## Secondary Allele
		secondary_dsp_ccg = secondary_allele.get_ccg(); secondary_fod_ccg = secondary_allele.get_fodccg()
		secondary_dsp_cag = secondary_allele.get_cag(); secondary_fod_cag = secondary_allele.get_fodcag()

		##
		## Double check fod peaks
		def dimension_checker(input_list):

			fod = input_list[0]
			dsp = input_list[1]
			allele = input_list[2]

			for i in range(0, len(fod)):
				if np.isclose([fod[i]], [dsp], atol=1.0):
					allele.set_fodcag(fod[i])

		for item in [[primary_fod_cag, primary_dsp_cag, primary_allele], [secondary_fod_cag, secondary_dsp_cag, secondary_allele]]:
			dimension_checker(item)
			primary_fod_cag = primary_allele.get_fodcag(); secondary_fod_cag = secondary_allele.get_fodcag()
		##
		## Check for potential homozygous haplotype/neighbouring peak
		if ccg_zygstate == 'HOMO' and abs(primary_dsp_cag-secondary_dsp_cag) == 1:
			primary_target = distribution_split['CCG{}'.format(primary_allele.get_ccg())]
			secondary_target = distribution_split['CCG{}'.format(secondary_allele.get_ccg())]
			primary_reads = primary_target[primary_allele.get_cag()-1]
			secondary_reads = secondary_target[secondary_allele.get_cag()-1]
			diff = abs(primary_reads-secondary_reads)
			pcnt = (diff/max([primary_reads, secondary_reads]))*100
			if pcnt <= 0.20:
				self.sequencepair_object.set_neighbouringpeaks(True)
				return pass_vld
			elif primary_fod_cag.all() == secondary_fod_cag.all():
				self.sequencepair_object.set_homozygoushaplotype(True)
				return pass_vld

		##
		## Double check zygosity..
		if not (primary_fod_ccg == secondary_fod_ccg) and ccg_zygstate == 'HOMO':
			raise Exception('CCG validity check failure')
		if (primary_fod_ccg == secondary_fod_ccg) and ccg_zygstate == 'HETERO':
			raise Exception('CCG validity check failure')

		def read_comparison(val1, val2):
			if np.isclose(val1,val2,atol=1):
				return val2
			else:
				return val1

		##
		## Ensure integrity
		if not primary_dsp_ccg == int(primary_fod_ccg):
			if read_comparison(primary_dsp_ccg, int(primary_fod_ccg)) == primary_fod_ccg:
				pass_vld = True
			else:
				pass_vld = False

		if not primary_dsp_cag == int(primary_fod_cag):
			if read_comparison(primary_dsp_cag, int(primary_fod_cag)) == primary_fod_cag:
				pass_vld = True
			else:
				pass_vld = False

		if not secondary_dsp_ccg == int(secondary_fod_ccg):
			if read_comparison(secondary_dsp_ccg, int(secondary_fod_ccg)) == secondary_fod_ccg:
				pass_vld = True
			else:
				pass_vld = False

		if not secondary_dsp_cag == int(secondary_fod_cag):
			if read_comparison(secondary_dsp_cag, int(secondary_fod_cag)) == secondary_fod_cag:
				pass_vld = True
			else:
				pass_vld = False

		return pass_vld

	def inspect_peaks(self):

		##
		## Constructs
		pass_vld = True

		for allele in [self.sequencepair_object.get_primaryallele(), self.sequencepair_object.get_secondaryallele()]:
			distribution_split = self.split_cag_target(allele.get_fwarray())
			target = distribution_split['CCG{}'.format(allele.get_ccg())]
			linspace = np.linspace(0,199,200)

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
			## Calc SNR ratio
			##TODO rework.. density estimation?
			# snr_ratio = stats.signaltonoise(target)
			# if snr_ratio > 0.5: pass_vld = False
			# else: allele.set_signaltonoise(snr_ratio)

			##
			## Calculate peak dropoff
			nminus = target[allele.get_cag()-2]; n = target[allele.get_cag()-1]; nplus = target[allele.get_cag()]
			nminus_overn = nminus/n; nplus_overn = nplus/n
			dropoff_list = [nminus_overn, nplus_overn]

			if not self.sequencepair_object.get_homozygoushaplotype() or self.sequencepair_object.get_neighbouringpeaks():
				allele.set_immediate_dropoff(dropoff_list)
				## inform user
				if np.isclose(nminus_overn, [0.25], atol=0.02):
					allele.set_nminuswarninglevel(1)
				if np.isclose(nplus_overn, [0.05], atol=0.02):
					allele.set_nminuswarninglevel(1)
				## warn user
				if np.isclose(nminus_overn, [0.35], atol=0.04):
					allele.set_nminuswarninglevel(2)
				if np.isclose(nplus_overn, [0.15], atol=0.03):
					allele.set_nminuswarninglevel(2)
				## severe warning
				if np.isclose(nminus_overn, [0.45], atol=0.05):
					allele.set_nminuswarninglevel(3)
				if np.isclose(nplus_overn, [0.27], atol=0.03):
					allele.set_npluswarninglevel(3)
				## fail sample
				if np.isclose(nminus_overn, [0.65], atol=0.05):
					allele.set_nminuswarninglevel(4)
					pass_vld = False
				if np.isclose(nplus_overn, [0.45], atol=0.05):
					allele.set_nminuswarninglevel(4)
					pass_vld = False
			else:
				allele.set_nminuswarninglevel(2)
				allele.set_npluswarninglevel(2)

			##
			## Somatic mosaicism
			#TODO calculate (N+1 to N+10) over N

			##
			## Slippage
			#TODO caluclate (N-3 to N-1) over N

			##
			## If we get here; alleles are valid
			allele.set_ccgvalid(True)
			allele.set_cagvalid(True)
			allele.set_genotypestatus(True)

		return pass_vld

	def render_graphs(self):

		def graph_subfunction(x, y, axis_labels, xticks, peak_index, predict_path, file_handle):
			x = np.linspace(x[0],x[1],x[2])
			plt.figure(figsize=(10, 6)); plt.title(self.sequencepair_object.get_label())
			plt.xlabel(axis_labels[0]); plt.ylabel(axis_labels[1])
			plt.xticks(np.arange(xticks[0][0], xticks[0][1], xticks[0][2])); plt.xlim(xticks[1][0], xticks[1][1])
			pplot(x, y, peak_index)
			peak_index = [i+1 for i in peak_index]
			plt.legend(['Genotype: {}'.format(peak_index)])
			plt.savefig(os.path.join(predict_path, file_handle), format='pdf')
			plt.close()

		##
		## CCG heterozygous example
		## i.e. CCG two peaks, one CAG dist per peak
		if self.zygosity_state == 'HETERO':

			##
			## Render CCG graph, append path to allele path list
			hetero_graphs = []
			ccg_peaks = [int(self.sequencepair_object.get_primaryallele().get_fodccg() - 1),
						 int(self.sequencepair_object.get_secondaryallele().get_fodccg() - 1)]
			concat = np.asarray([a + b for a, b in zip(self.sequencepair_object.get_primaryallele().get_rvarray(),
													   self.sequencepair_object.get_secondaryallele().get_rvarray())])
			graph_subfunction([0, 21, 20], concat,
							  ['CCG Value', 'Read Count'], ([1, 20, 1], [1, 20]),
							  ccg_peaks, self.sequencepair_object.get_predictpath(), 'CCGDetection.pdf')
			hetero_graphs.append(os.path.join(self.sequencepair_object.get_predictpath(), 'CCGDetection.pdf'))
			plt.close()

			##
			## For each CCG allele in this heterozygous sample
			for allele in [self.sequencepair_object.get_primaryallele(),
						   self.sequencepair_object.get_secondaryallele()]:
				distribution_split = self.split_cag_target(allele.get_fwarray())
				target_distro = distribution_split['CCG{}'.format(allele.get_ccg())]
				graph_subfunction([0, 199, 200], target_distro, ['CAG Value', 'Read Count'],
								  ([1, 200, 50], [1, 200]), [allele.get_fodcag() - 1],
								  self.sequencepair_object.get_predictpath(),
								  'CCG{}-CAGDetection.pdf'.format(allele.get_fodccg()))
				hetero_graphs.append(os.path.join(self.sequencepair_object.get_predictpath(),
												  'CCG{}-CAGDetection.pdf'.format(allele.get_fodccg())))
				plt.close()
			self.sequencepair_object.get_primaryallele().set_allelegraphs(hetero_graphs)
			self.sequencepair_object.get_secondaryallele().set_allelegraphs(hetero_graphs)

		##
		## CCG homozygous example
		## i.e. CCG one peak, one CAG dist per peak
		if self.zygosity_state == 'HOMO':
			homo_graphs = []
			ccg_peaks = [int(self.sequencepair_object.get_primaryallele().get_fodccg() - 1),
						 int(self.sequencepair_object.get_secondaryallele().get_fodccg() - 1)]
			cag_peaks = [int(self.sequencepair_object.get_primaryallele().get_fodcag() - 1),
						 int(self.sequencepair_object.get_secondaryallele().get_fodcag() - 1)]
			distribution_split = self.split_cag_target(self.sequencepair_object.get_primaryallele().get_fwarray())
			target_distro = distribution_split['CCG{}'.format(self.sequencepair_object.get_primaryallele().get_ccg())]
			graph_subfunction([0, 21, 20], self.sequencepair_object.get_primaryallele().get_rvarray(),
							  ['CCG Value', 'Read Count'], ([1, 20, 1], [1, 20]), ccg_peaks,
							  self.sequencepair_object.get_predictpath(), 'CCGDetection.pdf')
			plt.close()
			graph_subfunction([0, 199, 200], target_distro, ['CAG Value', 'Read Count'],
							  ([1, 200, 50], [1, 200]), cag_peaks, self.sequencepair_object.get_predictpath(),
							  'CCG{}-CAGDetection.pdf'.format(
								  self.sequencepair_object.get_primaryallele().get_fodccg()))
			plt.close()
			homo_graphs.append(os.path.join(self.sequencepair_object.get_predictpath(), 'CCGDetection.pdf'))
			homo_graphs.append(os.path.join(self.sequencepair_object.get_predictpath(), 'CCG{}-CAGDetection.pdf'.format(
				self.sequencepair_object.get_primaryallele().get_fodccg())))
			self.sequencepair_object.get_primaryallele().set_allelegraphs(homo_graphs)
			self.sequencepair_object.get_secondaryallele().set_allelegraphs(homo_graphs)

		##
		## Merge graphs into a single summary PDF
		sample_pdf_path = os.path.join(self.sequencepair_object.get_predictpath(), 'SampleSummary.pdf')
		c = canvas.Canvas(sample_pdf_path, pagesize=(500,250))
		header_string = '{}{}'.format('Sample header: ', self.sequencepair_object.get_label())
		primary_string = '{}({},{}) ({})'.format('Primary: ', self.sequencepair_object.get_primaryallele().get_fodcag(),
												 self.sequencepair_object.get_primaryallele().get_fodccg(),
												 self.sequencepair_object.get_primaryallele().get_allelestatus())
		secondary_string = '{}({},{}) ({})'.format('Secondary: ', self.sequencepair_object.get_secondaryallele().get_fodcag(),
												   self.sequencepair_object.get_secondaryallele().get_fodccg(),
												   self.sequencepair_object.get_secondaryallele().get_allelestatus())
		##
		## Set font colour
		## If invalid data -- atypical allele but no re-alignment, bad
		## If valid data -- atypical allele but re-alignment, so.. ok
		if self.invalid_data:
			c.setFillColorRGB(255, 0, 0)
			c.drawCentredString(250, 50, '!! Atypical alleles without re-alignment !!')
		if not self.invalid_data:
			c.setFillColorRGB(0, 0, 0)
		c.drawCentredString(250, 150, header_string)
		c.drawCentredString(250, 125, primary_string)
		c.drawCentredString(250, 100, secondary_string)
		c.save()

		##
		## Allele graphs
		## Ensure uniqueness of entries in primary/secondary (i.e. no duplicating CCG graph)
		primary_graphs = self.sequencepair_object.get_primaryallele().get_allelegraphs()[0]
		secondary_graphs = self.sequencepair_object.get_primaryallele().get_allelegraphs()[0]
		sample_graphs = primary_graphs + secondary_graphs
		sample_uniques = list(set(sample_graphs))
		target_pdflist = [sample_pdf_path] + sample_uniques

		##
		## Merge alleles together
		merger = PyPDF2.PdfFileMerger()
		for pdf in target_pdflist:
			merger.append(pdf)
		merger.write(sample_pdf_path)
		merger.close()

	def calculate_score(self):

		##
		## For both alleles
		for allele in [self.sequencepair_object.get_primaryallele(), self.sequencepair_object.get_secondaryallele()]:

			##
			## Start score high, deduct for questionable calls..
			allele_confidence = 100

			##
			## Sample based genotyping flags
			if self.sequencepair_object.get_recallcount() == 7: allele_confidence -= 25
			if 7 > self.sequencepair_object.get_recallcount() > 4: allele_confidence -= 15
			if 4 > self.sequencepair_object.get_recallcount() > 0: allele_confidence -= 5
			else: allele_confidence += 10

			if self.sequencepair_object.get_homozygoushaplotype(): allele_confidence -= 15
			elif self.sequencepair_object.get_neighbouringpeaks(): allele_confidence -= 25
			else: allele_confidence += 20

			##
			## Allele based genotyping flags
			## Allele typical/atypical structure
			if allele.get_allelestatus() == 'Atypical':
				allele_confidence -= 5
				if np.isclose([float(allele.get_atypicalpcnt())],[50.00],atol=5.00):
					allele_confidence -= 20
				if np.isclose([float(allele.get_atypicalpcnt())],[80.00],atol=20.00):
					allele_confidence += 15
			if allele.get_allelestatus() == 'Typical':
				allele_confidence += 10
				if np.isclose([float(allele.get_typicalpcnt())],[50.00],atol=5.00):
					allele_confidence -= 20
				if np.isclose([float(allele.get_typicalpcnt())],[80.00],atol=20.00):
					allele_confidence += 15

			##
			## Total reads in sample..
			if allele.get_totalreads() > 15000:
				allele_confidence += 5
			if allele.get_totalreads() < 3000:
				allele_confidence -= 10

			##
			## Peak Interpolation
			if allele.get_interpolation_warning():
				allele_confidence -= 5
				if 2.00 > allele.get_interpdistance() > 0.00:
					allele_confidence -= 10

			##
			## Signal to noise
			## TODO

			##
			## Slippage and somatic mosaicism
			## TODO

			##
			## Diminished peaks
			## TODO

			##
			## Peak calling thresholds
			for contig in [allele.get_ccgthreshold(), allele.get_cagthreshold()]:
				if contig != 0.5:
					if 0.5 > contig > 0.3: allele_confidence -= 5
					if 0.3 > contig > 0.0: allele_confidence -= 10
				else: allele_confidence += 10

			##
			## Peak dropoff warnings
			for peak_position_error in [allele.get_nminuswarninglevel(), allele.get_npluswarninglevel()]:
				if peak_position_error == 1: allele_confidence -= 5
				if 2 >= peak_position_error > 1: allele_confidence -= 10
				else: allele_confidence -= 15

			allele.set_alleleconfidence(sorted([0, allele_confidence, 100])[1])
			allele.set_allelegenotype('{}_{}_{}_{}_{}'.format(allele.get_fodcag(), allele.get_caacag(),
															  allele.get_ccgcca(), allele.get_fodccg(), allele.get_cct()))

	def set_report(self):

		for allele in [self.sequencepair_object.get_primaryallele(), self.sequencepair_object.get_secondaryallele()]:

			##
			## Report path for this allele
			allele_filestring = '{}{}'.format(allele.get_reflabel(), '_AlleleReport.txt')
			report_path = os.path.join(self.sequencepair_object.get_predictpath(), allele_filestring)
			allele.set_allelereport(report_path)
			report_string = '{}{}\n\n{}\n{}{}\n{}{}\n{}{}\n{}' \
							'{}\n{}{}\n{}{}\n\n{}\n{}{}\n{}{}\n' \
							'{}{}\n{}{}\n{}{}\n{}{}\n{}{}\n' \
							'{}{}\n\n{}\n{}{}\n{}{}\n{}{}\n{}' \
							'{}\n{}{}\n{}{}'.format('Allele Report>> ', self.sequencepair_object.get_label(),
													'Summary Information>>',
													'Genotype: ', allele.get_allelegenotype(),
													'Confidence: ', allele.get_alleleconfidence(),
													'Structure Status: ', allele.get_allelestatus(),
													'Typical Pcnt: ', allele.get_typicalpcnt(),
													'Atypical Pcnt: ', allele.get_atypicalpcnt(),
													'Total Reads: ', allele.get_totalreads(),
													'Flags>>',
													'Recall Count: ', self.sequencepair_object.get_recallcount(),
													'Homozygous Haplotype: ', self.sequencepair_object.get_homozygoushaplotype(),
													'Neighbouring Peaks: ', self.sequencepair_object.get_neighbouringpeaks(),
													'Diminished Peaks: ', 'TODO',
													'Backwards Slippage: ', 'TODO',
													'Somatic Mosaicism: ', 'TODO',
													'Peak Interpolation Warning: ', allele.get_interpolation_warning(),
													'Peak Interpolation Distance: ', allele.get_interpdistance(),
													'Data Quality>>',
													'Signal To Noise: ', 'TODO',
													'Immediate Dropoffs: ', allele.get_immediate_dropoff(),
													'N-1 Warning Level: ', allele.get_nminuswarninglevel(),
													'N+1 Warning Level: ', allele.get_npluswarninglevel(),
													'CCG Threshold: ', allele.get_ccgthreshold(),
													'CAG Threshold: ', allele.get_cagthreshold())

			##
			## Write to file
			with open(report_path, 'w') as outfi:
				outfi.write(report_string)
				outfi.close()

	def get_report(self):

		self.allele_report = [self.sequencepair_object.get_primaryallele().get_allelereport(),
							  self.sequencepair_object.get_secondaryallele().get_allelereport()]
		return self.allele_report