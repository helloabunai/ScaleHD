from __future__ import division

#/usr/bin/python
__version__ = 0.01
__author__ = 'alastair.maxwell@glasgow.ac.uk'

##
## Generic imports
import sys
import os
import csv
import peakutils
import numpy as np
import logging as log
from collections import Counter
from sklearn import svm
from sklearn import preprocessing
from sklearn.multiclass import OneVsOneClassifier
from peakutils.plot import plot as pplot
from matplotlib import pyplot
import pandas as pd

##
## Backend Junk
from ..__backend import Colour as clr
from ..__backend import DataLoader

class GenotypePrediction:
	def __init__(self, data_pair, prediction_path, training_data, instance_params):

		"""
		Prediction stage of the pipeline -- use of SVM, density estimation and First Order Differentials
		to automate calling of a sample's genotype based on data information derived from alignment read counts.
		Utilises forward reads for CAG information, reverse reads for CCG information

		General workflow of the process:
		--Take reverse reads, aggregate every CAG count for each CCG
		--Use unlabelled sample into CCG zygosity SVM for het/hom prediction
		--Data cleaning/normalisation/etc
		--Two pass algorithm searches over CCG distribution with SVM result in mind
		--Density estimation for peak estimate/location estimate/clarity estimate
		--F.O.D. for more specific peak calling in particular distribution
		--Get appropriate CAG distribution for the current CCG distribution
		--Repeat two-pass algorithm for this distribution
		--Return allele
		--Repeat for second allele
		--Combine and return genotype

		Work in progress
		Things to improve:
		--Error flagging for manual inspection
		--Alternate performant mechanisms?
		"""

		##
		## Paths to files/data required
		self.data_pair = data_pair
		self.prediction_path = prediction_path
		self.training_data = training_data
		self.instance_params = instance_params

		##
		## Arrays of the current sample's CCG (aggregated)
		## Used in SVM classifier for predicting unlabelled state
		self.classifier, self.encoder = self.build_zygosity_model()

		##
		## Confidence flags (subject to change)
		## Flag meanings:
		## --CCGZyg_disconnect: between forward/reverse zygosity, a difference is recorded.
		## 		We only care about reverse read zygosity, but it's still a good QOL metric
		## --CCGExpansion_skew: Of the 'major' (first) peak, n-1 was > in value than the 'minor' peak.
		## 		Not explicitly bad, as sometimes mosaicism means n-1 is >> expanded peak. Innocent informative warning.
		## --CCGPeak_ambiguous: Similar to density_ambiguous;; specifically lots of low density values NEXT to a peak
		## 		For an automated system we want clean peaks, so inform the user when this is not the case
		## --CCGDensity_ambiguous: A lot of similar (low) densities were detected in the distribution.
		## 		In our data, the sparsest value == peak. Lots of low values == messy data, which is bad.
		## --CCGRecall_warning: first/second passes on CCG genotyping mismatched, had to re-call with lower threshold
		## 		Peaks differed between first/second pass, accuracy may be an issue for this particular sample.
		## --CCGPeak_oob: more than two peaks were returned from the FOD function
		## 		If this happens, something fucked up.
		## --CAGConsensusSpread_warning: Forward reads were very broad peaks; consensus CAG from FW+RV used
		##  	If this is used, quality of prediction is going to be very poor. Brute force last ditch attempt at genotyping.
		## --FPSP_Disconnect: Prediction between first and second pass of algorithm differ
		## 		If this happens, minor warning. Expected to happen occasionally. Doesn't hurt accuracy.

		self.prediction_confidence = 0
		self.cag_intermediate = [0,0]
		self.genotype_flags = {'Primary_Allele':[0,0],
							   'Secondary_Allele':[0,0],
							   'PeakOverPopulation':False,
							   'CCGZyg_disconnect':False,
							   'CCGExpansion_skew':False,
							   'CCGPeak_ambiguous':False,
							   'CCGDensity_ambiguous':False,
							   'CCGRecall_warning':False,
							   'CCGPeak_oob':False,
							   'CAGRecall_warning':False,
							   'CAGConsensusSpread_warning':False,
							   'FPSP_Disconnect':False,
							   'Primary_Mosaicism':[],
							   'Secondary_Mosaicism':[]}

		##
		## Unlabelled distributions to utilise for prediction
		self.forward_distribution = self.scrape_distro(self.data_pair[0])
		self.reverse_distribution = self.scrape_distro(self.data_pair[1])

		##
		## !! STAGE ONE !!
		## Determine Zygosity of input distributions
		self.forwardccg_aggregate = self.distribution_collapse(self.forward_distribution)
		self.reverseccg_aggregate = self.distribution_collapse(self.reverse_distribution)
		self.zygosity_state = self.predict_zygstate()

		##
		## !! STAGE TWO !!
		## Determine CCG peak(s)/genotype(s)
		## Call generic two-pass algorithm: Density estimation/First order differentials
		ccg_pass, ccg_genotype = self.determine_ccg_genotype()
		while not ccg_pass:
			self.genotype_flags['CCGRecall_warning'] = True
			ccg_pass, ccg_genotype = self.determine_ccg_genotype(threshold_bias=True)

		self.genotype_flags['Primary_Allele'][1] = ccg_genotype[0]
		self.genotype_flags['Secondary_Allele'][1] = ccg_genotype[1]

		##
		## !! STAGE THREE !!
		## Now we have CCG identified, we can get the relevant CAG distribution for that CCG
		## And utilise the same generic functions to call the peaks for CAG
		## Once combined, we can call our genotype
		cag_pass, cag_genotype = self.determine_cag_genotype()
		while not cag_pass:
			self.genotype_flags['CAGRecall_warning'] = True
			cag_pass, cag_genotype = self.determine_cag_genotype(threshold_bias=True)

		if cag_genotype is ['err','err']:
			raise Exception

		self.genotype_flags['Primary_Allele'][0] = cag_genotype[0]
		self.genotype_flags['Secondary_Allele'][0] = cag_genotype[1]

		##
		## !! STAGE FOUR !!
		## Do simple somatic mosaicism calculations
		## Then append all information to the report to be returned to shd.sherpa
		self.genotype_flags['Primary_Mosaicism'] = self.somatic_calculations(self.genotype_flags['Primary_Allele'])
		self.genotype_flags['Secondary_Mosaicism'] = self.somatic_calculations(self.genotype_flags['Secondary_Allele'])
		self.gtype_report = self.generate_report()

	def build_zygosity_model(self):

		##
		## Classifier object to be wrapped into One vs One classifier
		svc_object = svm.LinearSVC(C=1.0, loss='squared_hinge', penalty='l2',
								   dual=False, tol=1e-4, multi_class='crammer_singer',
								   fit_intercept=True, intercept_scaling=1, verbose=0,
								   random_state=0, max_iter=-1)

		##
		## Take raw training data into object
		trainingdat_collapsed_ccg = self.training_data['CollapsedCCGZygosity']
		trainingdat_descriptionfi = self.training_data['GenericDescriptor']
		training_model = DataLoader(trainingdat_collapsed_ccg, trainingdat_descriptionfi).load_model()

		##
		## Model data taken into array, fit to classifier
		## return objects
		X = preprocessing.normalize(training_model.DATA)
		Y = training_model.TARGET
		ovo_svc = OneVsOneClassifier(svc_object).fit(X,Y)
		encoder = training_model.ENCDR

		return ovo_svc, encoder

	@staticmethod
	def scrape_distro(distribution_file):

		##
		## Scrapes the read distribution from the input file
		placeholder_array = []
		with open(distribution_file) as csv_file:
			source = csv.reader(csv_file, delimiter=',')
			next(source)  # skipping header
			for row in source:
				placeholder_array.append(int(row[2]))
			csv_file.close()

		novel_distro = np.array(placeholder_array)
		return novel_distro

	@staticmethod
	def distribution_collapse(distribution_array):

		##
		## Assume there is 20 CCG
		try:
			ccg_arrays = np.split(distribution_array, 20)
		except ValueError:
			log.critical('{}{}{}{}'.format(clr.red,'shd__ ',clr.end,'Repeat distribution not evenly divided by 20. Aligned to non CAG1-200/CCG1-20 reference?'))
			sys.exit(2)

		##
		## Sum the values for each ccg
		ccg_counter = 1
		collapsed_array = []
		for ccg_array in ccg_arrays:
			collapsed_array.append(np.sum(ccg_array))
			ccg_counter+=1

		return np.asarray(collapsed_array)

	def predict_zygstate(self):

		##
		## Reshape the input distributions so skl doesn't complain about 1D
		## normalise too, as in tandem with training data -- requires casting to float64
		forward_reshape = preprocessing.normalize(np.float64(self.forwardccg_aggregate.reshape(1,-1)))
		reverse_reshape = preprocessing.normalize(np.float64(self.reverseccg_aggregate.reshape(1,-1)))

		##
		## Predict the zygstate, then decode from hash into literal label value
		forward_zygstate = str(self.encoder.inverse_transform(self.classifier.predict(forward_reshape)))
		reverse_zygstate = str(self.encoder.inverse_transform(self.classifier.predict(reverse_reshape)))

		##
		## We only really care about the zygosity from reverse reads, but for a QOL metric
		## we can take the value from the forward reads in addition -- if match, great
		## if disconnect, fine, but set bool to be used in confidence algorithm later on..
		if not forward_zygstate == reverse_zygstate:
			self.genotype_flags['CCGZyg_disconnect'] = True
			return reverse_zygstate[2:-2]
		else:
			self.genotype_flags['CCGZyg_disconnect'] = False
			return reverse_zygstate[2:-2]

	def update_flags(self, updated_flags):

		"""
		Compares keys between instance-wide self.genotype_flags and updated values
		which can originate from various methods in the pipeline
		When a key is matched, the value is updated to represent the current state of play
		"""

		for updated_key, updated_value in updated_flags.iteritems():
			for initial_key, initial_value in self.genotype_flags.iteritems():
				if initial_key == updated_key:
					self.genotype_flags[initial_key] = updated_value

	def determine_ccg_genotype(self, genotype_pass=True, threshold_bias=False):

		"""
		Workflow director for the CCG genotype stage
		Fill this out..
		"""

		target_peak_count = 0
		if self.zygosity_state == 'HOMO':
			target_peak_count = 1
		if self.zygosity_state == 'HETERO':
			target_peak_count = 2

		##
		## Create object for two-pass algorithm
		graph_parameters = [20,'CCGDensityEstimation.png','CCG Density Distribution',['Read Count', 'Bin Density']]
		ccg_inspector = PredictionTwoPass(prediction_path=self.prediction_path,
										  input_distribution=self.reverseccg_aggregate,
										  target_peak_count=target_peak_count,
										  graph_parameters=graph_parameters)

		##
		## Density, first pass
		## Update instance dictionary of any error flags returned from density stage
		first_pass = ccg_inspector.density_estimation(plot_flag=True)
		density_warnings = ccg_inspector.get_warnings()
		self.update_flags(density_warnings)


		##
		## First order differentials
		## Update instance dictionary of any error flags returned from FOD stage
		fod_parameters = [[0,19,20],'CCG Peaks',['CCG Value','Read Count'], 'CCGPeakDetection.png']
		try:
			second_pass = ccg_inspector.differential_peaks(first_pass, fod_parameters, threshold_bias)
		except ValueError:
			return False, [0,0]

		fod_warnings = ccg_inspector.get_warnings()
		self.update_flags(fod_warnings)

		##
		## Check first pass results == second pass results
		## if not, genotype_pass = False
		first_pass_estimate = [first_pass['PrimaryPeak'],first_pass['SecondaryPeak']]
		second_pass_estimate = [second_pass['PrimaryPeak'],second_pass['SecondaryPeak']]

		if len(second_pass_estimate)>len(first_pass_estimate): self.genotype_flags['CCGTriplet_warning'] = True
		if not first_pass_estimate == second_pass_estimate or len(second_pass_estimate)>len(first_pass_estimate):
			genotype_pass = False

		##
		## Return when done
		return genotype_pass, second_pass_estimate

	@staticmethod
	def split_cag_target(input_distribution, ccg_target):

		##
		## We need to take the relevant information from forward HD sequence
		## as it is better quality for the target CAG repeat region
		## Split the entire distribution per sample into contigs for each CCG (4000 -> 200*20)

		cag_split = [input_distribution[i:i + 200] for i in xrange(0, len(input_distribution), 200)]
		distribution_dictionary = {}
		for i in range(0,len(cag_split)):
			distribution_dictionary['CCG'+str(i+1)] = cag_split[i]

		current_target_distribution = distribution_dictionary['CCG' + str(ccg_target)]
		return current_target_distribution

	def determine_cag_genotype(self, genotype_pass=True, threshold_bias=False):

		##
		## Target peak count
		## We're looking at CAG distributions for a specific CAG..
		## If CCG is homozygous (i.e., one CCG distribution), there will be two CAG peaks on that single distribution
		## If CCG is heterozygous (i.e., two CCG distribtuions), there will be one peak on each CCG distribution
		target_peak_count = 0
		target_distribution = {}
		if self.zygosity_state == 'HOMO':
			target_peak_count = 2
			cag_target = self.split_cag_target(self.forward_distribution, self.genotype_flags['Primary_Allele'][1])
			target_distribution[self.genotype_flags['Primary_Allele'][1]] = cag_target
		if self.zygosity_state == 'HETERO':
			target_peak_count = 1
			cag_target_major = self.split_cag_target(self.forward_distribution, self.genotype_flags['Primary_Allele'][1])
			cag_target_minor = self.split_cag_target(self.forward_distribution, self.genotype_flags['Secondary_Allele'][1])
			target_distribution[self.genotype_flags['Primary_Allele'][1]] = cag_target_major
			target_distribution[self.genotype_flags['Secondary_Allele'][1]] = cag_target_minor

		##
		## Iterate over distributions that we are looking at
		for cag_key, distro_value in target_distribution.iteritems():

			graph_parameters = [20, 'CAG'+str(cag_key)+'DensityEstimation.png','CAG Density Distribution',['Read Count','Bin Density']]
			cag_inspector = PredictionTwoPass(prediction_path=self.prediction_path,
											  input_distribution=distro_value,
											  target_peak_count=target_peak_count,
											  graph_parameters=graph_parameters)

			##
			## Check distribution spread
			## Pre-stage to 2-pass; CAG distributions are longer so spread is more likely
			## Check quality of peaks before progressing; if super-poor quality,
			## Combine FW+RV for this CCG; attempt brute force genotyping on consensus sequence
			pre_pass = cag_inspector.investigate_spread()
			if pre_pass:
				self.genotype_flags['CAGConsensusSpread_warning'] = True
				bruteforce_consensus = cag_inspector.generate_consensus(self.split_cag_target(self.forward_distribution, cag_key),
																		self.split_cag_target(self.reverse_distribution, cag_key))
				cag_inspector.set_distribution(bruteforce_consensus)

			##
			## Density, first pass
			first_pass = cag_inspector.density_estimation(plot_flag=False)

			if first_pass['PeakThreshold'] is None:
				log.error('{}{}{}{}'.format(clr.red,'shd__ ',clr.end,'Density Estimation failed. Alignment must be worse than low quality.'))
				return True, ['err','err']

			##
			## FOD, second pass
			fod_parameters = [[0,199,200],'CAG'+str(cag_key)+' Peaks',['CAG Value','Read Count'],'CAG'+str(cag_key)+'PeakDetection.png']
			try:
				second_pass = cag_inspector.differential_peaks(first_pass, fod_parameters, threshold_bias)
			except ValueError:
				return False, [0,0]

			##
			## Concat results
			first_pass_estimate = [first_pass['PrimaryPeak'],first_pass['SecondaryPeak']]
			second_pass_estimate = [second_pass['PrimaryPeak'],second_pass['SecondaryPeak']]
			if not first_pass_estimate == second_pass_estimate or len(second_pass_estimate)>len(first_pass_estimate):
				genotype_pass = False

			##
			## Ensure the right CAG is assigned to the right CCG
			if self.zygosity_state == 'HOMO':
				if cag_key == self.genotype_flags['Primary_Allele'][1]:
					self.cag_intermediate[0] = second_pass_estimate[0]
					self.cag_intermediate[1] = second_pass_estimate[1]
			if self.zygosity_state == 'HETERO':
				if cag_key == self.genotype_flags['Primary_Allele'][1]:
					self.cag_intermediate[0] = second_pass_estimate[0]
				if cag_key == self.genotype_flags['Secondary_Allele'][1]:
					self.cag_intermediate[1] = second_pass_estimate[0]

		cag_genotype = [self.cag_intermediate[0],self.cag_intermediate[1]]
		return genotype_pass, cag_genotype

	def somatic_calculations(self, genotype):

		"""
		Function to acquire N-1/N/N+1 distribution values for the acquired genotype
		Calculations carried out on these values and returned for report append
		"""

		##
		## Create object of mosaicism investigator to begin calculation prep
		## Takes raw 200x20 distribution and slices into 20 arrays of 200
		## Orders into a dataframe object with CCG<val> labels
		## Scrapes appropriate N-1/N/N+1 based on genotype for this instance
		mosaicism_object = MosaicismInvestigator(genotype, self.forward_distribution)
		ccg_slices = mosaicism_object.chunks(200)
		ccg_ordered = mosaicism_object.arrange_chunks(ccg_slices)
		allele_values = mosaicism_object.get_nvals(ccg_ordered, genotype)

		##
		## Now we have the values, just calculate and return
		allele_calculations = mosaicism_object.calculate_mosaicism(allele_values)

		##
		## Generate a padded/aligned distribution for the current instance
		## so that for each sample in the report, N will be at the same position
		padded_distribution = mosaicism_object.distribution_padder(ccg_ordered, genotype)

		##
		## Combine calculation dict and distribution into object
		## Return object
		mosaicism_values = [allele_calculations, padded_distribution]
		return mosaicism_values

	def generate_report(self):

		##
		## Mini function to generate the report to be written in sample folder
		## And later utilised in the 'instance wide' report

		report = [self.genotype_flags['Primary_Allele'],
				  self.genotype_flags['Secondary_Allele'],
				  self.genotype_flags['PeakOverPopulation'],
				  self.genotype_flags['CCGZyg_disconnect'],
				  self.genotype_flags['CCGExpansion_skew'],
				  self.genotype_flags['CCGPeak_ambiguous'],
				  self.genotype_flags['CCGDensity_ambiguous'],
				  self.genotype_flags['CCGRecall_warning'],
				  self.genotype_flags['CCGPeak_oob'],
				  self.genotype_flags['CAGRecall_warning'],
				  self.genotype_flags['CAGConsensusSpread_warning'],
				  self.genotype_flags['FPSP_Disconnect'],
				  self.genotype_flags['Primary_Mosaicism'],
				  self.genotype_flags['Secondary_Mosaicism']]

		return report

	def get_report(self):

		return self.gtype_report

class PredictionTwoPass:

	def __init__(self, prediction_path, input_distribution, target_peak_count, graph_parameters):
		"""
		Intro spiel for two-pass object
		fill this out..
		"""
		##
		## Variables for the instance of this object
		self.prediction_path = prediction_path
		self.input_distribution = input_distribution
		self.target_peak_count = target_peak_count
		self.bin_count = graph_parameters[0]
		self.filename = graph_parameters[1]
		self.graph_title = graph_parameters[2]
		self.axes = graph_parameters[3]
		self.instance_parameters = {}

		##
		## Warnings raised during processing of this instance
		self.density_ambiguity = False
		self.expansion_skew = False
		self.peak_ambiguous = False
		self.peak_overpopulated = False

	def histogram_generator(self, filename, graph_title, axes, plot_flag):

		density_histo, density_bins = np.histogram(self.input_distribution, bins=self.bin_count, density=True)
		if plot_flag:
			pyplot.figure(figsize=(10,6))
			bin_width = 0.7 * (density_bins[1] - density_bins[0])
			center = (density_bins[:-1] + density_bins[1:]) / 2
			pyplot.title(graph_title)
			pyplot.xlabel(axes[0])
			pyplot.ylabel(axes[1])
			pyplot.bar(center, density_histo, width=bin_width)
			pyplot.savefig(os.path.join(self.prediction_path, filename), format='png')
			pyplot.close()

		density_frequency = Counter(density_histo)
		for key, value in density_frequency.iteritems():
			if not key == np.float64(0.0) and value > 2:
				self.density_ambiguity = True

		return density_histo, density_bins

	@staticmethod
	def drop_calc(distribution, peak_estimate, peak_index):

		"""
		Submethod used to determine % drop around a peak_estimate
		if the drop of nm1/np1 is < 40%, we're not happy it's a clean peak
		"""

		nm1 = peak_index-1
		np1 = peak_index+1
		peak_nm1 = distribution[nm1]
		peak_np1 = distribution[np1]
		nm1_drop = ((peak_estimate - peak_nm1) / peak_estimate) * 100
		np1_drop = ((peak_estimate - peak_np1) / peak_estimate) * 100

		for drop in [nm1_drop, np1_drop]:
			if drop < 40.00:
				return 1
			else:
				return 0

	def investigate_spread(self, inspection_failure=False):

		##
		## Function to investigate how spread a CAG distribution is
		## Used as a pre-screen since CAG distributions are much longer (200 vs 20)
		## And mosaicism/etc is more prevalent in CAG distributions

		max_variance = np.var(self.input_distribution) / max(self.input_distribution)
		if max_variance < 20.00:
			inspection_failure = True

		return inspection_failure

	@staticmethod
	def generate_consensus(forward_dist, reverse_dist):

		##
		## For each position in both fw/rv
		## total up readcount, produce into 'consensus' distr
		consensus = []
		for f, r in zip(forward_dist, reverse_dist):
			c = f + r
			consensus.append(c)
		return np.asarray(consensus)

	def density_estimation(self, plot_flag):

		##
		## Take distro from input (to this object instance)
		## List for indexing functions; major = normal peak; minor = expanded peak
		distro_list = list(self.input_distribution)
		major_estimate = None
		minor_estimate = None
		peak_distance = None
		peak_threshold = None

		##
		## Create dictionary for variables to returned
		estimated_attributes = {'PrimaryPeak': major_estimate,
								'SecondaryPeak': minor_estimate,
								'PeakDistance': peak_distance,
								'PeakThreshold': peak_threshold}

		##
		## Begin estimating density for this distribution
		## Default behaviour is to do everything as if heterozygous distribution
		## If this is not required, only output is tailored (least computational work)
		major_estimate = max(self.input_distribution)
		major_index = distro_list.index(major_estimate)
		minor_estimate = max(n for n in distro_list if n!=major_estimate)
		minor_index = distro_list.index(minor_estimate)

		## Check that major(n-1) is not minor peak
		## Raise warning flag if this is the case
		if minor_index == major_index-1:
			real_minor_estimate = max(n for n in distro_list if n!=major_estimate and n!=minor_estimate)
			real_minor_index = distro_list.index(real_minor_estimate)
			minor_estimate = real_minor_estimate
			minor_index = real_minor_index
			self.expansion_skew = True

		##
		## Density Histogram
		hist, bins = self.histogram_generator(self.filename, self.graph_title, self.axes, plot_flag)
		histo_list = list(hist)

		##
		## Bins is len(hist+1) & 0 indexed.. correct
		major_estimate_bin = np.digitize(major_estimate, bins)-2
		minor_estimate_bin = np.digitize(minor_estimate, bins)-1

		##
		## If two peaks are identical sparsity, set it as such
		## Otherwise, major = most sparse, minor = 2nd most sparse
		if histo_list[major_estimate_bin] == histo_list[minor_estimate_bin]:
			major_estimate_sparsity = min(n for n in hist if n!=0)
			minor_estimate_sparsity = major_estimate_sparsity
		else:
			major_estimate_sparsity = min(n for n in hist if n!=0)
			minor_estimate_sparsity = min(n for n in hist if n!=0 and n!=major_estimate_sparsity)

		##
		## Peak Distance
		peak_distance = np.absolute(major_index-minor_index)

		##
		## Multiple low densities within distribution check
		fuzzy_count = 0
		for density in histo_list:
			if np.isclose(major_estimate_sparsity, density):
				fuzzy_count+=1
		if fuzzy_count > 3:
			self.density_ambiguity = True

		##
		## Peak clarity check
		## Slice around peaks, check fuzzy similarity of densities
		clarity_count = 0
		major_slice = histo_list[major_estimate_bin-2:major_estimate_bin+2]
		minor_slice = histo_list[minor_estimate_bin-2:minor_estimate_bin+2]
		for density in major_slice:
			if np.isclose(major_estimate_sparsity, density):
				clarity_count+=1
		for density in minor_slice:
			if np.isclose(minor_estimate_sparsity, density):
				clarity_count+=1
		if clarity_count > 4:
			self.peak_ambiguous = True

		##
		## Threshold setting
		## Peaks are unclear? Alter threshold accordingly
		major_drop_evaluation = self.drop_calc(self.input_distribution, major_estimate, major_index)
		minor_drop_evaluation = self.drop_calc(self.input_distribution, minor_estimate, minor_index)
		evaluation_total = major_drop_evaluation+minor_drop_evaluation
		if evaluation_total == 2: peak_threshold = 0.35
		if evaluation_total == 1: peak_threshold = 0.45
		if evaluation_total == 0: peak_threshold = 0.55

		##
		## if histogram[estimate's bin] == value found from sparsity, ok
		## increment index from 0 indexed into real CCG val
		## create objects of rounded -- precision not requried to 100s of 0.0
		histogram_derived_major_estimate = np.around(histo_list[major_estimate_bin], 15)
		sparsity_derived_major_estimate = np.around(major_estimate_sparsity, 15)
		histogram_derived_minor_estimate = np.around(histo_list[minor_estimate_bin], 15)
		sparsity_derived_minor_estimate = np.around(minor_estimate_sparsity, 15)

		if self.target_peak_count == 1:
			if histogram_derived_major_estimate == sparsity_derived_major_estimate:
				estimated_attributes['PrimaryPeak'] = major_index+1
				estimated_attributes['SecondaryPeak'] = major_index+1
				estimated_attributes['PeakDistance'] = 0
				estimated_attributes['PeakThreshold'] = peak_threshold
		if self.target_peak_count == 2:
			if (histogram_derived_major_estimate == sparsity_derived_major_estimate) and (histogram_derived_minor_estimate == sparsity_derived_minor_estimate):
				estimated_attributes['PrimaryPeak'] = major_index+1
				estimated_attributes['SecondaryPeak'] = minor_index+1
				estimated_attributes['PeakDistance'] = peak_distance-1
				estimated_attributes['PeakThreshold'] = peak_threshold

		return estimated_attributes

	def differential_peaks(self, attribute_dict, inherited_parameters, threshold_bias):

		##
		## Get relevance distribution info from previous density pass
		## If threshold_bias is true, we're in a re-call, reduce threshold
		## If threshold happens to go < 0, set to 0 (why would it ever get that far?)

		peak_distance = attribute_dict['PeakDistance']
		if not peak_distance:
			peak_distance = 1
		peak_threshold = attribute_dict['PeakThreshold']
		if threshold_bias:
			peak_threshold -= 0.35
			peak_threshold = min(peak_threshold, 0.25)

		##
		## Graph parameters
		linspace_dim = inherited_parameters[0]
		graph_title = inherited_parameters[1]
		axes = inherited_parameters[2]
		filename = inherited_parameters[3]

		##
		## Create plane for calculation/plotting
		x = np.linspace(linspace_dim[0],linspace_dim[1],linspace_dim[2])
		y = np.asarray(self.input_distribution)
		peak_indexes = peakutils.indexes(y, thres=peak_threshold, min_dist=peak_distance)

		if not len(peak_indexes) == self.target_peak_count:
			peak_indexes = peak_indexes[:1]
			self.peak_overpopulated = True
		fixed_indexes = peak_indexes+1

		##
		## Plot graph
		pyplot.figure(figsize=(10,6))
		pyplot.title(graph_title)
		pyplot.xlabel(axes[0])
		pyplot.ylabel(axes[1])
		formatted_label = '{}{}'.format('Allele(s): ', str(fixed_indexes))
		pyplot.plot(x,y, label=formatted_label, markevery=peak_indexes)
		pplot(x,y,peak_indexes)
		pyplot.legend(loc='upper right', shadow=True)
		pyplot.savefig(os.path.join(self.prediction_path, filename), format='png')

		if not len(peak_indexes) == self.target_peak_count:
			raise ValueError
		else:
			if self.target_peak_count == 1:
				attribute_dict['PrimaryPeak'] = fixed_indexes[0]
				attribute_dict['SecondaryPeak'] = fixed_indexes[0]
			if self.target_peak_count == 2:
				attribute_dict['PrimaryPeak'] = fixed_indexes[0]
				attribute_dict['SecondaryPeak'] = fixed_indexes[1]

		return attribute_dict

	def set_distribution(self, new_distribution):

		##
		## In the case where we needed to make a consensus sequence
		## the object's input distribution is changed to the new distribution here
		self.input_distribution = new_distribution

	def get_warnings(self):

		##
		## Method for returning warnings raised in this instance of the object
		return {'Expansion_skew': self.expansion_skew,
				'Peak_ambiguous':self.peak_ambiguous,
				'Density_ambiguous':self.density_ambiguity,
				'PeakOverPopulation':self.peak_overpopulated}

class MosaicismInvestigator:

	def __init__(self, genotype, distribution):

		"""
		Class with functions for investigation somatic mosaicism
		"""

		self.genotype = genotype
		self.distribution = distribution

	def chunks(self, n):

		"""
		Yield successive n-sized chunks from input distribution.
		"""

		for i in xrange(0, len(self.distribution), n):
			yield self.distribution[i:i + n]

	@staticmethod
	def arrange_chunks(ccg_slices):

		"""
		Take distribution lists that have been split into 20 CCG
		Organise into a pandas dataframe for later use
		"""

		arranged_rows = []
		for ccg_value in ccg_slices:
			column = []
			for i in range(0, len(ccg_value)):
				column.append(ccg_value[i])
			arranged_rows.append(column)

		df = pd.DataFrame({'CCG1': arranged_rows[0], 'CCG2': arranged_rows[1], 'CCG3': arranged_rows[2],
						   'CCG4': arranged_rows[3], 'CCG5': arranged_rows[4], 'CCG6': arranged_rows[5],
						   'CCG7': arranged_rows[6], 'CCG8': arranged_rows[7], 'CCG9': arranged_rows[8],
						   'CCG10': arranged_rows[9], 'CCG11': arranged_rows[10], 'CCG12': arranged_rows[11],
						   'CCG13': arranged_rows[12], 'CCG14': arranged_rows[13], 'CCG15': arranged_rows[14],
						   'CCG16': arranged_rows[15], 'CCG17': arranged_rows[16], 'CCG18': arranged_rows[17],
						   'CCG19': arranged_rows[18], 'CCG20': arranged_rows[19]})

		return df

	@staticmethod
	def get_nvals(df, input_allele):

		"""
		Take specific CCG sub-list from dataframe
		based on genotype taken from GenotypePrediction
		Extract n-1/n/n+1
		"""
		allele_nvals = {}
		cag_value = input_allele[0]
		ccgframe = df['CCG'+str(input_allele[1])]

		try:
			nminus = str(ccgframe[int(cag_value)-2])
			nvalue = str(ccgframe[int(cag_value)-1])
			nplus = str(ccgframe[int(cag_value)])
		except KeyError:
			log.info('{}{}{}{}'.format(clr.red,'shd__ ',clr.end,'N-Value scraping Out of Bounds.'))

		allele_nvals['NMinusOne'] = nminus
		allele_nvals['NValue'] = nvalue
		allele_nvals['NPlusOne'] = nplus

		return allele_nvals

	@staticmethod
	def calculate_mosaicism(allele_values):

		"""
		Add additional calculations here..
		"""

		nmo = allele_values['NMinusOne']
		n = allele_values['NValue']
		npo = allele_values['NPlusOne']
		nmo_over_n = 0
		npo_over_n = 0

		try:
			nmo_over_n = int(nmo) / int(n)
			npo_over_n = int(npo) / int(n)
		except ZeroDivisionError:
			log.info('{}{}{}{}'.format(clr.red,'shd__ ',clr.end,' Divide by 0 attempted in Mosaicism Calculation.'))

		calculations = {'NMinusOne':nmo,'NValue':n,'NPlusOne':npo,'NMinusOne-Over-N': nmo_over_n, 'NPlusOne-Over-N': npo_over_n}
		return calculations

	@staticmethod
	def distribution_padder(ccg_dataframe, genotype):

		"""
		Ensure that all distribution's N will be in the same position
		in a null-padded larger distribution
		for manual comparison of larger somatic mosaicism trends
		"""

		unpadded_distribution = list(ccg_dataframe['CCG'+str(genotype[1])])
		n_value = genotype[0]

		anchor = 203
		anchor_to_left = anchor - n_value
		anchor_to_right	= anchor_to_left + 200
		left_buffer = ['-'] * anchor_to_left
		right_buffer = ['-'] * (403-anchor_to_right)
		padded_distribution = left_buffer + unpadded_distribution + right_buffer

		return padded_distribution