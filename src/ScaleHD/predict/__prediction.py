from __future__ import division

#/usr/bin/python
__version__ = 0.01
__author__ = 'alastair.maxwell@glasgow.ac.uk'

##
## Generic imports
import os
import csv
import peakutils
import numpy as np
import logging as log
from collections import Counter
from sklearn import svm
from sklearn.multiclass import OutputCodeClassifier
from sklearn import preprocessing
import matplotlib
matplotlib.use('Agg') #servers/clients without x-11
from peakutils.plot import plot as pplot
import matplotlib.pyplot as plt
import pandas as pd

##
## Backend Junk
from ..__backend import Colour as clr
from ..__backend import DataLoader

class GenotypePrediction:
	def __init__(self, data_pair, prediction_path, training_data, instance_params):
		"""
		Prediction stage of the pipeline -- use of SVM, density estimation and first order differentials.
		Automates calling of sample's genotype based on data information dervied from aligned read counts.
		Utilises forward reads for CAG information, reverse reads for CCG information. Combine for genotype.

		General workflow overview:
		--Take reverse reads, aggregate every CAG for each CCG bin
		--Use unlabelled sample into CCG zygosity SVM for het/hom prediction
		--Data cleaning/normalisation/etc
		--Two Pass algorithm to determine genotype
		--1) Density Estimation on distribution to gauge where the peaks may be/peak distances
		--2) Peak Detection via first order differentials, taking into account density results for tailoring
		--Repeat process for relevant CAG distribution(s) taking CCG het/hom into account
		--Return genotype

		:param data_pair: Files to be used for scraping
		:param prediction_path: Output path to save all resultant files from this process
		:param training_data: Data to be used in building CCG SVM model
		:param instance_params: redundant parameters--unused in this build
		"""

		##
		## Paths and files to be used in this instance
		self.data_pair = data_pair
		self.prediction_path = prediction_path
		self.training_data = training_data
		self.instance_params = instance_params

		##
		## Build a classifier and class label hash-encoder for CCG SVM
		self.classifier, self.encoder = self.build_zygosity_model()

		"""
		Information/Error flags that exist within this class::
		--CCGZygDisconnect :: Forward and Reverse CCG SVM predictions differed
		--CCGExpansionSkew  :: CCG Peak - N-1 of Major Peak is > in (val) than N of Minor Peak
		--CCGPeakAmbiguous :: Multiple low densities surrounding a suspected peak density
		--CCGDensityAmbiguous :: Many low densities spread across the KDE histogram for a distribution
		--CCGRecallWarning :: CCG Detection had to be re-called for one reason or another
		--CCGPeakOOB :: Too many peaks were called (but threshold was lowered for a reason), truncated results are returned but confidence is low
		--CAGRecallWarning :: CAG Detection had to be re-called for one reason or another
		--CAGPeakOOB :: Too many peaks were called (but threshold was lowered for a reason), truncated results are returned but confidence is low
		--FPSPDisconnect :: First and Second passes differed in results.. not a big deal but w/e

		And data attributes that are used as output::
		PrimaryAllele = [CAGa, CCGb]
		SecondaryAllele = [CAGc, CCGd]
		PrimaryMosaicism = [<values>]
		SecondaryMosaicism = [<values>]
		"""
		self.prediction_confidence = 0
		self.cag_intermediate = [0,0]
		self.genotype_flags = {'PrimaryAllele':[0,0],
							   'SecondaryAllele':[0,0],
							   'PrimaryMosaicism':[],
							   'SecondaryMosaicism':[],
							   'ThresholdUsed':0,
							   'RecallCount':0,
							   'AlignmentPadding':False,
							   'SVMPossibleFailure':False,
							   'PotentialHomozygousHaplotype':False,
							   'PHHInterpDistance':0.0,
							   'NeighbouringPeaks':False,
							   'DiminishedPeak':False,
							   'DiminishedUncertainty':False,
							   'CCGZygDisconnect':False,
							   'CCGExpansionSkew':False,
							   'CCGPeakAmbiguous':False,
							   'CCGDensityAmbiguous':False,
							   'CCGRecallWarning':False,
							   'CCGPeakOOB':False,
							   'CAGExpansionSkew':False,
							   'CAGPeakAmbiguous':False,
							   'CAGDensityAmbiguous':False,
							   'CAGRecallWarning':False,
							   'CAGPeakOOB':False,
							   'CAGBackwardsSlippage':False,
							   'CAGForwardSlippage':False,
							   'FPSPDisconnect':False}

		##
		## Unlabelled distributions to utilise for SVM prediction
		## Padded distro = None, in case where SAM aligned to (0-100,0-20 NOT 1-200,1-20), use this
		self.forward_distribution = self.scrape_distro(self.data_pair[0])
		self.reverse_distribution = self.scrape_distro(self.data_pair[1])
		self.forward_distr_padded = None

		"""
		!! Stage one !!
		Determine Zygosity of CCG from input distribution
		-- Aggregate CCG reads from 200x20 to 1x20
		-- Feed into SVM
		-- Compare results between forward and reverse (reverse takes priority)
		"""

		self.forwardccg_aggregate = self.distribution_collapse(self.forward_distribution, st=True)
		self.reverseccg_aggregate = self.distribution_collapse(self.reverse_distribution)
		self.zygosity_state = self.predict_zygstate()

		"""
		!! Stage two !!
		Determine CCG Peak(s)/Genotype(s) via 2-Pass Algorithm
		Run first attempt with no clauses; if pass, continue to next stage
		However, if something fails, a loop will trigger until the function passes
		"""
		ccg_failstate, ccg_genotype = self.determine_ccg_genotype()
		while ccg_failstate:
			self.genotype_flags['CCGRecallWarning'] = True
			ccg_failstate, ccg_genotype = self.determine_ccg_genotype(threshold_bias=True)

		self.genotype_flags['PrimaryAllele'][1] = ccg_genotype[0]
		self.genotype_flags['SecondaryAllele'][1] = ccg_genotype[1]

		"""
		!! Stage three !!
		Now we have identified the CCG peaks (successfully), we can investigate the appropriate
		CAG distributions for these CCG distribution(s). The same generic functions will be called
		for CAG determination, and results from CCG and CAG are combined to produce a genotype for this sample
		"""
		cag_failstate, cag_genotype = self.determine_cag_genotype()
		while cag_failstate:
			self.genotype_flags['CAGRecallWarning'] = True
			cag_failstate, cag_genotype = self.determine_cag_genotype(threshold_bias=True)

		self.genotype_flags['PrimaryAllele'][0] = cag_genotype[0]
		self.genotype_flags['SecondaryAllele'][0] = cag_genotype[1]

		##
		## Ensure that normal and expanded alleles are located in the correct position
		## Based on CAG ordering value, NOT CCG ordering value
		if int(self.genotype_flags['PrimaryAllele'][0]) > int(self.genotype_flags['SecondaryAllele'][0]):
			intermediate = self.genotype_flags['PrimaryAllele']
			self.genotype_flags['PrimaryAllele'] = self.genotype_flags['SecondaryAllele']
			self.genotype_flags['SecondaryAllele'] = intermediate

		"""
		!! Stage four !!
		Simple Somatic Mosaicism calculations are done here
		"""
		self.genotype_flags['PrimaryMosaicism'] = self.somatic_calculations(self.genotype_flags['PrimaryAllele'])
		self.genotype_flags['SecondaryMosaicism'] = self.somatic_calculations(self.genotype_flags['SecondaryAllele'])

		"""
		!! Stage five !!
		Determine confidence in this genotype prediciton, and return for report
		"""
		self.confidence_calculation()
		self.gtype_report = self.generate_report()

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
		ovo_svc = OutputCodeClassifier(svc_object, code_size=2, random_state=0).fit(X,Y)
		encoder = traindat_model.ENCDR

		##
		## Return the fitted OvO(SVM) and Encoder
		return ovo_svc, encoder

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
			next(source) #skip header
			for row in source:
				placeholder_array.append(int(row[2]))
			dfi.close()
		unlabelled_distro = np.array(placeholder_array)
		return unlabelled_distro

	def distribution_collapse(self, distribution_array, st=False):
		"""
		Function to take a full 200x20 array (struc: CAG1-200,CCG1 -- CAG1-200CCG2 -- etc CCG20)
		and aggregate all CAG values for each CCG
		:param distribution_array: input dist (should be (1-200,1-20))
		:param st: flag for if we're in forward stage; i.e. input aligned to wrong ref -> assign padded to fwrd
		:return: 1x20D np(array)
		"""

		##
		## Object for CCG split
		ccg_arrays = None

		##
		## Hopefully the user has aligned to the right reference dimensions
		## Check.. if not, hopefully we can pad (and raise flag.. since it is not ideal)
		try:
			ccg_arrays = np.split(distribution_array, 20)
		##
		## User aligned to the wrong reference..
		except ValueError:
			self.genotype_flags['AlignmentPadding'] = True

			##
			## If the distro is this size, they aligned to (0-100,0-20)...
			## Split by 21, append with 99x1's to end of each CCG
			## Trim first entry in new list (CCG0.. lol who even studies that)
			## If we're on a forward distro collapse, assign padded distro (for mosaicism later)
			if len(distribution_array) == 2121:
				altref_split = np.split(distribution_array, 21)
				padded_split = []
				for ccg in altref_split:
					current_pad = np.append(ccg[1:], np.ones(100))
					padded_split.append(current_pad)
				ccg_arrays = padded_split[1:]

				if st:
					self.forward_distr_padded = np.asarray([item for sublist in ccg_arrays for item in sublist])

		##
		## Aggregate each CCG
		ccg_counter = 1
		collapsed_array = []
		for ccg_array in ccg_arrays:
			collapsed_array.append(np.sum(ccg_array))
			ccg_counter+=1
		return np.asarray(collapsed_array)

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
		forward_reshape = preprocessing.normalize(np.float64(self.forwardccg_aggregate.reshape(1,-1)))
		reverse_reshape = preprocessing.normalize(np.float64(self.reverseccg_aggregate.reshape(1,-1)))

		##
		## Predict the zygstate of these reshapen, noramlised 20D CCG arrays using SVM object earlier
		## Results from self.classifier are #encoded; so convert with our self.encoder.inverse_transform
		forward_zygstate = str(self.encoder.inverse_transform(self.classifier.predict(forward_reshape)))
		reverse_zygstate = str(self.encoder.inverse_transform(self.classifier.predict(reverse_reshape)))

		##
		## We only particularly care about the reverse zygosity (CCG reads are higher quality in reverse data)
		## However, for a QoL metric, compare fw/rv results. If match, good! If not, who cares!
		if not forward_zygstate == reverse_zygstate:
			self.genotype_flags['CCGZygDisconnect'] = True
		else:
			self.genotype_flags['CCGZyg_disconnect'] = False
		return reverse_zygstate[2:-2]

	def update_flags(self, target_updates):
		"""
		Function that will take a list of flags that were raised from the 2-Pass algorithm
		and update this pipeline's instance of self.genotype_flags accordingly.
		This allows us to keep a current state-of-play of this sample's prediction.
		:param target_updates: List of flags from the 2-Pass algorithm
		:return: None
		"""

		for update_key, update_value in target_updates.iteritems():
			for initial_key, initial_value in self.genotype_flags.iteritems():
				if initial_key == update_key:
					self.genotype_flags[initial_key] = update_value

	def determine_ccg_genotype(self, fail_state=False, threshold_bias=False):
		"""
		Function to determine the genotype of this sample's CCG alleles
		Ideally this function will be called one time, but where exceptions occur
		it may be re-called with a lower quality threshold -- inform user when this occurs
		:param fail_state: optional flag for re-calling when a previous call failed
		:param threshold_bias: optional flag for lowering FOD threshold when a previous called failed
		:return: failure state, CCG genotype data ([None,X],[None,Y])
		"""

		##
		## We're going to be paranoid and bootstrap a check on the zygstate call
		def paranoia(input_dist):
			## Top1/2 for CCG (top3 not applicable)
			major_estimate = max(input_dist)
			minor_estimate = max(n for n in input_dist if n!=major_estimate)
			## percentage drop from maj to min
			literal_drop = (abs(major_estimate-minor_estimate) / major_estimate) * 100

			## if minor is 60% of major, or less, maybe we are het and the SVM was wrong
			if literal_drop	<= 60.00:
				self.genotype_flags['SVMPossibleFailure'] = True
				self.zygosity_state = 'HETERO'
				return 2
			else:
				self.zygosity_state = 'HOMO'
				return 1
		peak_target = 0
		if self.zygosity_state == 'HOMO': peak_target = 1
		if self.zygosity_state == 'HETERO': peak_target = 2
		if peak_target == 1: peak_target = paranoia(self.reverseccg_aggregate)

		##
		## Create object for 2-Pass algorithm to use with CCG
		graph_parameters = [20, 'CCGDensityEstimation.png', 'CCG Density Distribution', ['Read Count', 'Bin Density']]
		ccg_inspector = SequenceTwoPass(prediction_path=self.prediction_path,
										input_distribution=self.reverseccg_aggregate,
										peak_target=peak_target,
										graph_parameters=graph_parameters,
										zygosity_state = self.zygosity_state,
										contig_stage='CCG')

		"""
		!! Sub-Stage one !!
		Now that we've made an object with the settings for this instance..
		Density estimation of the CCG distribution..
		Get warnings encountered by this instance of SequenceTwoPass
		Update equivalent warning flags within GenotypePrediction
		"""
		first_pass = ccg_inspector.density_estimation(plot_flag=False)
		density_warnings = ccg_inspector.get_warnings()
		self.update_flags(density_warnings)

		"""
		!! Sub-Stage two !!
		Now we have our estimates from the KDE sub-stage, we can use these findings
		in our FOD peak identification for more specific peak calling and thus, genotyping
		"""
		fod_param = [[0,20,21],'CCG Peaks',['CCG Value', 'Read Count'], 'CCGPeakDetection.png']
		fod_failstate, second_pass = ccg_inspector.differential_peaks(first_pass, fod_param, threshold_bias)
		while fod_failstate:
			fod_failstate, second_pass = ccg_inspector.differential_peaks(first_pass, fod_param, threshold_bias, fod_recall=True)
		differential_warnings = ccg_inspector.get_warnings()
		self.update_flags(differential_warnings)

		##
		## Check if First Pass Estimates == Second Pass Results
		## If there is a mismatch, genotype calling has failed and thus re-call will be required
		first_pass_estimate = [first_pass['PrimaryPeak'],first_pass['SecondaryPeak']]
		second_pass_estimate = [second_pass['PrimaryPeak'], second_pass['SecondaryPeak']]

		if not first_pass_estimate == second_pass_estimate or len(second_pass_estimate)>len(first_pass_estimate):
			self.genotype_flags['FPSPDisconnect'] = True
			fail_state = True

		##
		## Return whether this process passed or not, and the genotype
		return fail_state, second_pass_estimate

	@staticmethod
	def split_cag_target(input_distribution, ccg_target):
		"""
		Function to gather the relevant CAG distribution for the specified CCG value
		We gather this information from the forward distribution of this sample pair as CCG reads are
		of higher quality in the forward sequencing direction.
		We split the entire fw_dist into contigs/bins for each CCG (4000 -> 200*20)
		:param input_distribution: input forward distribution (4000d)
		:param ccg_target: target value we want to select the 200 values for
		:return: the sliced CAG distribution for our specified CCG value
		"""

		cag_split = [input_distribution[i:i+200] for i in xrange(0, len(input_distribution), 200)]
		distribution_dict = {}
		for i in range(0, len(cag_split)):
			distribution_dict['CCG'+str(i+1)] = cag_split[i]

		current_target_distribution = distribution_dict['CCG' + str(ccg_target)]
		return current_target_distribution

	def determine_cag_genotype(self, fail_state=False, threshold_bias=False):
		"""
		Function to determine the genotype of this sample's CAG alleles
		Ideally this function will be called one time, but where exceptions occur,
		it may be re-caled with a lower quality-threshold -- inform user when this occurs

		If CCG was homozygous (i.e. one CCG distro) -- there will be 2 CAG peaks to investigate
		If CCG was heterozygous (i.e. two CCG distro) -- there will be 1 CAG peak in each CCG to investigate

		:param fail_state: optional flag for re-calling when a previous call failed
		:param threshold_bias: optional flag for lowering FOD threshold when a previous call failed
		:return: failure state, CAG genotype data ([X,None],[Y,None])
		"""

		##
		## Check for padding status
		## If user aligned to wrong reference, we need to use self.forward_distr_padded instead of
		## self.forward_distribution; otherwise dimensionality is incorrect...

		if not self.genotype_flags['AlignmentPadding']: forward_utilisation = self.forward_distribution
		else: forward_utilisation = self.forward_distr_padded

		##
		## Set up distributions we require to investigate
		## If Homozygous, we have one CCG distribution that will contain 2 CAG peaks to investigate
		## If Heterozygous, we have two CCG distributions, each with 1 CAG peak to investigate
		peak_target = 0
		target_distribution = {}
		if self.zygosity_state == 'HOMO':
			peak_target = 2
			cag_target = self.split_cag_target(forward_utilisation, self.genotype_flags['PrimaryAllele'][1])
			target_distribution[self.genotype_flags['PrimaryAllele'][1]] = cag_target
		if self.zygosity_state == 'HETERO':
			peak_target = 1
			cag_target_major = self.split_cag_target(forward_utilisation, self.genotype_flags['PrimaryAllele'][1])
			cag_target_minor = self.split_cag_target(forward_utilisation, self.genotype_flags['SecondaryAllele'][1])
			target_distribution[self.genotype_flags['PrimaryAllele'][1]] = cag_target_major
			target_distribution[self.genotype_flags['SecondaryAllele'][1]] = cag_target_minor

		##
		## Now iterate over our scraped distributions with our 2 pass algorithm
		for cag_key, distro_value in target_distribution.iteritems():

			##
			## Generate KDE graph parameters
			## Generate CAG inspector Object for 2Pass-Algorithm
			graph_parameters = [20, '{}{}{}'.format('CAG',str(cag_key),'DensityEstimation.png'), 'CAG Density Distribution', ['Read Count', 'Bin Density']]
			cag_inspector = SequenceTwoPass(prediction_path=self.prediction_path,
											input_distribution=distro_value,
											peak_target=peak_target,
											graph_parameters=graph_parameters,
											zygosity_state=self.zygosity_state,
											contig_stage='CAG')

			"""
			!! Sub-stage one !!
			Now that we've made an object with the settings for this instance..
			Density estimation of the CCG distribution..
			Get warnings encountered by this instance of SequenceTwoPass
			Update equivalent warning flags within GenotypePrediction
			"""
			first_pass = cag_inspector.density_estimation(plot_flag=False)
			density_warnings = cag_inspector.get_warnings()
			self.update_flags(density_warnings)

			"""
			!! Sub-stage two !!
			Now we have our estimates from the KDE sub-stage, we can use these findings
			in our FOD peak identification for more specific peak calling and thus, genotyping
			"""
			fod_param = [[0,200,201],'{}{})'.format('CAG Peaks (CCG',str(cag_key)),['CAG Value', 'Read Count'], '{}{}{}'.format('CCG',str(cag_key),'-CAGPeakDetection.png')]
			fod_failstate, second_pass = cag_inspector.differential_peaks(first_pass, fod_param, threshold_bias)
			while fod_failstate:
				fod_failstate, second_pass = cag_inspector.differential_peaks(first_pass, fod_param, threshold_bias, fod_recall=True)
			differential_warnings = cag_inspector.get_warnings()
			self.update_flags(differential_warnings)

			##
			## Concatenate results into a sample-wide genotype format
			first_pass_estimate = [first_pass['PrimaryPeak'], first_pass['SecondaryPeak']]
			second_pass_estimate = [second_pass['PrimaryPeak'], second_pass['SecondaryPeak']]

			if not first_pass_estimate == second_pass_estimate or len(second_pass_estimate)>len(first_pass_estimate):
				self.genotype_flags['FPSPDisconnect'] = True
				fail_state = True

			##
			## Ensure the correct CAG is assigned to the appropriate CCG
			if self.zygosity_state == 'HOMO':
				if cag_key == self.genotype_flags['PrimaryAllele'][1]:
					self.cag_intermediate[0] = second_pass_estimate[0]
					self.cag_intermediate[1] = second_pass_estimate[1]
			if self.zygosity_state == 'HETERO':
				if cag_key == self.genotype_flags['PrimaryAllele'][1]:
					self.cag_intermediate[0] = second_pass_estimate[0]
				if cag_key == self.genotype_flags['SecondaryAllele'][1]:
					self.cag_intermediate[1] = second_pass_estimate[0]

		##
		## Generate object and return
		cag_genotype = [self.cag_intermediate[0],self.cag_intermediate[1]]
		return fail_state, cag_genotype

	def somatic_calculations(self, genotype):
		"""
		Function for basic somatic mosaicism calculations; featureset will be expanded upon later
		For now; N-1 / N, N+1 / N calculations are executed on arranged contigs where the N value is known
		(from genotype prediction -- assumed to be correct)
		In addition, the read count distribution for the forward and reverse reads in a sample pair are both
		aligned so that their N value is in the same position; lets end-user investigate manual distribution
		data quality etc.
		:param genotype: value of predicted genotype from the SVM/2PA stages of GenotypePrediction()
		:return: mosaicism_values; results of simple sommos calculations :))))))
		"""

		##
		## Create mosaicism investigator object to begin calculation prep
		## Takes raw 200x20 dist and slices into 20 discrete 200d arrays
		## Orders into a dataframe with CCG<val> labels
		## Scrapes appropriate values for SomMos calculations
		if not self.genotype_flags['AlignmentPadding']:
			mosaicism_object = MosaicismInvestigator(genotype, self.forward_distribution)
		else:
			mosaicism_object = MosaicismInvestigator(genotype, self.forward_distr_padded)
		ccg_slices = mosaicism_object.chunks(200)
		ccg_ordered = mosaicism_object.arrange_chunks(ccg_slices)
		allele_values = mosaicism_object.get_nvals(ccg_ordered, genotype)

		##
		## With these values, we can calculate and return
		allele_calcs = mosaicism_object.calculate_mosaicism(allele_values)

		##
		## Generate a padded distribution (aligned to N=GTYPE)
		padded_distro = mosaicism_object.distribution_padder(ccg_ordered, genotype)

		##
		## Combine calculation dictionary and distribution into object, return
		## TODO rework padded distributions?
		mosaicism_values = [allele_calcs, padded_distro]
		return mosaicism_values

	def confidence_calculation(self):
		"""
		Function that will calculate a score for confidence in the genotype
		prediction for the current sample being processed...
		Check which flags have been raised throughout the genotyping process
		Weight accordingly, return
		:return: None
		"""
		current_confidence = 100

		"""
		>> Highest severity <<
		If these flags are raised/over a threshold then there
		has probably been a strong effect on the precision and accuracy of
		the genotype prediction.. user should probably manually check results
		"""

		##
		## Threshold utilisation during FOD peak calling
		if self.genotype_flags['ThresholdUsed'] != 0.5:
			if self.genotype_flags['ThresholdUsed'] < 0.5: current_confidence -= 15
			elif self.genotype_flags['ThresholdUsed'] < 0.3: current_confidence -= 20
			else: current_confidence -= 30
		else: current_confidence += 20

		##
		## Recall count/CAG-CCG specific recall warnings
		if self.genotype_flags['RecallCount'] != 0:
			if self.genotype_flags['RecallCount'] < 3: current_confidence -= 10
			elif self.genotype_flags['RecallCount'] <= 4: current_confidence -= 15
			else: current_confidence -= 25
			if self.genotype_flags['CCGRecallWarning']: current_confidence -= 15
			if self.genotype_flags['CAGRecallWarning']: current_confidence -= 10
		else: current_confidence += 20

		##
		## SVM Possible failure
		if self.genotype_flags['SVMPossibleFailure']: current_confidence -= 30

		##
		## PeakOOB: >2 peaks returned per allele (i.e. results were sliced)
		for peakoob in [self.genotype_flags['CAGPeakOOB'],self.genotype_flags['CCGPeakOOB']]:
			if peakoob: current_confidence -= 15

		##
		## Sample wasn't aligned to CAG1-200/CCG1-20.. padded but raises questions...
		if self.genotype_flags['AlignmentPadding']: current_confidence -= 25

		"""
		>> Medium severity <<
		If these flags have been raised/over a threshold, then there
		has probably been a mild effect on the precision/accuracy, but not
		necessarily.. Different weights of flags alter the severity outcome..
		"""

		##
		## Homozygous Haplotype detection?
		if self.genotype_flags['PotentialHomozygousHaplotype']:
			current_confidence -= 5
			if self.genotype_flags['PHHInterpDistance'] > 1.0:
				current_confidence -= 5

		##
		## Peaks are neighbouring? (e.g. 16/17)
		if self.genotype_flags['NeighbouringPeaks']: current_confidence -= 15

		##
		## Diminished peak detection
		if self.genotype_flags['DiminishedPeak']:
			current_confidence -= 2.5
			if self.genotype_flags['DiminishedUncertainty']: current_confidence -= 7.5

		##
		## Peak / Density ambiguity (only matters if re-call occurred)
		if self.genotype_flags['RecallCount'] != 0:
			if self.genotype_flags['CAGPeakAmbiguous']: current_confidence -= 5
			if self.genotype_flags['CCGPeakAmbiguous']: current_confidence -= 10
			if self.genotype_flags['CAGDensityAmbiguous']: current_confidence -= 10
			if self.genotype_flags['CCGDensityAmbiguous']: current_confidence -= 15

		"""
		>> Lowest severity <<
		If these flags are raised, it is for information only
		It is highly doubtful that anything encountered here would negatively
		alter the accuracy/precision of a genotype prediction
		"""

		##
		## Slippage
		if self.genotype_flags['CAGBackwardsSlippage']: current_confidence -= 5
		if self.genotype_flags['CAGForwardSlippage']: current_confidence -= 2

		"""
		>> Mosaicism Investigation <<
		Based on how much somatic mosaicism we see around the peaks..
		That could maybe have altered the genotype prediction away from the
		true value and/or altered accuracy/precision
		"""

		##
		## Mosaicism!
		if self.genotype_flags['PrimaryMosaicism'][0]['NMinusOne-Over-N'] > 0.30: current_confidence -= 2.5
		if self.genotype_flags['PrimaryMosaicism'][0]['NPlusOne-Over-N'] > 0.25: current_confidence -= 4.25
		if self.genotype_flags['SecondaryMosaicism'][0]['NMinusOne-Over-N'] > 0.65: current_confidence -= 7.5
		if self.genotype_flags['SecondaryMosaicism'][0]['NPlusOne-Over-N'] > 0.70: current_confidence -= 10

		##
		## With all flags processed; assign confidence score for this instance
		## Limit output to 100%
		self.prediction_confidence = sorted([0, current_confidence, 100])[1]

	def generate_report(self):
		"""
		Function which will, eventually, calculate the confidence score of this genotype prediction
		by taking into account flags raised, and meta-data about the current sample distribution etc
		:return: For now, a list of report flags. eventually, probably a dictionary with more info within
		"""

		##
		## Hideous string based report for individual samples
		## Will get changed at a later date
		sample_name = self.prediction_path.split('/')[-2]
		sample_report_name = os.path.join(self.prediction_path, sample_name+'QuickReport.txt')
		sample_report = '{}: {}\n{}: {}\n' \
						'{}: {}\n{}: {}%\n' \
						'{}: {}\n{}: {}\n' \
						'{}: {}\n{}: {}\n' \
						'{}: {}\n{}: {}\n' \
						'{}: {}\n{}: {}\n' \
						'{}: {}\n{}: {}\n' \
						'{}: {}\n{}: {}\n' \
						'{}: {}\n{}: {}\n' \
						'{}: {}\n{}: {}\n' \
						'{}: {}\n{}: {}\n' \
						'{}: {}\n{}: {}\n' \
						'{}: {}\n{}: {}\n' \
						'{}: {}\n{}: {}\n' \
						'{}: {}'.format('File Name', sample_name,
										'Primary Allele', self.genotype_flags['PrimaryAllele'],
										'Secondary Allele', self.genotype_flags['SecondaryAllele'],
										'Prediction Confidence', self.prediction_confidence,
										'Threshold Used', self.genotype_flags['ThresholdUsed'],
										'Recall Count', self.genotype_flags['RecallCount'],
										'Alignment Padding', self.genotype_flags['AlignmentPadding'],
										'\nCCG Flags', '',
										'CCG Zygosity Disconnect', self.genotype_flags['CCGZygDisconnect'],
										'CCG Expansion Skew', self.genotype_flags['CCGExpansionSkew'],
										'CCG Peak Ambiguity', self.genotype_flags['CCGPeakAmbiguous'],
										'CCG Density Ambiguity', self.genotype_flags['CCGDensityAmbiguous'],
										'CCG Recall Warning', self.genotype_flags['CCGRecallWarning'],
										'CCG Peak OOB', self.genotype_flags['CCGPeakOOB'],
										'\nCAG Flags', '',
										'CAG Expansion Skew', self.genotype_flags['CAGExpansionSkew'],
										'CAG Peak Ambiguity', self.genotype_flags['CAGPeakAmbiguous'],
										'CAG Density Ambiguity', self.genotype_flags['CAGDensityAmbiguous'],
										'CAG Recall Warning', self.genotype_flags['CAGRecallWarning'],
										'CAG Peak OOB', self.genotype_flags['CAGPeakOOB'],
										'CAG Backwards Slippage', self.genotype_flags['CAGBackwardsSlippage'],
										'CAG Forwards Slippage', self.genotype_flags['CAGForwardSlippage'],
										'\nOther Flags', '',
										'FPSP Disconnect', self.genotype_flags['FPSPDisconnect'],
										'Homozygous Haplotype', self.genotype_flags['PotentialHomozygousHaplotype'],
										'Haplotype Interp Distance', self.genotype_flags['PHHInterpDistance'],
										'Neighbouring Peaks', self.genotype_flags['NeighbouringPeaks'],
										'Diminished Peak', self.genotype_flags['DiminishedPeak'],
										'Diminished Peak Uncertainty', self.genotype_flags['DiminishedUncertainty'])
		sample_file = open(sample_report_name, 'w')
		sample_file.write(sample_report)
		sample_file.close()

		report = {'PrimaryAllele':self.genotype_flags['PrimaryAllele'],
				  'SecondaryAllele':self.genotype_flags['SecondaryAllele'],
				  'PredictionConfidence':self.prediction_confidence,
				  'PrimaryMosaicism':self.genotype_flags['PrimaryMosaicism'],
				  'SecondaryMosaicism':self.genotype_flags['SecondaryMosaicism'],
				  'ThresholdUsed':self.genotype_flags['ThresholdUsed'],
				  'RecallCount':self.genotype_flags['RecallCount'],
				  'AlignmentPadding':self.genotype_flags['AlignmentPadding'],
				  'SVMPossibleFailure':self.genotype_flags['SVMPossibleFailure'],
				  'PotentialHomozygousHaplotype':self.genotype_flags['PotentialHomozygousHaplotype'],
				  'PHHIntepDistance':self.genotype_flags['PHHInterpDistance'],
				  'NeighbouringPeaks':self.genotype_flags['NeighbouringPeaks'],
				  'DiminishedPeak':self.genotype_flags['DiminishedPeak'],
				  'DiminishedUncertainty':self.genotype_flags['DiminishedUncertainty'],
				  'CCGZygDisconnect':self.genotype_flags['CCGZygDisconnect'],
				  'CCGExpansionSkew':self.genotype_flags['CCGExpansionSkew'],
				  'CCGPeakAmbiguous':self.genotype_flags['CCGPeakAmbiguous'],
				  'CCGDensityAmbiguous':self.genotype_flags['CCGDensityAmbiguous'],
				  'CCGRecallWarning':self.genotype_flags['CCGRecallWarning'],
				  'CCGPeakOOB':self.genotype_flags['CCGPeakOOB'],
				  'CAGExpansionSkew':self.genotype_flags['CAGExpansionSkew'],
				  'CAGPeakAmbiguous':self.genotype_flags['CAGPeakAmbiguous'],
				  'CAGDensityAmbiguous':self.genotype_flags['CAGDensityAmbiguous'],
				  'CAGRecallWArning':self.genotype_flags['CAGRecallWarning'],
				  'CAGPeakOOB':self.genotype_flags['CAGPeakOOB'],
				  'CAGBackwardsSlippage':self.genotype_flags['CAGBackwardsSlippage'],
				  'CAGForwardSlippage':self.genotype_flags['CAGForwardSlippage'],
				  'FPSPDisconnect':self.genotype_flags['FPSPDisconnect']}

		return report

	def get_report(self):
		"""
		Function to just return the report for this class object from the point of calling
		:return: a report. wow
		"""
		return self.gtype_report

class SequenceTwoPass:
	def __init__(self, prediction_path, input_distribution, peak_target, graph_parameters, zygosity_state, contig_stage):
		"""
		Class that will be used as an object for each genotyping stage of the GenotypePrediction pipe
		Each function within this class has it's own doctstring for further explanation
		This class is called into an object for each of CCG/CAG deterministic stages

		:param prediction_path: Output path to save all resultant files from this process
		:param input_distribution: Distribution to put through the two-pass (CAG or CCG..)
		:param peak_target: Number of peaks we expect to see in this current distribution
		:param graph_parameters: Parameters (names of axes etc..) for saving results to graph
		:param zygosity_state: hetero/homozygous (labelling indexing issues require this)
		"""

		##
		## Variables for this instance of this object
		self.prediction_path = prediction_path
		self.input_distribution = input_distribution
		self.peak_target = peak_target
		self.bin_count = graph_parameters[0]
		self.filename = graph_parameters[1]
		self.graph_title = graph_parameters[2]
		self.axes = graph_parameters[3]
		self.zygosity_state = zygosity_state
		self.contig_stage = contig_stage
		self.instance_parameters = {}

		##
		## Potential warnings raised in this instance/useful variables
		self.CCGDensityAmbiguous = False
		self.CCGExpansionSkew = False
		self.CCGPeakAmbiguous = False
		self.CCGPeakOOB = False
		self.PotentialHomozygousHaplotype = False
		self.PHHInterpDistance = 0.0
		self.NeighbouringPeaks = False
		self.DiminishedPeak = False
		self.DiminishedUncertainty = False
		self.CAGDensityAmbiguous= False
		self.CAGExpansionSkew = False
		self.CAGPeakAmbiguous = False
		self.CAGPeakOOB = False
		self.CAGBackwardsSlippage = False
		self.CAGForwardSlippage = False
		self.ThresholdUsed = 0
		self.RecallCount = 0

	def histogram_generator(self, filename, graph_title, axes, plot_flag):
		"""
		Generate histogrm of kernel density estimation for this instance of 2PA
		:param filename: Filename for graph to be saved as..
		:param graph_title: self explanatory
		:param axes: self explanatory
		:param plot_flag: CCG? Plot KDE. CAG? Don't.
		:return: histogram, bins
		"""

		##
		## Generate KDE histogram and plot to graph
		hist, bins = np.histogram(self.input_distribution, bins=self.bin_count, density=True)
		if plot_flag:
			plt.figure(figsize=(10,6))
			bin_width = 0.7 * (bins[1] - bins[0])
			center = (bins[:-1] + bins[1:]) / 2
			plt.title(graph_title)
			plt.xlabel(axes[0])
			plt.ylabel(axes[1])
			plt.bar(center, hist, width=bin_width)
			plt.savefig(os.path.join(self.prediction_path, filename), format='png')
			plt.close()

		##
		## Check the number of densities that exist within our histogram
		## If there are many (>2) values that are very low density (i.e. relevant)
		## then raise the flag for density ambiguity -- there shouldn't be many values so issue with data
		density_frequency = Counter(hist)
		for key, value in density_frequency.iteritems():
			if not key == np.float64(0.0) and value > 2:
				if self.contig_stage == 'CCG': self.CCGDensityAmbiguous = True
				if self.contig_stage == 'CAG': self.CAGDensityAmbiguous = True

		##
		## Return histogram and bins to where this function was called
		return hist, bins

	@staticmethod
	def peak_clarity(peak_target, hist_list, major_bin, major_sparsity, minor_bin=None, minor_sparsity=None):
		"""
		Function to determine how clean a peak is (homo/hetero)
		Look at densities around each supposed peak, and if the value is close then increment a count
		if the count is above a threshold, return False to indicate failure (and raise a flag)
		:param peak_target: hetero/homo
		:param hist_list: list of histogram under investigation
		:param major_bin: bin of hist of major peak
		:param major_sparsity: sparsity value of that
		:param minor_bin: bin of hist of minor peak
		:param minor_sparsity: sparsity value of that
		:return: True/False
		"""
		clarity_count = 0
		if peak_target == 1:
			major_slice = hist_list[int(major_bin) - 2:int(major_bin) + 2]
			for density in major_slice:
				if np.isclose(major_sparsity, density):
					clarity_count += 1
			if clarity_count > 3:
				return False
		if peak_target == 2:
			major_slice = hist_list[int(major_bin) - 2:int(major_bin) + 2]
			minor_slice = hist_list[int(minor_bin) - 2:int(minor_bin) + 2]
			for density in major_slice:
				if np.isclose(major_sparsity, density):
					clarity_count += 1
			for density in minor_slice:
				if np.isclose(minor_sparsity, density):
					clarity_count += 1
			if clarity_count > 5:
				return False
		return True

	def density_estimation(self, plot_flag):
		"""
		Denisity estimate for a given input distribution (self.input_distribution)
		Use KDE to determine roughly where peaks should be located, peak distances, etc
		Plot graphs for visualisation, return information to origin of call
		:param plot_flag: do we plot a graph or not? (CCG:Yes,CAG:No)
		:return: {dictionary of estimated attributes for this input}
		"""

		##
		## Set up variables for this instance's run of density estimation
		## and generate a dictionary to be modified/returned
		distro_list = list(self.input_distribution)
		major_estimate = None; minor_estimate = None
		peak_distance = None; peak_threshold = None
		estimated_attributes = {'PrimaryPeak':major_estimate,
								'SecondaryPeak':minor_estimate,
								'PeakDistance':peak_distance,
								'PeakThreshold':peak_threshold}

		##
		## Begin density estimation!
		## By default, runs in heterozygous assumption
		## If instance requires homozygous then tailor output instead of re-running
		major_estimate = max(self.input_distribution); major_index = distro_list.index(major_estimate)
		minor_estimate = max(n for n in distro_list if n!=major_estimate); minor_index = distro_list.index(minor_estimate)

		##
		## Check that N-1 of <MAJOR> is not <MINOR> (i.e. slippage)
		## If so, correct for minor == 3rd highest and raise error flag!
		if minor_index == major_index-1:
			literal_minor_estimate = max(n for n in distro_list if n!=major_estimate and n!=minor_estimate)
			literal_minor_index = distro_list.index(literal_minor_estimate)
			minor_estimate = literal_minor_estimate
			minor_index = literal_minor_index
			if self.contig_stage == 'CCG': self.CCGExpansionSkew = True
			if self.contig_stage == 'CAG': self.CAGExpansionSkew = True

		##
		## Actual execution of the Kernel Density Estimation histogram
		hist, bins = self.histogram_generator(self.filename, self.graph_title, self.axes, plot_flag)
		hist_list = list(hist)

		##
		## Determine which bin in the density histogram our estimate values reside within
		## -1 because for whatever reason np.digitize adds one to the literal index
		major_estimate_bin = np.digitize([major_estimate], bins)-2
		minor_estimate_bin = np.digitize([minor_estimate], bins)-1

		##
		## Relevant densities depending on zygosity of the current sample
		major_estimate_sparsity = None; minor_estimate_sparsity = None
		if self.peak_target == 1:
			major_estimate_sparsity = min(n for n in hist if n!=0)
			minor_estimate_sparsity = min(n for n in hist if n!=0)
			peak_distance = 0
			if not self.peak_clarity(self.peak_target, hist_list, major_estimate_bin, major_estimate_sparsity):
				if self.contig_stage == 'CCG': self.CCGPeakAmbiguous = True
				if self.contig_stage == 'CAG': self.CAGPeakAmbiguous = True
		if self.peak_target == 2:
			major_estimate_sparsity = min(n for n in hist if n!=0)
			minor_estimate_sparsity = min(n for n in hist if n!=0 and n!=major_estimate_sparsity)
			peak_distance = np.absolute(major_index - minor_index)
			if not self.peak_clarity(self.peak_target, hist_list, major_estimate_bin, major_estimate_sparsity, minor_estimate_bin, minor_estimate_sparsity):
				if self.contig_stage == 'CCG': self.CCGPeakAmbiguous = True
				if self.contig_stage == 'CAG': self.CAGPeakAmbiguous = True

		##
		## Check for multiple low densities in distribution
		fuzzy_count = 0
		for density in hist_list:
			if np.isclose(major_estimate_sparsity, density):
				fuzzy_count+=1
		if fuzzy_count > 3:
			if self.contig_stage == 'CCG': self.CCGDensityAmbiguous = True
			if self.contig_stage == 'CAG': self.CAGDensityAmbiguous = True

		##
		## Determine Thresholds for this instance sample
		peak_threshold = 0.50
		if self.contig_stage == 'CCG':
			if self.CCGExpansionSkew:
				if self.peak_target == 2:
					peak_threshold -= 0.05
			if self.CCGDensityAmbiguous:
				peak_threshold -= 0.075
			if self.CCGPeakAmbiguous:
				peak_threshold -= 0.10

		if self.contig_stage == 'CAG':
			if self.CAGExpansionSkew:
				if self.peak_target == 2:
					peak_threshold -= 0.05
			if self.CAGDensityAmbiguous:
				peak_threshold -= 0.075
			if self.CAGPeakAmbiguous:
				peak_threshold -= 0.10

		##
		## Preparing estimated attributes for return
		if self.peak_target == 1:
			estimated_attributes['PrimaryPeak'] = major_index+1
			estimated_attributes['SecondaryPeak'] = major_index+1
		if self.peak_target == 2:
			estimated_attributes['PrimaryPeak'] = major_index+1
			estimated_attributes['SecondaryPeak'] = minor_index+1
		estimated_attributes['PeakDistance'] = peak_distance
		estimated_attributes['PeakThreshold'] = peak_threshold

		return estimated_attributes

	def differential_peaks(self, first_pass, fod_params, threshold_bias, fail_state=False, fod_recall=False):
		"""
		Function which takes in parameters gathered from density estimation
		and applies them to a First Order Differential peak detection algorithm
		to more precisely determine the peak (and thus, genotype) of a sample
		:param first_pass: Dictionary of results from KDE
		:param fod_params: Parameters for graphs made in this function
		:param threshold_bias: Bool for whether this call is a re-call or not (lower threshold if True)
		:param fail_state: did this FOD fail or not?
		:param fod_recall: do we need to do a local re-call?
		:return: dictionary of results from KDE influenced FOD
		"""

		##
		## Get Peak information from the KDE dictionary
		## If threshold_bias == True, this is a recall, so lower threshold
		## but ensure the threshold stays within the expected ranges
		peak_distance = first_pass['PeakDistance']
		peak_threshold = first_pass['PeakThreshold']
		if threshold_bias or fod_recall:
			self.RecallCount+=1
			if self.RecallCount > 5:
				raise AttributeError('Re-called 5+ times. Unable to accurately predict genotype.')
			first_pass['PeakThreshold'] -= 0.075
			peak_threshold -= 0.075
			peak_threshold = max(peak_threshold,0.05)

			first_pass['PeakDistance'] += 1
			peak_distance += 1
			peak_distance = max(peak_distance, 4)

		self.ThresholdUsed = peak_threshold

		##
		## Graph Parameters expansion
		linspace_dimensionality = fod_params[0]
		graph_title = fod_params[1]
		axes = fod_params[2]
		filename = fod_params[3]

		##
		## Create planar space for plotting
		## Send paramters to FOD
		## Increment results by 1 (to resolve 0 indexing)
		x = np.linspace(linspace_dimensionality[0], linspace_dimensionality[1], linspace_dimensionality[2])
		buffered_y = np.asarray([0] + list(self.input_distribution))
		y = self.input_distribution
		peak_indexes = peakutils.indexes(y, thres=peak_threshold, min_dist=peak_distance-1)
		fixed_indexes = np.array(peak_indexes+1)

		##
		## Check that we didn't get too many peaks..
		if self.contig_stage == 'CCG':
			if len(fixed_indexes) > 2:
				self.CCGPeakOOB = True
		if self.contig_stage == 'CAG':
			if len(fixed_indexes) > 2:
				self.CAGPeakOOB = True

		##
		## Plot Graph!
		## Set up dimensions for plotting
		plt.figure(figsize=(10,6))
		plt.title(graph_title)
		plt.xlabel(axes[0])
		plt.ylabel(axes[1])

		##
		## Set appropriate range size for CCG/CAG graph dimension
		if self.contig_stage == 'CCG':
			plt.xticks(np.arange(0,21,1))
			plt.xlim(1,20)
		if self.contig_stage == 'CAG':
			plt.xticks(np.arange(0,201,50))
			plt.xlim(1,200)

		##
		## Try to assign peaks to the appropriate indexes
		## If we're expecting one peak (CCG Het) then it's simple to do so..
		## Index Error = too few peaks were called, so we can fail and re-call
		if self.peak_target == 1:
			try:
				first_pass['PrimaryPeak'] = fixed_indexes.item(0)
				first_pass['SecondaryPeak'] = fixed_indexes.item(0)
			except IndexError:
				fail_state = True
		## However, if we're expecting 2 peaks (CCG Hom) then it's more complicated
		## If there are not enough peaks called (1 instead of 2) then the current sample may fall into:
		## -- Homozygous Haplotype.. a true homozygote; CAGxCCGyCAGxCCGy
		## -- Diminished Peak.. there are two peaks but one is significantly smaller than the main, and
		##    lowering threshold to detect it would just introduce noise / worsen results
		## -- Neighbouring Peak.. there are two peaks within this sample, they are just right next to each other
		## And thus the following deterministic functions figure out which rare-case our sample fits into
		if self.peak_target == 2 and self.contig_stage == 'CAG':
			try:
				first_pass['PrimaryPeak'] = fixed_indexes.item(0)
				first_pass['SecondaryPeak'] = fixed_indexes.item(1)
			except IndexError:
				haplotype_failure, interp_dist = self.haplotype_deterministic([x,y,peak_indexes])
				if haplotype_failure:
					diminished_failure, diminished_peaks = self.diminished_deterministic([x,y,buffered_y,peak_indexes])
					if diminished_failure:
						neighbour_failure, neighbour_peaks = self.neighbour_deterministic([x,y,buffered_y,peak_indexes])
						if neighbour_failure:
							fail_state = True
						else:
							self.NeighbouringPeaks = True
							first_pass['PrimaryPeak'] = neighbour_peaks.item(0)
							first_pass['SecondaryPeak'] = neighbour_peaks.item(1)
					else:
						self.DiminishedPeak = True
						first_pass['PrimaryPeak'] = diminished_peaks.item(0)
						first_pass['SecondaryPeak'] = diminished_peaks.item(1)
				else:
					self.PotentialHomozygousHaplotype = True
					self.PHHInterpDistance = interp_dist
					first_pass['PrimaryPeak'] = fixed_indexes.item(0)
					first_pass['SecondaryPeak'] = fixed_indexes.item(0)

		##
		## CCG search dimensions are so small that the only issue ever will be neighbouring peaks
		## So if we expect 2 peaks, and don't have 2, pass to neighbouring_deterministic()
		## If that fails.. CCG unsalvagable and thus CAG can't be acquired -- sample failure
		if self.peak_target == 2 and self.contig_stage == 'CCG':
			try:
				first_pass['PrimaryPeak'] = fixed_indexes.item(0)
				first_pass['SecondaryPeak'] = fixed_indexes.item(1)
			except IndexError:
				neighbour_failure, neighbour_peaks = self.neighbour_deterministic([x,y,buffered_y,peak_indexes])
				if neighbour_failure:
					fail_state = True
				else:
					self.NeighbouringPeaks = True
					first_pass['PrimaryPeak'] = neighbour_peaks.item(0)
					first_pass['SecondaryPeak'] = neighbour_peaks.item(1)

		##
		## Re-create indexes incase that we had a haplotype/neighbouring flag issue
		fixed_indexes = np.array([first_pass['PrimaryPeak'], first_pass['SecondaryPeak']])

		##
		## Execute actual plotting last, incase of homozyg haplotype/neighbouring peaks
		## Plot graph and identified peaks; label appropriately based on size of fixed_indexes
		pplot(x,buffered_y,fixed_indexes)
		if fixed_indexes.item(0) == fixed_indexes.item(1): plt.legend(['Genotype: {}'.format(fixed_indexes.item(0))])
		else: plt.legend(['Genotype: {},{}'.format(fixed_indexes.item(0),fixed_indexes.item(1))])
		plt.savefig(os.path.join(self.prediction_path,filename), format='png')
		plt.close()

		return fail_state, first_pass

	def haplotype_deterministic(self, peak_info, fail_state=False):
		"""
		Function to determine whether a potential heterozygous haplotype is legitimate
		Peak Clarity... sliding window around peak && comparison of raw read counts to gauge
		whether a single peak in an expected heterozygous context is legitimate
		:param peak_info: [x, y, genotype]
		:param fail_state: whether we pass/fail this deterministic
		:return: fail_state (boolean for success of homozygous haplotype), interp_distance (gaussian interp peak distance to prediction)
		"""
		##
		## Calculate the percentage drop-off of our suspected peak
		## If the thresholds are meant, do an interpolation on the peak
		nmt = self.input_distribution[peak_info[2]-2]
		nmo = self.input_distribution[peak_info[2]-1]
		n = self.input_distribution[peak_info[2]]
		npo = self.input_distribution[peak_info[2]+1]
		npt = self.input_distribution[peak_info[2]+2]
		nmt_n = nmt / n; nmo_n = nmo / n
		npo_n = npo / n; npt_n = npt / n

		##
		## For N Minus/Plus One/Two, determine whether the dropoff is acceptable for a clean peak
		## If it is, add to a pass total; if 3/4 tests pass, peak is clean enough to interpolate with confidence
		pass_total = 0
		point_tests = [[[nmt_n],[0.35],[0.150]],
					   [[nmo_n],[0.65],[0.075]],
					   [[npo_n],[0.50],[0.300]],
					   [[npt_n],[0.25],[0.250]]]
		for test in point_tests:
			if np.isclose(test[0],test[1],atol=test[2]):
				pass_total += 1
		##
		## When the passrate was 3/4, we can interpolate our peak with confidence
		## Attempt to fit a guassian near our suspected haplotype peak
		interp_distance = 0.0
		if pass_total >= 3:
			peaks_interp = peakutils.interpolate(peak_info[0], peak_info[1], ind=peak_info[2])
			if np.isclose([peaks_interp],[float(peak_info[2])], atol=1.25):
				interp_distance = abs(peaks_interp - float(peak_info[2]))
				pass
			else:
				fail_state = True
		else:
			fail_state = True

		##
		## Return whether we think this is a homozygous haplotype or not
		return fail_state, interp_distance

	def diminished_deterministic(self, peak_info, fail_state=False):
		"""
		Function to do a low-pass inspection along the entire distribution, in chunks
		Allows us to gauge whether this input distribution has a secondary peak which is considerably much lower
		than the current main peak prediction in terms of literal read count value, but is still definitely a peak
		within it's own right. (e.g. main peak = 20k reads, sub peak = 500). We do this as continually lowering
		the threshold of the main peak calling algorithm will just introduce noise beyond a certain point; this
		function allows us to comb over the details of the entire distribution and ensure that we get the correct
		second peak...
		:param peak_info: list of peak information to be utilised within this investigation
		:param fail_state: whether or not this function hit the fail conditions
		:return: fail_state, new_peaks [x, y]
		"""

		##
		## Get the relevant information for this investigation
		n = self.input_distribution[peak_info[3]]
		peak_resolutions = [[25,8,6],[40,5,3]]
		resolution_peaks = []

		##
		## Loop over two different slicing contexts; to ensure we don't get too many peaks in one slice
		for resolution in peak_resolutions:
			##
			## Slice the distribution we're working with so we can check lower threshold areas specifically
			sliced_distribution = np.split(self.input_distribution, resolution[0])
			slice_index_dict = {}
			for i in range(0,len(sliced_distribution)):
				slice_indexes = peakutils.indexes(sliced_distribution[i], thres=0.025, min_dist=resolution[2])
				slice_index_dict[i+1] = {'SlicePeakIndex':slice_indexes,
										 'SlicePeakValue':sliced_distribution[i][slice_indexes],
										 'SliceDistro':sliced_distribution[i]}

			##
			## Loop over every low-pass distribution slice to determine the 'candidate' for the diminished peak
			## if the current Slice Peak Value is > the previously seen max candidate, and not = n (real peak)
			## then current maximum candidate for diminished peak = current Slice Peak Value
			## Set objects accordingly..
			current_candidate = 0
			current_candidate_slice = 0
			current_candidate_index = 0
			for dkey, ddata in sorted(slice_index_dict.iteritems()):
				if ddata['SlicePeakValue'] > current_candidate and ddata['SlicePeakValue'] != n:
					current_candidate = ddata['SlicePeakValue']
					current_candidate_slice = int(dkey)
					current_candidate_index = ddata['SlicePeakIndex']

			##
			## Calculate the literal (self.input_distribution based) index of this supposed candidate diminished peak
			## Then return an array of new peaks, after checking the index int value is in the right order (low, high)..
			subpeak_literal = np.array(((current_candidate_slice-1) * resolution[1]) + current_candidate_index)
			if int(subpeak_literal) < int(peak_info[3]):
				new_peaks = np.array([subpeak_literal.item(0), peak_info[3].item(0)])
			else:
				new_peaks = np.array([peak_info[3].item(0), subpeak_literal.item(0)])

			##
			## If the two resolutions produced a differing prediction.. check that the values make sense
			## I.E. even if a resolution detected a peak, ensure that the literal read count value is likely (T1/2/3/4)
			topone = max(self.input_distribution)
			toptwo = max(n for n in self.input_distribution if n!=topone)
			topthree = max(n for n in self.input_distribution if n!=topone and n!=toptwo)
			topfour = max(n for n in self.input_distribution if n!=topone and n!=toptwo and n!=topthree)

			##
			## Once the top4 have been discerned; loop over the current resolution's predicitons
			## If both peaks in this allele are determined to confine within the desired logic, score 2
			## If one, score 1; if none, score 0 (fail)
			peak_register_score = 0
			for value in new_peaks:
				if self.input_distribution[value] not in [topone, toptwo, topthree, topfour]:
					list(new_peaks).remove(value); continue
				else: peak_register_score += 1
			if peak_register_score == 0: fail_state = True
			else: resolution_peaks.append([new_peaks, peak_register_score])

		##
		## Now that we have the two slice resolution peak predictions; check and take the best scoring one
		## If the scores are both 1... hmm.. raise flag for inspection
		current_resolution_candidate = None
		for i in range(0,1):
			try:
				firstreso_peak = resolution_peaks[i][0]; firstreso_score = resolution_peaks[i][1]
				secndreso_peak = resolution_peaks[i+1][0]; secndreso_score = resolution_peaks[i+1][1]
				if firstreso_score > secndreso_score: current_resolution_candidate = firstreso_peak
				if secndreso_score > firstreso_score: current_resolution_candidate = secndreso_peak
				if firstreso_score == secndreso_score:
					current_resolution_candidate = firstreso_peak
					self.DiminishedUncertainty = True
			except IndexError:
				pass

		##
		## Rudimentary fail state catcher so that we don't force a neighbouring peak
		## into a diminished peak determination..
		## If the literal read count value is so small that it couldn't possibly be a diminished peak
		## Then we fail.. pass onto neighbouring deterministic
		for peak in current_resolution_candidate:
			if np.isclose([list(self.input_distribution)[peak]], [0], atol=25):
				fail_state = True

		##
		## Return our fail_state and new peaks array
		## Adding one to the peaks to resolve inner (0-index) and literal(1-index) issue
		new_peaks = np.array(current_resolution_candidate + 1)
		return fail_state, new_peaks

	def neighbour_deterministic(self, peak_info, fail_state=False):
		"""
		Function to determine whether a sample which was not a homozygous haplotype is infact
		a heterozygous normal, with target peaks being neighbours to N; either N-1 or N+1
		However, before assuming we are in a neighbouring situation, we do a low-read count pass
		of the distribution incase there are any legitimate discrete peaks which are very small
		compared to the major peak; this is the real 'secondary' peak and will be used if found.
		Otherwise, we follow the neighbouring peak deterministic...
		:param peak_info: [x, y, buffered_y, genotype]
		:param fail_state: whether we pass or fail..
		:return: fail_state (boolean for success), new_peaks (determined neighbour or mini-peak)
		"""

		##
		## Get the relevant information for this investigation
		nmo = self.input_distribution[peak_info[3]-1]
		n = self.input_distribution[peak_info[3]]
		npo = self.input_distribution[peak_info[3]+1]
		new_peaks = np.array([0,0])

		##
		## Discrete check on raw read value of bins
		## Return arrays (indexes corrected for 0-indexing of np, but 1-indexing of real)
		if nmo > npo:
			new_peaks = np.array([peak_info[3],peak_info[3]+1])
		elif npo > nmo:
			new_peaks = np.array([peak_info[3]+1,peak_info[3]+2])
			if n > npo:
				self.CAGBackwardsSlippage = True
		else:
			fail_state = True ## if here.. sample uncallable

		##
		## Return our new genotypes and failstate
		return fail_state, new_peaks

	def get_warnings(self):
		"""
		Function which generates a dictionary of warnings encountered in this instance of SequenceTwoPass
		Dictionary is later sorted into the GenotypePrediction equivalency for returning into a report file
		:return: {warnings}
		"""

		return {'CCGDensityAmbiguous':self.CCGDensityAmbiguous,
				'CCGExpansionSkew':self.CCGExpansionSkew,
				'CCGPeakAmbiguous':self.CCGPeakAmbiguous,
				'CCGPeakOOB':self.CCGPeakOOB,
				'PotentialHomozygousHaplotype':self.PotentialHomozygousHaplotype,
				'PHHInterpDistance':self.PHHInterpDistance,
				'NeighbouringPeaks':self.NeighbouringPeaks,
				'DiminishedPeak':self.DiminishedPeak,
				'DiminishedUncertainty':self.DiminishedUncertainty,
				'CAGDensityAmbiguous':self.CAGDensityAmbiguous,
				'CAGExpansionSkew':self.CAGExpansionSkew,
				'CAGPeakAmbiguous':self.CAGPeakAmbiguous,
				'CAGPeakOOB':self.CAGPeakOOB,
				'CAGBackwardsSlippage':self.CAGBackwardsSlippage,
				'CAGForwardSlippage':self.CAGForwardSlippage,
				'ThresholdUsed':self.ThresholdUsed,
				'RecallCount':self.RecallCount}

class MosaicismInvestigator:
	def __init__(self, genotype, distribution):
		"""
		A class which is called when the functions within are required for somatic mosaicism calculations
		As of now there is only a basic implementation of somatic mosaicism studies but it's WIP
		"""

		self.genotype = genotype
		self.distribution = distribution

	def chunks(self, n):
		"""
		Function which takes an entire sample's distribution (200x20) and split into respective 'chunks'
		I.E. slice one distribution into contigs for each CCG (200x1 x 20)
		:param n: number to slice the "parent" distribution into
		:return: CHUNKZ
		"""

		for i in xrange(0, len(self.distribution), n):
			yield self.distribution[i:i + n]

	@staticmethod
	def arrange_chunks(ccg_slices):
		"""
		Function which takes the sliced contig chunks and orders them into a dataframe for ease of
		interpretation later on in the application. Utilises pandas for the dataframe class since
		it's the easiest to use.
		:param ccg_slices: The sliced CCG contigs
		:return: df: a dataframe which is ordered with appropriate CCG labels.
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
		Function to take specific CCG contig sub-distribution from dataframe
		Extract appropriate N-anchored values for use in sommos calculations
		:param df: input dataframe consisting of all CCG contig distributions
		:param input_allele: genotype derived from GenotypePrediction (i.e. scrape target)
		:return: allele_nvals: dictionary of n-1/n/n+1
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
		Function to execute the actual calculations
		Perhaps float64 precision is better? Might not matter for us
		Also required to add additional calculations here to make the 'suite' more robust
		:param allele_values: dictionary of this sample's n-1/n/n+1
		:return: dictionary of calculated values
		"""

		nmo = allele_values['NMinusOne']
		n = allele_values['NValue']
		npo = allele_values['NPlusOne']
		nmo_over_n = 0
		npo_over_n = 0

		try:
			nmo_over_n = int(nmo) / int(n)
			npo_over_n = int(npo) / int(n)
		except ValueError:
			nmo_over_n = float(nmo) / float(n)
			npo_over_n = float(npo) / float(n)
		except ZeroDivisionError:
			log.info('{}{}{}{}'.format(clr.red,'shd__ ',clr.end,' Divide by 0 attempted in Mosaicism Calculation.'))

		calculations = {'NMinusOne':nmo,'NValue':n,'NPlusOne':npo,'NMinusOne-Over-N': nmo_over_n, 'NPlusOne-Over-N': npo_over_n}
		return calculations

	@staticmethod
	def distribution_padder(ccg_dataframe, genotype):
		"""
		Function to ensure all distribution's N will be anchored to the same position in a file
		This is to allow the end user manual insight into the nature of the data (requested for now)
		E.G. larger somatic mosaicism spreads/trends in a distribution + quick sample-wide comparison
		:param ccg_dataframe: dataframe with all CCG contigs
		:param genotype: genotype derived from GenotypePrediction i.e. scrape target
		:return: distribution with appropriate buffers on either side so that N is aligned to same position as all others
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