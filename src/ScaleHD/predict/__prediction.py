from __future__ import division

#/usr/bin/python
__version__ = 0.01
__author__ = 'alastair.maxwell@glasgow.ac.uk'

##
## Generic imports
import sys
import csv
import numpy as np
import logging as log
from sklearn import svm
from sklearn import preprocessing
from sklearn.multiclass import OneVsOneClassifier

##
## Backend Junk
from ..__backend import Colour as clr
from ..__backend import DataLoader

PRD_REPORT = []

class GenotypePrediction:
	def __init__(self, data_pair, prediction_path, training_data, instance_params):

		"""
		Prediction stage of the pipeline -- use of SVM and other applications to automatically
		determine the genotype of a sample. Utilise forward distribution for CAG information,
		utilise reverse reads for CCG/intervening information

		General workflow of functions within this class:

		-- Start with CCG determination
		-- Regardless of origin, scrape CCG distribution, then sum each CAG for every CCG aka collapse
		-- Once collapsed, normalise and feed into SVM model to determine HOMO/HETEROZYGOUS
			-- If disconnect between forward/reverse zygosity, attempt brute force determination
		-- From there, tailor behaviour to latter data

		Still a work in progress...

		"""

		##
		## Paths to files/data
		self.data_pair = data_pair
		self.prediction_path = prediction_path
		self.training_data = training_data
		self.instance_params = instance_params

		##
		## Arrays of this sample pair's CCG (collapsed)
		## Classifier used for predicting these arrays
		self.classifier, self.encoder = self.build_predictive_model()

		##
		## Determine zygosity of input distributions
		self.collapsed_forward_ccg = self.collapse_ccg(self.scrape_distro(self.data_pair[0]))
		self.collapsed_reverse_ccg = self.collapse_ccg(self.scrape_distro(self.data_pair[1]))
		self.ccgzygosity = self.predict_ccg_zygstate()

		##
		## Confidence scalars (subject to change)
		## Disconnect bool will act as negative confidence anchor
		## in the event that zygosity information is executed after disconnect analysis...
		self.fwdist_spread_scalar = 0
		self.rvdist_spread_scalar = 0
		self.prediction_confidence = 0
		self.disconnect_bool = False

		##
		## We've determined the zygstate, so we can tailor our behaviour from here
		self.workflow_guide()

	def workflow_guide(self):

		"""
		Flow function which directs behaviour based on zygosity of CCG
		-- removed from __init__ as required to be re-called in disconnect status
		"""

		if self.ccgzygosity == 'HOMO':
			self.homozygous_workflow()
		if self.ccgzygosity == 'HETERO':
			self.heterozygous_workflow()
		if self.ccgzygosity == 'DISCONNECT':
			self.disconnect_bool = True
			self.analyse_disconnect()

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
	def collapse_ccg(distribution_array):

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

	def build_predictive_model(self):

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

	def predict_ccg_zygstate(self):

		##
		## Reshape the input distributions so skl doesn't complain about 1D
		## normalise too, as in tandem with training data -- requires casting to float64
		forward_reshape = preprocessing.normalize(np.float64(self.collapsed_forward_ccg.reshape(1,-1)))
		reverse_reshape = preprocessing.normalize(np.float64(self.collapsed_reverse_ccg.reshape(1,-1)))

		##
		## Predict the zygstate, then decode from hash into literal label value
		forward_zygstate = str(self.encoder.inverse_transform(self.classifier.predict(forward_reshape)))
		reverse_zygstate = str(self.encoder.inverse_transform(self.classifier.predict(reverse_reshape)))

		if not forward_zygstate == reverse_zygstate:
			print 'FW: ', forward_zygstate
			print 'RV: ', reverse_zygstate
			return 'DISCONNECT'
		else:
			return forward_zygstate[2:-2]

	def homozygous_workflow(self):

		##
		## Working on a homozygous patient sample, thus...
		print 'Homozygous function...'
		print self.data_pair[0].split('/')[7]
		print self.collapsed_forward_ccg
		print self.collapsed_reverse_ccg

	def heterozygous_workflow(self):

		##
		## Working on a heterozygous patient sample...
		print 'Heterozygous function...'
		print self.data_pair[0].split('/')[7]
		print self.collapsed_forward_ccg
		print self.collapsed_reverse_ccg

	def analyse_disconnect(self):

		def get_indices(input_distribution, n):

			index = input_distribution.tolist().index(n)
			index_nminus = index-1
			index_nplus = index+1

			return index, index_nminus, index_nplus

		def check_peak_distance(maxpeak_index, subpeak_index):

			distance = maxpeak_index-subpeak_index
			if distance < 0:
				pass
			else:
				return distance

		def edit_minor_estimate(input_distribution, max_peak, sub_peak):

			sub_estimate = max(n for n in input_distribution if (n!=max_peak and n!=sub_peak))
			sub_index, sub_nminus, sub_nplus = get_indices(input_distribution, sub_estimate)
			return sub_estimate, sub_index, sub_nminus, sub_nplus

		def percentage_shift(input_distribution, n_minus=None, n=None, n_plus=None, distro_mean=None):

			##
			## If distro mean is present, we're calculating distro_mean to N
			if not distro_mean is None:

				n_val = input_distribution[n]
				distro_mean = distro_mean

				mean_to_n = ((n_val - distro_mean) / n_val) * 100
				return mean_to_n

			##
			## Otherwise, we're calculating the percentage shift from n-1/+1 to N
			else:
				n_val = input_distribution[n]
				nm1_val = input_distribution[n_minus]
				np1_val = input_distribution[n_plus]

				nm1_to_n = ((n_val - nm1_val) / n_val) * 100
				np1_to_n = ((n_val - np1_val) / n_val) * 100

				return nm1_to_n, np1_to_n

		##
		## Working on a sample with two 'different' zygstates
		## manually check the distributions to see what's up
		print 'Disconnect analysis function...'
		print self.data_pair[0].split('/')[7]

		##
		## Model predicted that zygstate == HOMO
		## Manual estimate check..
		for distribution in [self.collapsed_forward_ccg, self.collapsed_reverse_ccg]:

			print distribution

			##
			## Highest and second highest value in distribution
			major_estimate = max(distribution)
			minor_estimate = max(n for n in distribution if n!=major_estimate)

			##
			## Indexes for N-1, N and N+1 of each 'peak'
			major_index, major_nminus, major_nplus = get_indices(distribution, major_estimate)
			minor_index, minor_nminus, minor_nplus = get_indices(distribution, minor_estimate)

			##
			## If the second peak is closer to origin in distro index
			## !!POSSIBLE HETEROZYGOTE WITH LOW READS AROUND FIRST PEAK!!
			if major_index > minor_index:

				##
				## Get intrapeak distance...
				## If the peaks are neighbours (1), check minor_estimate is not n-1 of major_estimate
				## If it is n-1 of major_estimate, modify to the third value in distro.. (assume)
				peak_distance = check_peak_distance(major_index, minor_index)
				if peak_distance == 1:
					if minor_index == major_nminus:
						minor_estimate, minor_index, minor_nminus, minor_nplus = edit_minor_estimate(distribution, major_estimate, minor_estimate)

				##
				## Now assuming that the major and minor estimates are ""OK"" at best
				## Look at spread surrounding these peaks for discrete properties
				print 'Peak1: ', major_estimate, ' -- N-1/N/N+1 (0 indexed): ', major_nminus, major_index, major_nplus
				print 'Peak2: ', minor_estimate, ' -- N-1/N/N+1 (0 indexed): ', minor_nminus, minor_index, minor_nplus

				peak_distance = check_peak_distance(major_index, minor_index)
				if peak_distance == 1:
					if minor_index == major_nminus:
						minor_estimate, minor_index, minor_nminus, minor_nplus = edit_minor_estimate(distribution, major_estimate, minor_estimate)

				##
				## Take mean of distro and compare mean to N
				distro_meanval = np.mean(distribution)
				print 'Variance: ', distro_meanval
				major_mean_change = percentage_shift(distribution, n=major_index, distro_mean=distro_meanval)
				minor_mean_change = percentage_shift(distribution, n=minor_index, distro_mean=distro_meanval)
				print 'Maj mean/n: ', major_mean_change
				print 'Min mean/n: ', minor_mean_change

				##
				## Look at n-1 and n+1 and gauge what % they are of N
				major_posterior_change, major_anterior_change = percentage_shift(distribution, n_minus=major_nminus, n=major_index, n_plus=major_nplus)
				minor_posterior_change, minor_anterior_change = percentage_shift(distribution, n_minus=minor_nminus, n=minor_index, n_plus=minor_nplus)
				print 'Maj pos/ant: ', major_posterior_change, ' ', major_anterior_change
				print 'Min pos/ant: ', minor_posterior_change, ' ', minor_anterior_change

				##
				## Percentage change between major/minor
				major_minor_change = ((major_estimate - minor_estimate) / major_estimate) * 100
				print 'Maj to min: ', major_minor_change

				##
				## Determine behaviour based on these results..
				## Different threshold priorities results in flow change
				## 1) overrall spread of distro
				## 2) pcnt of mean to 'peaks'
				## 3) pcnt of n-1/n+1 to 'peaks'
				## 4) failure

				##TODO HEHE WRITE THAT
				##TODO adjust scalars, re-call workflow_guide()

			##
			## maybe check spread around the points.. just for confidence
			if not major_index > minor_index:
				print 'Assuming distribution OK'
				##TODO from red-lined function here, but less confidence
				##TODO adjust scalars, re-call workflow_guide()

			print '\n'

def get_predictionreport():
	return PRD_REPORT