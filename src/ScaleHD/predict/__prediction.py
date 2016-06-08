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
		## We've determined the zygstate, so we can tailor our behaviour from here
		if self.ccgzygosity == 'HOMO':
			self.homozygous_workflow()
		if self.ccgzygosity == 'HETERO':
			self.heterozygous_workflow()
		if self.ccgzygosity == 'ZYG_DISCONNECT_ERROR':
			PRD_REPORT.append('Zygosity Fw/Rv disconnect error')

	@staticmethod
	def scrape_distro(distribution_file):

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
		## normalise too, as in tandem with training data
		forward_reshape = preprocessing.normalize(np.float64(self.collapsed_forward_ccg.reshape(1,-1)))
		reverse_reshape = preprocessing.normalize(np.float64(self.collapsed_reverse_ccg.reshape(1,-1)))

		##
		## Predict the zygstate, then decode from hash into literal
		forward_zygstate = str(self.encoder.inverse_transform(self.classifier.predict(forward_reshape)))
		reverse_zygstate = str(self.encoder.inverse_transform(self.classifier.predict(reverse_reshape)))

		if not forward_zygstate == reverse_zygstate:
			return 'ZYG_DISCONNECT_ERROR'
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


def get_predictionreport():
	return PRD_REPORT