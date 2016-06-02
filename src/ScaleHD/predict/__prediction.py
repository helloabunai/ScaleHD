#/usr/bin/python
__version__ = 0.01
__author__ = 'alastair.maxwell@glasgow.ac.uk'

##
## Generic imports
import os
import sys
import csv
import numpy as np
import logging as log
from sklearn import svm

##
## Backend Junk
from ..__backend import Colour as clr
from ..__backend import DataLoader

PRD_REPORT = []

class GenotypePrediction:
	def __init__(self, data_pair, prediction_path, training_data, instance_params):
		self.data_pair = data_pair
		self.prediction_path = prediction_path
		self.training_data = training_data
		self.instance_params = instance_params

		self.ccg_collapsed_forward = []
		self.ccg_collapsed_reverse = []
		self.ccgzygosity = self.ccg_zygosity()

		print self.ccgzygosity

		if self.ccgzygosity == 'HOMO':
			self.homozygous_workflow()
		if self.ccgzygosity == 'HETERO':
			self.heterozygous_workflow()

	def ccg_zygosity(self):

		##
		## Workflow of stages for determining zygosity (CCG)
		zygosity_classifier = self.build_classifier()

		##
		## Distribution files for the data_pair we've been passed into an object of GenotypePrediction
		self.ccg_collapsed_forward = self.collapse_ccg(self.scrape_distro(self.data_pair[0]))
		self.ccg_collapsed_reverse = self.collapse_ccg(self.scrape_distro(self.data_pair[1]))

		##
		## Training data for model(s)
		hd_generic_desc = self.training_data['GenericDescriptor']
		collapsed_ccg_data = self.training_data['CollapsedCCGZygosity']

		##
		## Building predictive model(s)
		## Slice return to remove extra formatting from sklearn stdout
		ccgzyg_model = self.build_model(zygosity_classifier, collapsed_ccg_data, hd_generic_desc)
		samplepair_zygosity = self.predict_CCG_zygosity(zygosity_classifier, ccgzyg_model, self.ccg_collapsed_forward, self.ccg_collapsed_reverse)
		return samplepair_zygosity[2:-2]

	@staticmethod
	def build_classifier():

		##
		## Classifier flags will go here when we determine what we want

		cls_SVC = svm.SVC(cache_size=200,coef0=0.0,
						  degree=3,gamma=0.0,kernel='rbf',max_iter=-1,
						  probability=True,random_state=None,
						  shrinking=True,tol=0.001,verbose=0)
		return cls_SVC

	@staticmethod
	def scrape_distro(distribution_file):

		placeholder_array = []
		with open(distribution_file) as csv_file:
			source = csv.reader(csv_file, delimiter=',')
			next(source)  # skipping header
			for row in source:
				placeholder_array.append(int(row[2]))

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

	@staticmethod
	def build_model(classifier, data, desc):

		##
		## Build model with desired input data/to desired classifier
		training_model = DataLoader(data, desc).load_model()

		X = training_model.DATA
		Y = training_model.TARGET
		hash_encoder = training_model.ENCDR
		classifier.fit(X, Y)

		return hash_encoder

	@staticmethod
	def predict_CCG_zygosity(classifier, encoder, forward_reads, reverse_reads):

		"""
		Function which does the actual prediction versus the zyg model
		"""

		hashed_forward_zygstate = classifier.predict(forward_reads)
		hashed_reverse_zygstate = classifier.predict(reverse_reads)

		forward_zygstate = str(encoder.inverse_transform(hashed_forward_zygstate))
		reverse_zygstate = str(encoder.inverse_transform(hashed_reverse_zygstate))

		if not forward_zygstate == reverse_zygstate:
			PRD_REPORT.append('CCG Zygosity difference. Error.')
		else:
			return forward_zygstate

	def homozygous_workflow(self):

		##
		## Working on a homozygous patient sample, thus...
		print self.ccg_collapsed_forward
		print self.ccg_collapsed_reverse

	def heterozygous_workflow(self):
		print 'working on a CCG heterozygous individual'

def get_predictionreport():
	return PRD_REPORT