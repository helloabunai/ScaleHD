##
## Generics
import os
import sys
import csv
import logging as log
import numpy as np
from sklearn import svm

##
## Backend junk
from ..backpack import Colour as clr
from ..backpack import DataLoader
from ..align.alignment import extract_repeat_distributions

class SeqPredict:

	def __init__(self, instance_rundir, instance_params, model_data, model_descriptor, assembly_files=None, distribution_files=None):

		self.instance_rundir = instance_rundir
		self.instance_params = instance_params
		self.assembly_files = assembly_files
		self.distribution_files = distribution_files
		self.model_data = model_data
		self.model_desc = model_descriptor

		log.info('{}{}{}{}'.format(clr.bold,'shd__ ',clr.end,'Loading training vectors and building classifier..'))
		self.classifier = self.build_classifier()
		self.label_encoder = self.build_model()
		log.info('{}{}{}{}'.format(clr.green,'shd__ ',clr.end,'Done!'))

		if self.assembly_files is not None:
			self.initialise_sam_output()
		if distribution_files is not None:
			self.initialise_csv_output()

	@staticmethod
	def distribution_information(distribution_file):

		placeholder_array = []
		with open(distribution_file) as csv_file:
			source = csv.reader(csv_file, delimiter=',')
			next(source)  # skipping header
			for row in source:
				placeholder_array.append(int(row[2]))

		novel_distro = np.array(placeholder_array)
		return novel_distro

	def build_classifier(self):

		try:
			probability_estimate = self.instance_params.config_dict['prediction_flags']['@probability_estimate']
			max_iteration = int(self.instance_params.config_dict['prediction_flags']['@max_iteration'])
		except AttributeError:
			probability_estimate = self.instance_params['probability_estimate']
			max_iteration = int(self.instance_params['max_iteration'])

		cls_SVC = svm.SVC(cache_size=200,coef0=0.0,
						  degree=3,gamma=0.0,kernel='linear',max_iter=max_iteration,
						  probability=bool(probability_estimate),random_state=None,
						  shrinking=True,tol=0.001,verbose=0)
		return cls_SVC

	def build_model(self):

		training_model = DataLoader(self.model_data, self.model_desc).load_model()

		X = training_model.DATA
		Y = training_model.TARGET
		hash_encoder = training_model.ENCDR
		self.classifier.fit(X, Y)

		return hash_encoder

	def initialise_sam_output(self):

		csvdatalen = len(self.assembly_files)
		for i in range(0, len(self.assembly_files)):

			##
			## User feedback if verbose
			log.info('{}{}{}{}{}{}{}'.format(clr.bold,'shd__ ',clr.end,'Predicting genotype from repeat distribution.. ',str(i+1),'/',str(csvdatalen)))

			##
			## If here, only predict stage ran (-i/-b)
			## so no chance of existing pipeline-stage folders in sample output folder
			sample_root = self.assembly_files[i].split('/')[-1].split('.')[0] ##absolutely_disgusting.jpg
			prediction_outdir = os.path.join(self.instance_rundir, sample_root, 'Prediction')
			if not os.path.exists(prediction_outdir): os.makedirs(prediction_outdir)

			##
			## SAM input so distributions haven't been extracted
			## Call on functions from alignment stage which do that
			distribution_file = extract_repeat_distributions(sample_root, prediction_outdir, self.assembly_files[i])
			literal_distribution = self.distribution_information(distribution_file)

			##
			## Genotype prediction
			label_encoder = self.build_model()
			try:
				hash_encoded_prediction = self.classifier.predict(literal_distribution)
			except ValueError:
				log.critical('{}{}{}{}{}{}'.format(clr.red,'shd__ ',clr.end,'Input distribution dimensions invalid! (Input: ', str(len(literal_distribution)) ,' | Training: 2020)'))
				log.critical('{}{}{}{}'.format(clr.red,'shd__ ',clr.end,'Please check your input sequence assembly. Aligned to wrong reference?'))
				log.critical('{}{}{}{}'.format(clr.red,'shd__ ',clr.end,'Cannot progress with this data. Exiting.'))
				sys.exit(2)
			normalised_prediction = str(label_encoder.inverse_transform(hash_encoded_prediction))
			decision_function = self.classifier.decision_function(literal_distribution)[0]
			decision_funclist = decision_function.tolist()

			temporary_output = os.path.join(prediction_outdir,'TemporaryPredictionOutput.txt')
			temporary_file = open(temporary_output, 'w')
			temporary_file.write(hash_encoded_prediction)
			temporary_file.write(normalised_prediction)
			temporary_file.write(decision_funclist)
			temporary_file.close()

			sys.stdout.flush()

	def initialise_csv_output(self):

		csvdatalen = len(self.distribution_files)
		for i in range(0, len(self.distribution_files)):

			##
			## User feedback if verbose
			log.info('{}{}{}{}{}{}{}'.format(clr.bold,'shd__ ',clr.end,'Predicting genotype from repeat distribution.. ',str(i+1),'/',str(csvdatalen)))

			##
			## I/O so we have appropriate sample folder for results subfolder
			sample_root = self.distribution_files[i].split('/')[-3].split('.')[0] ##absolutely_disgusting.jpg
			prediction_outdir = os.path.join(self.instance_rundir, sample_root, 'Prediction')
			if not os.path.exists(prediction_outdir): os.makedirs(prediction_outdir)

			##TODO call genotyping functions here

			sys.stdout.flush()

