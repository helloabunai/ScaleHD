##
## Generics
import os
import sys
import logging as log

##
## Backend junk
from ..backpack import Colour as clr
from ..align.alignment import extract_repeat_distributions

class SeqPredict:

	def __init__(self, instance_rundir, instance_params, model_data, model_descriptor, assembly_files=None, distribution_files=None):

		self.instance_rundir = instance_rundir
		self.instance_params = instance_params
		self.assembly_files = assembly_files
		self.distribution_files = distribution_files
		self.model_data = model_data
		self.model_desc = model_descriptor

		##TODO build model once, call upon when required (rather than loading for each unlabelled sample)

		if self.assembly_files is not None:
			self.initialise_sam_output()
		if distribution_files is not None:
			self.initialise_csv_output()

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

			##TODO call genotyping functions here

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