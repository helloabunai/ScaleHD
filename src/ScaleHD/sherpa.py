from __future__ import division

#/usr/bin/python
__version__ = 0.01
__author__ = 'alastair.maxwell@glasgow.ac.uk'

##
## Python libraries
import os
import sys
import csv
import gc
import PyPDF2
import argparse
import StringIO
import pkg_resources
import logging as log
from reportlab.pdfgen import canvas
from multiprocessing import cpu_count

##
## Backend junk
from __backend import ConfigReader
from __backend import Colour as clr
from __backend import initialise_libraries
from __backend import mkdir_p
from __backend import sanitise_inputs
from __backend import extract_data
from __backend import sequence_pairings
from __backend import sanitise_outputs
from __backend import collate_peaks
from __backend import generate_atypical_xml
from __backend import generate_reference
from __allelecontainer import SequenceSample

##
## Package stages
from . import seq_qc
from . import align
from . import predict

##
## HTML for instance report
from django.template import Context, Template

##
## Globals
THREADS = cpu_count()

class ScaleHD:
	def __init__(self):
		"""
		ScaleHD: Automated triplet repeat genotyping for Huntington Disease
		ScaleHD has two modes of usage; sequence and batch
		Sequence mode consists of a pipeline behaviour for genome sequence QC, alignment and genotyping
		Batch mode consists of a linear behaviour only for genotyping (from pre-aligned files)
		If you want a full explanation of the ways in which ScaleHD can be run; scalehd --help
		"""

		##
		## Package data
		self.generic_descriptor = pkg_resources.resource_filename(__name__, 'train/long_descr.rst')
		self.collapsed_ccg_zygosity = pkg_resources.resource_filename(__name__, 'train/polyglutamine.csv')
		self.likelihood_matrix = pkg_resources.resource_filename(__name__, 'train/likelihood_matrix.csv')
		self.raw_matrix = pkg_resources.resource_filename(__name__, 'train/raw_matrix.csv')
		self.training_data = {'GenericDescriptor': self.generic_descriptor, 'CollapsedCCGZygosity': self.collapsed_ccg_zygosity}

		##
		## Argument parser from CLI
		self.parser = argparse.ArgumentParser(prog='scalehd', description='ScaleHD: Automated DNA micro-satellite genotyping.')
		self.parser.add_argument('-v', '--verbose', help='Verbose output mode. Setting this flag enables verbose output. Default: off.', action='store_true')
		input_group = self.parser.add_mutually_exclusive_group(required=True)
		input_group.add_argument('-b', '--batch', help='Input batch. Folder of multiple .sam files. For FastQ processing, use -c.', nargs=1)
		input_group.add_argument('-c', '--config', help='Pipeline config. Specify a directory to your ArgumentConfig.xml file.', nargs=1)
		self.parser.add_argument('-t', '--threads', help='Thread utilisation. Typically only alters third party alignment performance. Default: system max.', type=int, choices=xrange(1, THREADS+1), default=THREADS)
		self.parser.add_argument('-p', '--purgesam', help='Remove non-unique reads in sam files the application will work with.', action='store_true')
		self.parser.add_argument('-j', '--jobname', help='Customised folder output name. If not specified, defaults to normal output naming schema.', type=str)
		self.parser.add_argument('-o', '--output', help='Output path. Specify a directory you wish output to be directed towards.', metavar='output', nargs=1, required=True)
		self.args = self.parser.parse_args()

		##
		## Set verbosity for CLI output
		if self.args.verbose:
			log.basicConfig(format='%(message)s', level=log.DEBUG)
			log.info('{}{}{}{}'.format(clr.bold, 'shd__ ', clr.end, 'ScaleHD: Automated DNA micro-satellite genotyping.'))
			log.info('{}{}{}{}'.format(clr.bold, 'shd__ ', clr.end, 'alastair.maxwell@glasgow.ac.uk\n'))
		else:
			log.basicConfig(format='%(message)s')

		##
		## Check inputs, generate outputs
		if sanitise_inputs(self.args):
			log.error('{}{}{}{}'.format(clr.red, 'shd__ ', clr.end, 'Error with specified input(s) configuration. Exiting.'))
			sys.exit(2)
		self.instance_rundir = sanitise_outputs(self.args.jobname, self.args.output)
		self.purge_flag = self.args.purgesam
		self.instance_summary = {}; self.instance_graphs = ''

		##
		## Set up config dictionary of all params.
		## if -c used, read from XML. Else, use 'defaults' in set_params().
		script_path = os.path.dirname(__file__)
		if self.args.config:
			self.configfile = self.args.config[0]
			self.instance_params = ConfigReader(script_path, self.configfile)

		##
		## Placeholder dictionary for use in failure-cases
		self.placeholder_dict = {'PrimaryAllele': 'Fail', 'PrimaryAlleleStatus': 'Fail',
								 'PrimaryAlleleReference': 'Fail', 'PrimaryAlleleOriginal': 'Fail',
								 'SecondaryAllele': 'Fail', 'SecondaryAlleleStatus': 'Fail',
								 'SecondaryAlleleReference': 'Fail', 'SecondaryAlleleOriginal': 'Fail',
								 'PredictionConfidence': 0, 'ForwardDistribution': [0] * 4000,
								 'ReverseDistribution': [0] * 4000}

		##
		## Check third party libraries before continuing
		if initialise_libraries(self.instance_params):
			log.error('{}{}{}{}'.format(clr.red, 'shd__ ', clr.end, 'Detected missing library from system/$PATH. Exiting.'))
			sys.exit(2)
		else:
			log.info('{}{}{}{}'.format(clr.green, 'shd__ ', clr.end, 'Required libraries present, assuming OK!\n'))

		##
		## Set-up instance wide applicable files
		self.index_path = ''; self.reference_indexes = []; self.typical_indexes = []
		self.instance_results = ''; self.instance_matrix = ''; self.instance_graphs = ''
		self.instance_data()

		##
		## Depending on input mode, direct flow of functions
		## -b == multiple files, loop files to class
		## -c == config, do as config parsed flags
		if not self.args.config: self.assembly_workflow()
		else: self.sequence_workflow()

		##
		## Instance wide output results
		## A simple report file is appended after each sample pair, currently..
		## In the future, replace with HTML based web-app, generated here?
		## For now, just exit
		log.info('{}{}{}{}'.format(clr.green, 'shd__ ', clr.end, 'ScaleHD pipeline completed; exiting.'))

	def instance_data(self):

		##
		## Reference indexes
		if self.instance_params.config_dict['instance_flags']['@sequence_alignment']:
			log.info('{}{}{}{}'.format(clr.bold,'shd__ ',clr.end,'Indexing reference(s) before initialising sample pair cycle..'))
			self.index_path = os.path.join(self.instance_rundir,'Indexes'); mkdir_p(self.index_path)
			forward_reference = self.instance_params.config_dict['@forward_reference']
			reverse_reference = self.instance_params.config_dict['@reverse_reference']
			forward_index = align.ReferenceIndex(forward_reference, self.index_path).get_index_path()
			reverse_index = align.ReferenceIndex(reverse_reference, self.index_path).get_index_path()
			self.typical_indexes = [forward_index, reverse_index]
			self.reference_indexes = [forward_index, reverse_index]

		##
		## Instance results (genotype table)
		self.instance_results = os.path.join(self.instance_rundir, 'InstanceReport.csv')
		header = '{},{},{},{},{},{},{},{},{},{},\n'.format('SampleName', 'Primary GTYPE', 'Status', 'Label', 'Original', 'Secondary GTYPE', 'Status', 'Label', 'Original', 'Confidence')
		with open(self.instance_results, 'w') as outfi: outfi.write(header); outfi.close()

		##
		## Instance matrix (sequence distributions)
		self.instance_matrix = os.path.join(self.instance_rundir, 'InstanceMatrix.csv')
		header = ['SampleName']
		for i in range(0,20):
			for j in range(0,200):
				header.append('CAG{}CCG{}'.format(j+1,i+1))
		with open(self.instance_matrix, 'w') as neofi: wr = csv.writer(neofi); wr.writerow(header); neofi.close()

		##
		## Instance graphs
		instance_string = '{} -- Job: {} -- {}'.format('ScaleHD',self.args.jobname,'DateGoesHere')
		packet = StringIO.StringIO()
		can = canvas.Canvas(packet, pagesize=(350, 150))
		can.drawCentredString(175, 75, instance_string)
		can.save()

		packet.seek(0)
		header_page = PyPDF2.PdfFileReader(packet).getPage(0)
		header_output = PyPDF2.PdfFileWriter()
		header_output.addPage(header_page)
		self.instance_graphs = os.path.join(self.instance_rundir, 'InstanceGraphs.pdf')
		header_output.write(file(self.instance_graphs, 'wb'))

	def assembly_workflow(self):
		"""
		Workflow for when ScaleHD is being ran in batch mode..
		I.e. Genotyping only!
		This mode still requires a few methods which were already written in the alignment stage;
		mainly the repeat count distribution method -- it's called here to gather such information
		General overview:
		- For each filepair we need to work on
		- Get paths, extract repeat count distributions
		- Send distributions to genotyping class
		- Append results to report
		:return: None
		"""
		##
		## Input path, create pairs of data in said path
		instance_inputdata = self.args.batch[0]
		assembly_pairs = sequence_pairings(instance_inputdata,self.instance_rundir,'assembly')

		##
		## Temporary report file to be appended to for each instance
		master_results_file = os.path.join(self.instance_rundir, 'InstanceReport.csv')
		header = '{},{},{},{},{},{},{},{},{},{}\n'.format('SampleName', 'Primary GTYPE', 'Status','Label', 'Original',
														  'Secondary GTYPE','Status', 'Label', 'Original', 'Confidence')
		with open(master_results_file, 'w') as outfi: outfi.write(header); outfi.close()

		##
		## Distribution matrix file for each instance
		master_matrix_file = os.path.join(self.instance_rundir,'DistributionMatrix.csv')
		header = ['SampleName']
		for i in range(0,20):
			for j in range(0,200):
				header.append('CAG{}CCG{}'.format(j+1,i+1))
		with open(master_matrix_file, 'w') as neofi: wr = csv.writer(neofi); wr.writerow(header); neofi.close()

		##
		## Executing the workflow for this SHD instance
		for i in range(len(assembly_pairs)):
			for assembly_label, assembly_data in assembly_pairs[i].iteritems():
				log.info('{}{}{}{}{}/{} ({})'.format(clr.bold, 'shd__ ', clr.end, 'Processing assembly pair: ', str(i + 1), str(len(assembly_pairs)),assembly_label))

				##
				## Required data to process
				forward_assembly = assembly_data[0]
				reverse_assembly = assembly_data[1]
				predict_path = assembly_data[2]
				bayesian_path = os.path.join(predict_path,'Bayesian')
				instance_params = self.set_prediction_params()
				placeholder_dict = {'PrimaryAllele': 'Fail', 'PrimaryAlleleStatus': 'Fail',
							   'PrimaryAlleleReference': 'Fail', 'PrimaryAlleleOriginal': 'Fail',
							   'SecondaryAllele': 'Fail', 'SecondaryAlleleStatus': 'Fail',
							   'SecondaryAlleleReference': 'Fail', 'SecondaryAlleleOriginal': 'Fail',
							   'PredictionConfidence': 0, 'ForwardDistribution':[0]*4000, 'ReverseDistribution':[0]*4000}

				##
				## Specific paths to pass to distribution scraper
				## In instance_workflow these would've been created for alignment, but since we don't align here
				## They have to be made in this location instead
				forward_filename = forward_assembly.split('/')[-1].split('.')[0]  ##absolutely_disgusting.jpg
				reverse_filename = reverse_assembly.split('/')[-1].split('.')[0]  ##absolutely_disgusting.jpg
				forward_path = os.path.join(predict_path,forward_filename)
				reverse_path = os.path.join(predict_path,reverse_filename)
				if not os.path.exists(forward_path): os.makedirs(forward_path)
				if not os.path.exists(reverse_path): os.makedirs(reverse_path)

				#############################################
				## Stage One!! Scan for atypical alleles.. ##
				#############################################
				try:
					log.info('{}{}{}{}'.format(clr.bold,'shd__ ',clr.end,'Scanning for atypical alleles..'))
					scanner_object = align.ScanAtypical((forward_path,forward_assembly))
					atypical_count = scanner_object.get_allele_status()
					atypical_info = scanner_object.get_atypical_info()
					processed_atypical, legitimate_atypicals = scrape_atypical(atypical_info)
					if atypical_count != 0: log.info('{}{}{}{}'.format(clr.yellow,'shd__ ',clr.end,'Scanning complete! Atypical allele(s) present.'))
					else: log.info('{}{}{}{}'.format(clr.green, 'shd__ ', clr.end, 'Scanning complete! No atypical allele(s) present.'))
					gc.collect()
				except Exception, e:
					self.write_failure([assembly_label, placeholder_dict, master_results_file, master_matrix_file])
					log.info('{}{}{}{}{}: {}\n'.format(clr.red, 'shd__ ', clr.end, 'Atypical scanning failure on ', assembly_label, str(e)))
					continue

				#######################################
				## Stage two!! Extract distributions ##
				#######################################
				log.info('{}{}{}{}'.format(clr.yellow,'shd__ ',clr.end,'Extracting repeat distributions from pre-assembled data..'))
				if self.args.purgesam:
					purged_fw = align.SeqAlign.purge_alignment_map(forward_path, forward_assembly)
					purged_rv = align.SeqAlign.purge_alignment_map(reverse_path, reverse_assembly)
					assembly_data[0] = align.SeqAlign.extract_repeat_distributions(assembly_label, forward_path, purged_fw)
					assembly_data[1] = align.SeqAlign.extract_repeat_distributions(assembly_label, reverse_path, purged_rv)
				else:
					assembly_data[0] = align.SeqAlign.extract_repeat_distributions(assembly_label,forward_path,forward_assembly)
					assembly_data[1] = align.SeqAlign.extract_repeat_distributions(assembly_label,reverse_path,reverse_assembly)
				log.info('{}{}{}{}'.format(clr.green,'shd__ ',clr.end,'Repeat distribution extraction complete!'))

				###########################################
				## Stage three!! Process allele status.. ##
				###########################################
				if atypical_count != 0:
					try:
						log.info('{}{}{}{}'.format(clr.yellow,'shd__ ',clr.end,'Currently in --batch mode, cannot realign. Trusting DSP genotyping.'))
						genotype_report = predict.DSPResultsGenerator(assembly_data, predict_path, processed_atypical).get_report()
						self.instance_graphs[assembly_label] = collate_peaks(predict_path, assembly_label)
						gc.collect()
						log.info('{}{}{}{}'.format(clr.green,'shd__ ',clr.end,'Genotyping workflow complete!'))
						workflow_dictionary = {'GenotypeReport': genotype_report, 'SampleLabel': assembly_label, 'MasterResultsFile': master_results_file, 'MasterMatrixFile': master_matrix_file}
						self.write_to_report(workflow_dictionary)
						del genotype_report  # rework when reporting overhaul
						log.info('{}{}{}{}'.format(clr.green, 'shd__ ', clr.end, 'Assembly pair workflow complete!\n'))
						continue
					except Exception, e:
						self.write_failure([assembly_label, placeholder_dict, master_results_file, master_matrix_file])
						log.info('{}{}{}{}{}: {}\n'.format(clr.red, 'shd__ ', clr.end, 'Processing atypical allele failure on ', assembly_label, str(e)))
						continue

				############################################
				## Stage four!! Distribution Genotyping.. ##
				############################################
				try:
					log.info('{}{}{}{}'.format(clr.yellow,'shd__ ',clr.end,'Executing genotyping workflow..'))
					genotype_report = predict.GenotypePrediction(assembly_data, predict_path, self.training_data, instance_params, processed_atypical).get_report()
					self.instance_graphs[assembly_label] = collate_peaks(predict_path, assembly_label)
					gc.collect()
					log.info('{}{}{}{}'.format(clr.green,'shd__ ',clr.end,'Genotyping workflow complete!'))
				except Exception, e:
					self.write_failure([assembly_label, placeholder_dict, master_results_file, master_matrix_file])
					log.info('{}{}{}{}{}: {}'.format(clr.red,'shd__ ',clr.end,'Genotyping failure on ',assembly_label,str(e)))
					continue

				########################################
				## Stage five!! Bayesian Genotyping.. ##
				########################################
				#try:
				# log.info('{}{}{}{}'.format(clr.yellow,'shd__ ',clr.end,'Experimental Bayesian workflow..'))
				# predict.BayesianLikelihood(bayesian_path, assembly_data, self.likelihood_matrix, self.raw_matrix)
				# gc.collect()
				# log.info('{}{}{}{}'.format(clr.green,'shd__ ',clr.end,'Experimental Bayesian workflow complete!'))
				# except Exception, e:
				# 	self.instance_summary[assembly_label] = {'SampleGenotype':{'BayesPrimaryAllele':'Fail',
				# 															   'BayesSecondaryAllele':'Fail',
				# 															   'BayesPredictionLikelihood':0}}
				# 	log.info('{}{}{}{}{}: {}\n'.format(clr.red,'shd__ ',clr.end,'Bayesian failure on ',assembly_label,str(e)))
				# 	continue

				##
				## Collating the required information for this data pair into a summary dictionary
				## Add dictionary to instance parent dictionary (dict of dicts for all data pairs in run...)
				workflow_dictionary = {'GenotypeReport': genotype_report, 'SampleLabel': assembly_label,
									   'MasterResultsFile':master_results_file, 'MasterMatrixFile':master_matrix_file}
				self.write_to_report(workflow_dictionary)

				##
				## Finished all desired stages for this file pair, inform user if -v
				log.info('{}{}{}{}'.format(clr.green, 'shd__ ', clr.end, 'Assembly pair workflow complete!\n'))

	def sequence_workflow(self):
		"""
		Workflow for when ScaleHD is being ran in config mode..
		Behaviours are tailored based on information extracted from the specified config XML file
		General overview:
		-- If align; index references beforehand (instead of each time we call __alignment)
		-- For each sample-pair from input:
		-- For Current pair, if quality control: run QC, modifies sequencepair_data with updated target files
			(replaces raw read files with trimmed/dmpx read files)
		-- if alignment, run alignment, modifies sequencepair_data with updated target files
			(replaces unaligned file in struct with extracted data distribution)
		-- if genotyping, pass data, returns genotyping report with relevant results within
		-- Append this sample-pair's results to report, continue with loop
		-- Process reporting
		:return: None
		"""
		##
		## Config generics
		instance_inputdata = self.instance_params.config_dict['@data_dir']

		##
		## Pre-stage: check for compressed data, extract
		if not extract_data(instance_inputdata):
			log.error('{}{}{}{}'.format(clr.red, 'shd__ ', clr.end, 'Error during file extraction. Please check data!'))

		##
		## Executing the workflow for this SHD instance
		## Ensure there are even amount of files for forward/reverse sequence pairings
		data_pairs = sequence_pairings(instance_inputdata, self.instance_rundir, 'sequence')
		for i in range(len(data_pairs)):
			for seqpair_lbl, seqpair_dat in data_pairs[i].iteritems():

				################################################
				## Pre stage! Sample object/Tree generation.. ##
				################################################
				log.info('{}{}{}{}{}/{} ({})'.format(clr.bold, 'shd__ ', clr.end, 'Processing sequence pair: ',
													 str(i + 1), str(len(data_pairs)), seqpair_lbl))
				current_seqpair = SequenceSample()
				current_seqpair.set_label(seqpair_lbl)
				current_seqpair.set_qcpath(seqpair_dat[2])
				current_seqpair.set_alignpath(seqpair_dat[3])
				current_seqpair.set_predictpath(seqpair_dat[4])
				current_seqpair.set_bayespath(seqpair_dat[5])
				current_seqpair.set_purgeflag(self.purge_flag)
				current_seqpair.set_fwidx(self.reference_indexes[0])
				current_seqpair.set_rvidx(self.reference_indexes[1])
				current_seqpair.set_fwreads(seqpair_dat[0])
				current_seqpair.set_rvreads(seqpair_dat[1])
				current_seqpair.generate_sampletree()

				############################################
				## Stage one!! Sequence quality control.. ##
				############################################
				try:
					self.quality_control(current_seqpair)
				except Exception, e:
					self.write_failure([seqpair_lbl, self.placeholder_dict, self.instance_results, self.instance_matrix])
					log.info('{}{}{}{}{}: {}\n'.format(clr.red,'shd__ ',clr.end,'SeqQC failure on ',seqpair_lbl,str(e)))
					continue

				##############################################
				## Stage two!! Sequence alignment via bwa.. ##
				##############################################
				try:
					self.sequence_alignment(current_seqpair)
				except Exception, e:
					self.write_failure([seqpair_lbl, self.placeholder_dict, self.instance_results, self.instance_matrix])
					log.info('{}{}{}{}{}: {}\n'.format(clr.red,'shd__ ',clr.end,'Alignment failure on ',seqpair_lbl,str(e)))
					continue

				###############################################
				## Stage three!! Scan for atypical alleles.. ##
				###############################################
				try:
					self.atypical_scanning(current_seqpair)
				except Exception, e:
					self.write_failure([seqpair_lbl, self.placeholder_dict, self.instance_results, self.instance_matrix])
					log.info('{}{}{}{}{}: {}\n'.format(clr.red,'shd__ ',clr.end,'Atypical scanning failure on ',seqpair_lbl,str(e)))
					continue

				##########################################
				## Stage four!! Process allele status.. ##
				##########################################
				realign_flag = self.instance_params.config_dict['instance_flags']['@atypical_realignment']
				for allele in [current_seqpair.get_primaryallele(), current_seqpair.get_secondaryallele()]:
					if allele.get_allelestatus() == 'Atypical':
						if realign_flag == 'True':
							try:
								self.sequence_realignment(current_seqpair, allele)
							except Exception, e:
								self.write_failure([seqpair_lbl, self.placeholder_dict, self.instance_results, self.instance_matrix])
								log.info('{}{}{}{}{}: {}\n'.format(clr.red,'shd__ ',clr.end,'Realignment failure on ',seqpair_lbl,str(e)))
								continue
						else:
							print '$DBG: formalise_dsp'
							##TODO dsp genotyping (when genotyping re-written)
							## try formalise_dsp_results, if flag is false, will execute
							## continue (next sequence pair in loop)
							## fail, raise formalisedsp exception (and continue with sequence loop)
					if allele.get_allelestatus() == 'Typical':
						allele.set_fwdist(current_seqpair.get_fwdist())
						allele.set_rvdist(current_seqpair.get_rvdist())
						allele.set_fwassembly(current_seqpair.get_fwassembly())
						allele.set_rvassembly(current_seqpair.get_rvassembly())

				###########################################
				## Stage five!! Genotype distributions.. ##
				###########################################
				# #try:
				# genotyping_flag = self.instance_params.config_dict['instance_flags']['@genotype_prediction']
				# if genotyping_flag == 'True':
				# 	log.info('{}{}{}{}'.format(clr.yellow,'shd__ ',clr.end,'Executing genotyping workflow..'))
				# 	genotype_report = predict.GenotypePrediction(sequencepair_data, predict_path, self.training_data, self.instance_params, processed_atypical).get_report()
				# 	self.instance_graphs[sequence_label] = collate_peaks(predict_path, sequence_label)
				# 	gc.collect()
				# 	log.info('{}{}{}{}'.format(clr.green,'shd__ ',clr.end,'Genotyping workflow complete!'))
				# #except Exception, e:
				# #	self.write_failure([sequence_label, self.placeholder_dict, master_results_file, master_matrix_file])
				# #	log.info('{}{}{}{}{}: {}\n'.format(clr.red,'shd__ ',clr.end,'Genotyping failure on ',sequence_label,str(e)))
				# #	continue

				#try:
					self.allele_genotyping(current_seqpair)
				#except Exception, e:
				#	self.write_failure([seqpair_lbl, self.placeholder_dict, self.instance_results, self.instance_matrix])
				#	log.info('{}{}{}{}{}: {}\n'.format(clr.red, 'shd__ ', clr.end, 'Genotyping failure on ',seqpair_lbl, str(e)))
				#	continue

				#######################################
				## Stage six!! Bayesian Genotyping.. ##
				#######################################

				# TODO take distro and work into R functions
				# try:
				# 	log.info('{}{}{}{}'.format(clr.yellow,'shd__',clr.end,'Experimental Bayesian workflow..'))
				# 	predict.BayesianLikelihood(bayesian_path, self.likelihood_matrix, self.raw_matrix)
				# 	gc.collect()
				# 	log.info('{}{}{}{}'.format(clr.green, 'shd__ ', clr.end, 'Experimental Bayesian workflow complete!'))
				# except Exception, e:
				# 	self.instance_summary[sequence_label] = {'SampleGenotype':{'BayesPrimaryAllele':'Fail',
				# 															   'BayesSecondaryAllele':'Fail',
				# 															   'BayesPredictionLikelihood':0}}
				# 	log.info('{}{}{}{}{}: {}\n'.format(clr.red,'shd__ ',clr.end,'Bayesian failure on ',sequence_label,str(e)))
				# 	continue


				#############################
				## Finished! File output.. ##
				#############################

				##
				## Make graph for each allele
				## render CCG graph
				## if ccg_het, render one graph per allele
				## if ccg_hom, render one graph per sample_pair
				## Take into account typical/atypical allele distribution as basis for graph
				##
				## Append workflow_dictionary to self.instance_results
				## Append workflow_matrix to self.instance_matrix
				## Append summary_graphs to self.instance_graphs

				# workflow_dictionary = {'GenotypeReport': genotype_report, 'SampleLabel': seqpair_lbl,
				# 					   'MasterResultsFile':master_results_file, 'MasterMatrixFile':master_matrix_file}
				# self.write_to_report(workflow_dictionary)
				# del genotype_report #rework when reporting overhaul
				# reference_indexes = typical_indexes # Reset reference indexes to typical format
				log.info('{}{}{}{}'.format(clr.green,'shd__ ',clr.end,'Sequence pair workflow complete!\n'))

	def quality_control(self, sequencepair_object):

		seq_qc_flag = self.instance_params.config_dict['instance_flags']['@quality_control']
		if seq_qc_flag == 'True':
			log.info('{}{}{}{}'.format(clr.yellow,'shd__ ',clr.end,'Executing sequence quality control workflow..'))
			if seq_qc.SeqQC(sequencepair_object, self.instance_params, 'validate'):
				log.info('{}{}{}{}'.format(clr.bold,'shd__ ',clr.end,'Initialising trimming..'))
				sequencepair_object.set_trimreport(seq_qc.SeqQC(sequencepair_object,self.instance_params,'trim').get_trimreport())
				gc.collect()
				log.info('{}{}{}{}'.format(clr.green,'shd__ ',clr.end,'Trimming complete!'))

	def sequence_alignment(self, sequencepair_object):

		alignment_flag = self.instance_params.config_dict['instance_flags']['@sequence_alignment']
		if alignment_flag == 'True':
			log.info('{}{}{}{}'.format(clr.yellow,'shd__ ',clr.end,'Executing alignment workflow..'))
			sequencepair_object.set_alignreport(align.SeqAlign(sequencepair_object, self.instance_params).get_alignreport())
			gc.collect()
			log.info('{}{}{}{}'.format(clr.green,'shd__ ',clr.end,'Sequence alignment workflow complete!'))

	def atypical_scanning(self, sequencepair_object):

		log.info('{}{}{}{}'.format(clr.bold,'shd__ ',clr.end,'Scanning for atypical alleles..'))
		align.ScanAtypical(sequencepair_object, self.instance_params)
		sequencepair_object.set_atypicalreport(align.ScanAtypical(sequencepair_object, self.instance_params).get_atypicalreport())
		atypical_count = sequencepair_object.get_atypicalcount()
		if atypical_count != 0:
			log.info('{}{}{}{}'.format(clr.yellow,'shd__ ',clr.end,'Scanning complete! Atypical allele(s) present.'))
		else:
			log.info('{}{}{}{}'.format(clr.green,'shd__ ',clr.end,'Scanning complete! No atypical allele(s) present.'))
		gc.collect()

	def sequence_realignment(self, sequencepair_object, individual_allele):

		log.info('{}{}{}{}'.format(clr.yellow,'shd__ ',clr.end,'User specified sequence re-alignment. Generating custom reference..'))
		atypical_xml = generate_atypical_xml(individual_allele, self.index_path)

		fwfasta = generate_reference(atypical_xml, self.index_path)
		rvfasta = generate_reference(atypical_xml, self.index_path)

		fwidx = align.ReferenceIndex(fwfasta, self.index_path).get_index_path()
		rvidx = align.ReferenceIndex(rvfasta, self.index_path).get_index_path()
		individual_allele.set_fwidx(fwidx)
		individual_allele.set_rvidx(rvidx)

		log.info('{}{}{}{}'.format(clr.yellow,'shd__ ',clr.end,'Re-aligning to custom reference..'))
		align.SeqAlign(sequencepair_object, self.instance_params, individual_allele)
		gc.collect()

		log.info('{}{}{}{}'.format(clr.green,'shd__ ',clr.end,'Sequence re-alignment workflow complete!'))

	def formalise_dsp_results(self, sequencepair_object):
			#log.info('{}{}{}{}'.format(clr.yellow, 'shd__ ', clr.end,
									#   'User specified no sequence re-alignment. Trusting DSP genotyping.'))

		# 			genotype_report = predict.DSPResultsGenerator(sequencepair_data, predict_path, processed_atypical).get_report()
		# 			self.instance_graphs[sequence_label] = collate_peaks(predict_path, sequence_label)
		# 			gc.collect()
		# 			log.info('{}{}{}{}'.format(clr.green,'shd__ ',clr.end,'Genotyping workflow complete!'))
		# 			workflow_dictionary = {'GenotypeReport': genotype_report, 'SampleLabel': sequence_label,
		# 								   'MasterResultsFile': master_results_file,
		# 								   'MasterMatrixFile': master_matrix_file}
		# 			self.write_to_report(workflow_dictionary)
		# 			del genotype_report  # rework when reporting overhaul
		# 			log.info('{}{}{}{}'.format(clr.green, 'shd__ ', clr.end, 'Sequence pair workflow complete!\n'))
		# 			continue

		atypical_flag = self.instance_params.config_dict['instance_flags']['@atypical_realignment']
		if atypical_flag == 'False':
			log.info('{}{}{}{}'.format(clr.yellow,'shd__',clr.end,'User specified no sequence re-alignment. Trusting DSP.'))
			##TODO genotype report from DSP results generation
			##TODO graphs
			##TODO workflow_dictionary, write_to_report
			log.info('{}{}{}{}'.format(clr.green,'shd__ ',clr.end,'Sequence pair workflow complete!\n'))

	def allele_genotyping(self, sequencepair_object):

		genotype_report = []
		for allele in [(sequencepair_object.get_primaryallele(), 'primary'), (sequencepair_object.get_secondaryallele(), 'secondary')]:
			log.info('{}{}{}{}{}{}'.format(clr.yellow,'shd__ ',clr.end,'Genotyping ', allele[1], ' allele..'))
			allele_report = predict.AlleleGenotyping(sequencepair_object, allele[0], self.instance_params, self.training_data).get_report()
			genotype_report.append(allele_report)
		log.info('{}{}{}{}'.format(clr.green,'shd__ ',clr.end,'Genotyping workflow complete!'))

	def write_failure(self, fail_list):

		sequence_label = fail_list[0]
		genotype_report = fail_list[1]
		master_results_file = fail_list[2]
		master_matrix_file = fail_list[3]

		self.instance_summary[sequence_label] = {
			'SampleGenotype': {'PrimaryAllele': 'Fail', 'PrimaryAlleleStatus': 'Fail',
							   'PrimaryAlleleReference': 'Fail', 'PrimaryAlleleOriginal': 'Fail',
							   'SecondaryAllele': 'Fail', 'SecondaryAlleleStatus': 'Fail',
							   'SecondaryAlleleReference': 'Fail', 'SecondaryAlleleOriginal': 'Fail',
							   'PredictionConfidence': 0, 'ForwardDistribution':[0]*4000, 'ReverseDistribution':[0]*4000}}
		workflow_dictionary = {'GenotypeReport': genotype_report, 'SampleLabel': sequence_label,
							   'MasterResultsFile': master_results_file,
							   'MasterMatrixFile': master_matrix_file}
		self.write_to_report(workflow_dictionary)

	def write_to_report(self, workflow_dictionary):

		"""
		Function to take information gathered during the specified workflow
		Turns it into (temporary) reports.. More stages added etc in future
		Eventually HTML5 based web-app report will replace this entire thing
		:param workflow_dictionary: Dict of info
		:param allele_status: typical or atypical
		"""

		genotype_report = workflow_dictionary['GenotypeReport']
		assembly_label = workflow_dictionary['SampleLabel']
		master_results_file = workflow_dictionary['MasterResultsFile']
		master_matrix_file = workflow_dictionary['MasterMatrixFile']
		datapair_summary = {'R1Trimming': 'n/a', 'R1Alignment': 'n/a',
							'R2Trimming': 'n/a', 'R2Alignment': 'n/a',
							'SampleGenotype': genotype_report}
		self.instance_summary[assembly_label] = datapair_summary

		##
		## Write the current sample's results to the temporary instance results file
		a1 = '"{}"'.format(genotype_report['PrimaryAllele'])
		a1status = genotype_report['PrimaryAlleleStatus']
		a1label = genotype_report['PrimaryAlleleReference']
		a1original = genotype_report['PrimaryAlleleOriginal']
		a2 = '"{}"'.format(genotype_report['SecondaryAllele'])
		a2status = genotype_report['SecondaryAlleleStatus']
		a2label = genotype_report['SecondaryAlleleReference']
		a2original = genotype_report['SecondaryAlleleOriginal']
		conf = genotype_report['PredictionConfidence']

		indie_row = '{},{},{},{},{},{},{},{},{},{}\n'.format(assembly_label, a1, a1status, a1label, a1original, a2, a2status, a2label, a2original, conf)
		with open(master_results_file, 'a') as outfi: outfi.write(indie_row); outfi.close()

		##
		## Write the current sample's aggregate distribution to the matrix file
		forward_dist = genotype_report['ForwardDistribution']
		reverse_dist = genotype_report['ReverseDistribution']
		aggregate_dist = [assembly_label] + [x + y for x, y in zip(forward_dist, reverse_dist)]
		aggregate_dist.append('{}{}\n'.format(genotype_report['PrimaryAllele'], genotype_report['SecondaryAllele']))
		with open(master_matrix_file, 'a') as neofi: wr = csv.writer(neofi); wr.writerow(aggregate_dist); neofi.close()

def main():
	try:
		ScaleHD()
	except KeyboardInterrupt:
		log.error('{}{}{}{}'.format(clr.red,'shd__ ',clr.end,'Fatal: Keyboard Interrupt detected. Exiting.'))
		sys.exit(2)