from __future__ import division

#/usr/bin/python
__version__ = 0.249
__author__ = 'alastair.maxwell@glasgow.ac.uk'

##
## Python libraries
import os
import sys
import gc
import PyPDF2
import argparse
import pkg_resources
import logging as log
import datetime as dt
from shutil import copyfile
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
from __backend import generate_atypical_xml
from __backend import generate_reference
from __allelecontainer import SequenceSample

##
## Package stages
from . import seq_qc
from . import align
from . import predict

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
		self.parser = argparse.ArgumentParser(prog='scalehd', description='ScaleHD (ALSPAC): Automated DNA micro-satellite genotyping.')
		self.parser.add_argument('-v', '--verbose', help='Verbose output mode. Setting this flag enables verbose output. Default: off.', action='store_true')
		self.parser.add_argument('-c', '--config', help='Pipeline config. Specify a directory to your ArgumentConfig.xml file.', nargs=1, required=True)
		self.parser.add_argument('-t', '--threads', help='Thread utilisation. Typically only alters third party alignment performance. Default: system max.', type=int, choices=xrange(1, THREADS+1), default=THREADS)
		self.parser.add_argument('-e', '--enshrine', help='Do not remove non-uniquely mapped reads from SAM files.', action='store_true')
		self.parser.add_argument('-b', '--broadscope', help='Do not subsample fastq data in the case of high read-count.', action='store_true')
		self.parser.add_argument('-g', '--groupsam', help='Outputs all sorted SAM files into one instance-wide output folder, rather than sample subfolders.', action='store_true')
		self.parser.add_argument('-j', '--jobname', help='Customised folder output name. If not specified, defaults to normal output naming schema.', type=str)
		self.parser.add_argument('-o', '--output', help='Output path. Specify a directory you wish output to be directed towards.', metavar='output', nargs=1, required=True)
		self.args = self.parser.parse_args()
		self.header = ''

		##
		## Set verbosity for CLI output
		if self.args.verbose:
			log.basicConfig(format='%(message)s', level=log.DEBUG)
			log.info('{}{}{}{}'.format(clr.bold, 'shd__ ', clr.end, 'ScaleHD: Automated DNA micro-satellite genotyping.'))
			log.info('{}{}{}{}'.format(clr.bold, 'shd__ ', clr.end, '! ALSPAC PheWAS branch !'))
			log.info('{}{}{}{}'.format(clr.bold, 'shd__ ', clr.end, 'alastair.maxwell@glasgow.ac.uk\n'))
		else:
			log.basicConfig(format='%(message)s')

		##
		## Check inputs, generate outputs
		if sanitise_inputs(self.args):
			log.error('{}{}{}{}'.format(clr.red, 'shd__ ', clr.end, 'Error with specified input(s) configuration. Exiting.'))
			sys.exit(2)
		try:
			self.instance_rundir = sanitise_outputs(self.args.jobname, self.args.output)
		except Exception, e:
			log.error('{}{}{}{}'.format(clr.red, 'shd__ ', clr.end, e))
			sys.exit(2)
		self.enshrine_assembly = self.args.enshrine
		self.group_flag = self.args.groupsam
		self.broad_flag = self.args.broadscope
		self.instance_summary = {}; self.instance_graphs = ''

		##
		## Set up config dictionary of all params.
		## Copy configuration file to instance output folder (for reproducability)
		script_path = os.path.dirname(__file__)
		self.configfile = self.args.config[0]
		instance_configuration = os.path.join(self.instance_rundir, 'UtilisedConfiguration.xml')
		copyfile(self.configfile, instance_configuration)
		self.instance_params = ConfigReader(script_path, self.configfile)
		self.instance_params.config_dict['JobName'] = self.args.jobname
		##
		## Check libraries for stages specified in config
		if initialise_libraries(self.instance_params):
			log.error('{}{}{}{}'.format(clr.red, 'shd__ ', clr.end, 'Detected missing library from system/$PATH. Exiting.'))
			sys.exit(2)
		else:
			log.info('{}{}{}{}'.format(clr.green, 'shd__ ', clr.end, 'Required libraries present, assuming OK!\n'))

		##
		## Set-up instance wide applicable files
		self.index_path = ''; self.reference_indexes = []; self.typical_indexes = []
		self.instance_results = ''; self.instance_matrix = ''; self.instance_graphs = ''
		self.padded_distributions = ''
		self.instance_data()

		##
		## Workflow time!
		## -c == config, do as config parsed flags
		self.sequence_workflow()

		##
		## Instance wide output results
		## A simple report file is appended after each sample pair, currently..
		## In the future, replace with HTML based web-app, generated here?
		## For now, just exit
		log.info('{}{}{}{}'.format(clr.green, 'shd__ ', clr.end, 'ScaleHD (ALSPAC) pipeline completed; exiting.'))

	def instance_data(self):

		##
		## Reference indexes
		if self.args.config:
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
		self.padded_distributions = os.path.join(self.instance_rundir, 'AlignedDistributions.csv')
		self.header = '{},{},{},{},{},{},{},{},{},{},{},{},' \
					  '{},{},{},{},{},{},{},{},{},{},{},{},{},' \
					  '{},{},{},{},{},{},{},{},{},{},{},{},{},{}\n'.format(
			'SampleName', '' ,'Primary GTYPE', 'Status', 'Map (FW)', 'Map% (FW)', 'Map (RV)', 'Map% (RV)', 'BSlippage',
			'Somatic Mosaicism', 'Intervening Sequence', 'Confidence', '', 'Secondary GTYPE', 'Status', 'Map (FW)',
			'Map% (FW)', 'Map (RV)','Map% (RV)', 'BSlippage', 'Somatic Mosaicism', 'Intervening Sequence', 'Confidence',
			'', 'Exception Raised', 'Homozygous Haplotype', 'Neighbouring Peaks', 'Diminished Peaks', 'Novel Atypical',
			'Alignment Warning', 'Atypical Alignment Warning', 'CCG Rewritten', 'CCG Zygosity Rewritten',
			'CCT Uncertainty', 'SVM Failure', 'Differential Confusion', 'Peak Inspection Warning',
			'Low Distribution Reads', 'Low Peak Reads'
		)
		padded_header = '{},{},{},{},{},N-VAL\n'.format('Filename','Allele','CCGVal','Dist',' ,'*200)
		with open(self.instance_results, 'w') as outfi: outfi.write(self.header); outfi.close()
		with open(self.padded_distributions, 'w') as padfi: padfi.write(padded_header); padfi.close()

		##
		## Instance matrix (sequence distributions)
		## >> Do we really want the matrix generation anymore? <<
		# self.instance_matrix = os.path.join(self.instance_rundir, 'InstanceMatrix.csv')
		# header = ['SampleName']
		# for i in range(0,20):
		# 	for j in range(0,200):
		# 		header.append('CAG{}CCG{}'.format(j+1,i+1))
		# with open(self.instance_matrix, 'w') as neofi: wr = csv.writer(neofi); wr.writerow(header); neofi.close()

		##
		## Instance graphs
		if not self.args.jobname: job_string = 'No Jobname specified'
		else: job_string = self.args.jobname
		date_string = dt.datetime.today().strftime("%d/%m/%Y")
		self.instance_graphs = os.path.join(self.instance_rundir, 'InstanceGraphs.pdf')
		c = canvas.Canvas(self.instance_graphs, pagesize=(500,250))
		first_string = 'ScaleHD: Automated Huntington Disease Genotyping'
		second_string = 'University of Glasgow: alastair.maxwell@glasgow.ac.uk'
		third_string = '{}{}'.format('Job Name: ', job_string)
		fourth_string = '{}{}'.format('Run Date: ', date_string)
		c.drawCentredString(250, 150, first_string)
		c.drawCentredString(250, 125, second_string)
		c.drawCentredString(250, 100, third_string)
		c.drawCentredString(250, 75, fourth_string)
		c.save()

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
				current_seqpair.set_instancepath(seqpair_dat[2])
				current_seqpair.set_qcpath(seqpair_dat[3])
				current_seqpair.set_alignpath(seqpair_dat[4])
				current_seqpair.set_predictpath(seqpair_dat[5])
				current_seqpair.set_bayespath(seqpair_dat[6])
				current_seqpair.set_enshrineflag(self.enshrine_assembly)
				current_seqpair.set_broadflag(self.broad_flag)
				current_seqpair.set_groupflag(self.group_flag)
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
					current_seqpair.set_exceptionraised('SeqQC')
					self.append_report(current_seqpair)
					log.info('{}{}{}{}{}: {}\n'.format(clr.red,'shd__ ',clr.end,'SeqQC failure on ',seqpair_lbl,str(e)))
					continue
				##############################################
				## Stage two!! Sequence alignment via bwa.. ##
				##############################################
				try:
					self.sequence_alignment(current_seqpair)
				except Exception, e:
					current_seqpair.set_exceptionraised('SeqALN')
					self.append_report(current_seqpair)
					log.info('{}{}{}{}{}: {}\n'.format(clr.red,'shd__ ',clr.end,'Alignment failure on ',seqpair_lbl,str(e)))
					continue
				###############################################
				## Stage three!! Scan for atypical alleles.. ##
				###############################################
				try:
					self.atypical_scanning(current_seqpair)
				except Exception, e:
					current_seqpair.set_exceptionraised('DSP')
					self.append_report(current_seqpair)
					log.info('{}{}{}{}{}: {}\n'.format(clr.red, 'shd__ ', clr.end, 'Atypical scanning failure on ', seqpair_lbl, str(e)))
					continue
				##########################################
				## Stage four!! Process allele status.. ##
				##########################################
				realign_flag = self.instance_params.config_dict['instance_flags']['@atypical_realignment']
				invalid_data = False
				for allele in [current_seqpair.get_primaryallele(), current_seqpair.get_secondaryallele()]:
					if allele.get_allelestatus() == 'Atypical':
						if realign_flag == 'True':
							try:
								self.sequence_realignment(current_seqpair, allele)
							except Exception, e:
								current_seqpair.set_exceptionraised('SeqRE-ALN')
								self.append_report(current_seqpair)
								log.info('{}{}{}{}{}: {}'.format(clr.red,'shd__ ',clr.end,'Realignment failure on ',seqpair_lbl,str(e)))
								continue
						else:
							log.info('{}{}{}{}'.format(clr.yellow,'shd__ ',clr.end,'Atypical realignment not selected. Brute-force genotyping on inaccurate data.'))
							invalid_data = True
							allele.set_fwdist(current_seqpair.get_fwdist())
							allele.set_rvdist(current_seqpair.get_rvdist())
							allele.set_fwassembly(current_seqpair.get_fwassembly())
							allele.set_rvassembly(current_seqpair.get_rvassembly())

					if allele.get_allelestatus() == 'Typical':
						allele.set_fwidx(current_seqpair.get_fwidx())
						allele.set_fwdist(current_seqpair.get_fwdist())
						allele.set_fwassembly(current_seqpair.get_fwassembly())
						allele.set_rvidx(current_seqpair.get_rvidx())
						allele.set_rvdist(current_seqpair.get_rvdist())
						allele.set_rvassembly(current_seqpair.get_rvassembly())

				## tidy up seq files
				for seqfi in [current_seqpair.get_fwreads(), current_seqpair.get_rvreads()]:
					os.remove(seqfi)
				#####################################################
				## Stage five!! Genotype distributions/SNP Calling ##
				#####################################################
				try:
					self.allele_genotyping(current_seqpair, invalid_data)
				except Exception, e:
					current_seqpair.set_exceptionraised('Genotype')
					self.append_report(current_seqpair)
					log.info('{}{}{}{}{}: {}\n'.format(clr.red, 'shd__ ', clr.end, 'Genotyping failure on ',seqpair_lbl, str(e)))
					continue
				#######################################
				## Stage six!! Bayesian Genotyping.. ##
				#######################################
				# try:
				# 	self.bayesian_analyses(current_seqpair)
				# except Exception, e:
				#	current_seqpair.set_exceptionraised('Bayes')
				# 	self.append_report(current_seqpair, "FAIL")
				# 	log.info('{}{}{}{}{}: {}\n'.format(clr.red, 'shd__ ', clr.end, 'Bayesian failure on ', seqpair_lbl, str(e)))
				# 	continue
				#############################
				## Finished! File output.. ##
				#############################
				try:
					self.collate_graphs(current_seqpair)
					current_seqpair.set_exceptionraised('N/A')
					self.append_report(current_seqpair)
				except Exception, e:
					current_seqpair.set_exceptionraised('Report/Graph')
					self.append_report(current_seqpair)
					log.info('{}{}{}{}{}: {}'.format(clr.red, 'shd__ ', clr.end, 'Report/Graphing failure on ', seqpair_lbl, str(e)))
				gc.collect()
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

		alignment_flag = self.instance_params.config_dict['instance_flags']['@sequence_alignment']
		if alignment_flag == 'True':
			log.info('{}{}{}{}'.format(clr.bold, 'shd__ ', clr.end, 'Scanning for atypical alleles..'))
			sequencepair_object.set_atypicalreport(align.ScanAtypical(sequencepair_object, self.instance_params).get_atypicalreport())
			atypical_count = sequencepair_object.get_atypicalcount()
			if atypical_count != 0:
				log.info('{}{}{}{}{}{}'.format(clr.yellow, 'shd__ ', clr.end, 'Scanning complete! ',str(sequencepair_object.get_atypicalcount()),' atypical allele(s) present.'))
			else:
				log.info('{}{}{}{}'.format(clr.green, 'shd__ ', clr.end,'Scanning complete! No atypical alleles present.'))
			gc.collect()

	def sequence_realignment(self, sequencepair_object, individual_allele):

		log.info('{}{}{}{}'.format(clr.yellow,'shd__ ',clr.end,'User specified sequence re-alignment. Generating custom reference..'))

		atypical_index_path = os.path.join(sequencepair_object.get_alignpath(), 'AtypicalIndexes')
		if not os.path.exists(atypical_index_path):	mkdir_p(atypical_index_path)

		fw_xml = generate_atypical_xml(sequencepair_object.get_label(), individual_allele, atypical_index_path, 'fw')
		rv_xml = generate_atypical_xml(sequencepair_object.get_label(), individual_allele, atypical_index_path, 'rv')

		fwfasta = generate_reference(fw_xml, atypical_index_path, self.reference_indexes, 'fw')
		rvfasta = generate_reference(rv_xml, atypical_index_path, self.reference_indexes, 'rv')

		fwidx = align.ReferenceIndex(fwfasta, atypical_index_path).get_index_path()
		rvidx = align.ReferenceIndex(rvfasta, atypical_index_path).get_index_path()

		individual_allele.set_fwidx(fwidx)
		individual_allele.set_rvidx(rvidx)

		log.info('{}{}{}{}'.format(clr.yellow,'shd__ ',clr.end,'Re-aligning to custom reference..'))
		align.SeqAlign(sequencepair_object, self.instance_params, individual_allele)
		gc.collect()

		log.info('{}{}{}{}'.format(clr.green,'shd__ ',clr.end,'Allele re-alignment complete!'))

	def allele_genotyping(self, sequencepair_object, invalid_data):

		genotyping_flag = self.instance_params.config_dict['instance_flags']['@genotype_prediction']
		snpcall_flag = self.instance_params.config_dict['instance_flags']['@snp_calling']

		## genotyping
		if genotyping_flag == 'True':
			log.info('{}{}{}{}'.format(clr.yellow,'shd__ ',clr.end,'Genotyping alleles.. '))
			sequencepair_object.set_genotypereport(predict.AlleleGenotyping(sequencepair_object, self.instance_params, self.training_data, atypical_logic=invalid_data, padded_target=self.padded_distributions).get_report())

		## snp calling
		# if snpcall_flag == 'True':
		# 	log.info('{}{}{}{}'.format(clr.yellow,'shd__ ',clr.end,'Calling SNPs.. '))
		# 	sequencepair_object.set_snpreport(predict.SNPCalling(sequencepair_object, self.instance_params).get_report())

		## tidy up
		gc.collect()
		log.info('{}{}{}{}'.format(clr.green,'shd__ ',clr.end,'Genotyping workflow complete!'))

	def bayesian_analyses(self, sequencepair_object):

		print self.args.jobname
		# TODO take distro and work into R functions
		print 'bayes in development'

	def collate_graphs(self, sequencepair_object):

		##
		## Paths required for merging
		sample_pdf_path = os.path.join(sequencepair_object.get_predictpath(), '{}{}'.format(sequencepair_object.get_label(),'.pdf'))
		instance_path = os.path.join(self.instance_rundir, 'InstanceGraphs.pdf')

		##
		## Merge sample summary PDF with instance-wide PDF
		merger = PyPDF2.PdfFileMerger()
		for filename in [self.instance_graphs, sample_pdf_path]:
			merger.append(PyPDF2.PdfFileReader(file(filename, 'rb')))
		merger.write(instance_path)

	def append_report(self, sequencepair_object):

		primary_allele = sequencepair_object.get_primaryallele()
		secondary_allele = sequencepair_object.get_secondaryallele()

		def call_object_scraper(input_list):
			rep_str = ''
			for obj_pair in input_list:
				seq_object = obj_pair[0]
				func_call = obj_pair[1]
				if not seq_object == 'NULL':
					try:
						func = getattr(seq_object, func_call)
						func_output = func()
					except AttributeError:
						func_output = 'FAIL'
					if func_output is None:
						func_output += 'FAIL'
				else:
					func_output = ' '
				rep_str += '{},'.format(func_output)
			return rep_str

		unparsed_info = [[sequencepair_object, 'get_label'], ['NULL', 'NULL'], [primary_allele, 'get_reflabel'],
						 [primary_allele, 'get_allelestatus'], [primary_allele, 'get_fwalncount'],
						 [primary_allele, 'get_fwalnpcnt'], [primary_allele, 'get_rvalncount'],
						 [primary_allele, 'get_rvalnpcnt'], [primary_allele, 'get_backwardsslippage'],
						 [primary_allele, 'get_somaticmosaicism'], [primary_allele, 'get_intervening'],
						 [primary_allele, 'get_alleleconfidence'], ['NULL', 'NULL'], [secondary_allele, 'get_reflabel'],
						 [secondary_allele, 'get_allelestatus'], [secondary_allele, 'get_fwalncount'],
						 [secondary_allele, 'get_fwalnpcnt'], [secondary_allele, 'get_rvalncount'],
						 [secondary_allele, 'get_rvalnpcnt'], [secondary_allele, 'get_backwardsslippage'],
						 [secondary_allele, 'get_somaticmosaicism'], [secondary_allele, 'get_intervening'],
						 [secondary_allele, 'get_alleleconfidence'], ['NULL', 'NULL'],
						 [sequencepair_object, 'get_exceptionraised'],[sequencepair_object, 'get_homozygoushaplotype'],
						 [sequencepair_object, 'get_neighbouringpeaks'], [sequencepair_object, 'get_diminishedpeaks'],
						 [sequencepair_object, 'get_novel_atypical_structure'], [sequencepair_object, 'get_alignmentwarning'],
						 [sequencepair_object, 'get_atypical_alignmentwarning'], [sequencepair_object, 'get_atypical_ccgrewrite'],
						 [sequencepair_object, 'get_atypical_zygrewrite'], [sequencepair_object, 'get_cctuncertainty'],
						 [sequencepair_object, 'get_peakinspection_warning'], [sequencepair_object, 'get_svm_failure'],
						 [sequencepair_object, 'get_differential_confusion'],
						 [sequencepair_object, 'get_distribution_readcount_warning'],
						 [sequencepair_object, 'get_fatalreadallele']]

		report_string = call_object_scraper(unparsed_info)
		report_string += '\n'

		try:
			with open(self.instance_results, 'a') as outfi:
				outfi.write(report_string)
				outfi.close()
		except IOError:
			from os.path import expanduser; home = expanduser("~")
			log.error('{}{}{}{}'.format(clr.red, 'shd__ ', clr.end, 'InstanceReport.csv resource LOCKED. Open in excel?'))
			log.info('{}{}{}{}{}'.format(clr.yellow, 'shd__ ', clr.end, 'Cannot write while locked. Writing to: ', home))
			with open(os.path.join(home, 'InstanceReport.csv'), 'w') as newoutfi:
				newoutfi.write(self.header); newoutfi.close()
			with open(os.path.join(home, 'InstanceReport.csv'), 'a') as newappfi:
				newappfi.write(report_string); newappfi.close()

def main():
	try:
		ScaleHD()
	except KeyboardInterrupt:
		log.error('{}{}{}{}'.format(clr.red,'shd__ ',clr.end,'Fatal: Keyboard Interrupt detected. Exiting.'))
		sys.exit(2)