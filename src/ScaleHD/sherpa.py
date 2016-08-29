#/usr/bin/python
__version__ = 0.01
__author__ = 'alastair.maxwell@glasgow.ac.uk'

import os
import sys
import argparse
import pkg_resources
import logging as log
from multiprocessing import cpu_count

##
## Backend junk
from __backend import ConfigReader
from __backend import Colour as clr
from __backend import initialise_libraries
from __backend import sanitise_inputs
from __backend import extract_data
from __backend import sequence_pairings
from __backend import sanitise_outputs
from __backend import seek_target
from __backend import sanitise_trimming_output
from __backend import sanitise_alignment_output

##
## Package stages
from . import seq_qc
from .seq_qc.__quality_control import get_trimreport
from . import align
from .align.__alignment import get_alignreport
from . import predict

##
## Globals
THREADS = cpu_count()

class ScaleHD:

	def __init__(self):

		"""
		ScaleHD -- automated Huntington Disease genotyping pipeline
		~~haha fill this out~~
		"""

		##
		## Package data
		self.generic_descriptor = pkg_resources.resource_filename(__name__, 'train/long_descr.rst')
		self.collapsed_ccg_zygosity = pkg_resources.resource_filename(__name__, 'train/polyglu_zygosity.csv')
		self.training_data = {'GenericDescriptor': self.generic_descriptor, 'CollapsedCCGZygosity': self.collapsed_ccg_zygosity}

		##
		## Argument parser from CLI
		self.parser = argparse.ArgumentParser(prog='scalehd', description='ScaleHD: Automated DNA micro-satellite genotyping.')
		self.parser.add_argument('-v', '--verbose', help='Verbose output mode. Setting this flag enables verbose output. Default: off.', action='store_true')
		input_group = self.parser.add_mutually_exclusive_group(required=True)
		input_group.add_argument('-b', '--batch', help='Input batch. Folder of multiple .sam files. For FastQ processing, use -c.', nargs=1)
		input_group.add_argument('-c', '--config', help='Pipeline config. Specify a directory to your ArgumentConfig.xml file.', nargs=1)
		self.parser.add_argument('-t', '--threads', help='Thread utilisation. Typically only alters third party alignment performance. Default: system max.', type=int, choices=xrange(1, THREADS+1), default=THREADS)
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
		self.instance_rundir = sanitise_outputs(self.args.output)
		self.instance_summary = {}

		##
		## Set up config dictionary of all params.
		## if -c used, read from XML. Else, use 'defaults' in set_params().
		script_path = os.path.dirname(__file__)
		if not self.args.config:
			self.instance_params = self.set_prediction_params()
		else:
			self.configfile = self.args.config[0]
			self.instance_params = ConfigReader(script_path, self.configfile)

		##
		## Check third party libraries before continuing
		if initialise_libraries(self.instance_params):
			log.error('{}{}{}{}'.format(clr.red, 'shd__ ', clr.end, 'Detected missing library from system/$PATH. Exiting.'))
			sys.exit(2)
		else:
			log.info('{}{}{}{}'.format(clr.green, 'shd__ ', clr.end, 'Required libraries present, assuming OK!\n'))

		##
		## Depending on input mode, direct flow of functions
		## -b == multiple files, loop files to class
		## -c == config, do as config parsed flags
		if not self.args.config: self.assembly_workflow()
		else: self.sequence_workflow()

		##
		## Print all the information from this
		## whole instance of the application -- 'master summary'
		self.process_report()
		log.info('{}{}{}{}'.format(clr.green, 'shd__ ', clr.end, 'ScaleHD pipeline completed; exiting.'))

	@staticmethod
	def set_prediction_params():

		param_dict = {'quality_control':'False',
					  'sequence_alignment':'False',
					  'genotype_prediction':'True',
					  'decision_function_shape':'ovr',
					  'probability_estimate':'True',
					  'max_iteration':'-1',}
		return param_dict

	def assembly_workflow(self):

		##
		## Input path, create pairs of data in said path
		instance_inputdata = self.args.batch[0]
		assembly_pairs = sequence_pairings(instance_inputdata, self.instance_rundir, 'assembly')

		##
		## Executing the workflow for this SHD instance
		for i in range(len(assembly_pairs)):
			log.info('{}{}{}{}{}{}{}'.format(clr.bold, 'shd__ ', clr.end, 'Processing assembly pair: ', str(i+1), '/', str(len(assembly_pairs))))
			for assembly_label, assembly_data in assembly_pairs[i].iteritems():

				##
				## Required data to process
				forward_assembly = assembly_data[0]
				reverse_assembly = assembly_data[1]
				predict_path = assembly_data[2]
				instance_params = self.set_prediction_params()

				##
				## List of paths to report files which may or may not be written to
				## Used to scrape later on for instance summary
				trim_report = []
				align_report = []

				##
				## Specific paths to pass to distribution scraper
				## In instance_workflow these would've been created for alignment, but since we don't align here
				## They have to be made in this location instead
				forward_filename = forward_assembly.split('/')[-1].split('.')[0] ##absolutely_disgusting.jpg
				reverse_filename = reverse_assembly.split('/')[-1].split('.')[0] ##absolutely_disgusting.jpg
				forward_path = os.path.join(predict_path,forward_filename)
				reverse_path = os.path.join(predict_path,reverse_filename)
				if not os.path.exists(forward_path): os.makedirs(forward_path)
				if not os.path.exists(reverse_path): os.makedirs(reverse_path)

				##
				## Pre stage -- extract distributions from input sam files
				## Update assembly_data list; replacing forward/reverse assemblys with respective repeat distributions
				log.info('{}{}{}{}'.format(clr.yellow,'shd__ ',clr.end,'Extracting repeat distributions from pre-assembled data..'))
				assembly_data[0] = align.SeqAlign.extract_repeat_distributions(assembly_label,forward_path,forward_assembly)
				assembly_data[1] = align.SeqAlign.extract_repeat_distributions(assembly_label,reverse_path,reverse_assembly)
				log.info('{}{}{}{}'.format(clr.green,'shd__ ',clr.end,'Repeat distribution extraction complete!'))

				## Prediction Stage
				log.info('{}{}{}{}'.format(clr.yellow,'shd__ ',clr.end,'Executing genotyping workflow..'))
				report = predict.GenotypePrediction(assembly_data, predict_path, self.training_data, instance_params).get_report()
				log.info('{}{}{}{}'.format(clr.green,'shd__ ',clr.end,'Genotyping workflow complete!'))

				##
				## Collating the required information for this data pair into a summary dictionary
				## Add dictionary to instance parent dictionary (dict of dicts for all data pairs in run...)
				r1_trimming = ''
				r2_trimming = ''
				r1_align = ''
				r2_align = ''

				if trim_report:
					r1_trimming = self.scrape_summary_data('trim', trim_report[0])
					r2_trimming = self.scrape_summary_data('trim', trim_report[1])
				if align_report:
					r1_align = self.scrape_summary_data('align', align_report[0])
					r2_align = self.scrape_summary_data('align', align_report[1])

				##
				## Collating the required information for this data pair into a summary dictionary
				## Add dictionary to instance parent dictionary (dict of dicts for all data pairs in run...)
				datapair_summary = {'R1_Trimming':r1_trimming,
									'R1_Alignment':r1_align,
									'R2_Trimming':r2_trimming,
									'R2_Alignment':r2_align,
									'Sample_Genotype':report}
				self.instance_summary[assembly_label] = datapair_summary

				##
				## Finished all desired stages for this file pair, inform user if -v
				log.info('{}{}{}{}'.format(clr.green, 'shd__ ', clr.end, 'Assembly pair workflow complete!\n'))

	def sequence_workflow(self):

		##
		## Config generics
		instance_inputdata = self.instance_params.config_dict['@data_dir']

		##
		## Pre-stage: check for compressed data, extract
		if not extract_data(instance_inputdata):
			log.error('{}{}{}{}'.format(clr.red, 'shd__ ', clr.end, 'Error during file extraction. Please check your input data.'))

		##
		## If the user wants to align, it is more intuitive to index references (which will be used in all pairs)
		## beforehand, so as to not repeat the computation for each cycle; thus the indexes are created outside of the
		## main "workflow" so to speak :: check flag and then get reference indexes
		reference_indexes = []
		if self.instance_params.config_dict['instance_flags']['@sequence_alignment']:
			log.info('{}{}{}{}'.format(clr.bold,'shd__ ',clr.end,'Indexing reference(s) before initialising sample pair cycle..'))
			index_path = os.path.join(self.instance_rundir,'Indexes')
			if not os.path.exists(index_path): os.makedirs(index_path)

			##
			## Ref path
			forward_reference = self.instance_params.config_dict['@forward_reference']
			reverse_reference = self.instance_params.config_dict['@reverse_reference']

			##
			## Return all bt2-index indexed files for the input reference(s)
			forward_index = align.ReferenceIndex(forward_reference, index_path).getIndexPath()
			reverse_index = align.ReferenceIndex(reverse_reference, index_path).getIndexPath()
			reference_indexes = [forward_index, reverse_index]

		##
		## Executing the workflow for this SHD instance
		## Ensure there are even amount of files for forward/reverse sequence pairings
		data_pairs = sequence_pairings(instance_inputdata, self.instance_rundir, 'sequence')
		for i in range(len(data_pairs)):
			log.info('{}{}{}{}{}{}{}'.format(clr.bold, 'shd__ ', clr.end, 'Processing sequence pair: ', str(i+1), '/', str(len(data_pairs))))
			for sequence_label, sequencepair_data in data_pairs[i].iteritems():

				##
				## For the Sequence Pair dictionary we're currently in
				## create object of the desired stage paths..
				qc_path = sequencepair_data[2]
				align_path = sequencepair_data[3]
				predict_path = sequencepair_data[4]

				##
				## List of paths to report files which may or may not be written to
				## Used to scrape later on for instance summary
				trim_report = []
				align_report = []

				##
				## Stage 1: QC and subflags
				seq_qc_flag = self.instance_params.config_dict['instance_flags']['@quality_control']
				seq_qc_trim = self.instance_params.config_dict['trim_flags']['@trim_data']

				if seq_qc_flag == 'True':
					log.info('{}{}{}{}'.format(clr.yellow,'shd__ ',clr.end,'Executing sequence quality control workflow..'))
					if seq_qc.SeqQC(sequencepair_data, qc_path, 'valid', self.instance_params):

						if seq_qc_trim == 'True':
							log.info('{}{}{}{}'.format(clr.bold, 'shd__ ', clr.end, 'Initialising trimming.'))
							seq_qc.SeqQC(sequencepair_data, qc_path, 'trim', self.instance_params)
							trim_report = get_trimreport()
							log.info('{}{}{}{}'.format(clr.green, 'shd__ ', clr.end, 'Trimming complete!'))

				##
				## Stage 2: Alignment flags
				alignment_flag = self.instance_params.config_dict['instance_flags']['@sequence_alignment']
				if alignment_flag == 'True':
					log.info('{}{}{}{}'.format(clr.yellow,'shd__ ',clr.end,'Executing alignment workflow..'))
					align.SeqAlign(sequence_label, sequencepair_data, align_path, reference_indexes, self.instance_params)
					align_report = get_alignreport()
					log.info('{}{}{}{}'.format(clr.green, 'shd__ ', clr.end, 'Sequence alignment workflow complete!'))

				##
				## Stage 3: Genotyping flags
				genotyping_flag = self.instance_params.config_dict['instance_flags']['@genotype_prediction']
				if genotyping_flag == 'True':
					log.info('{}{}{}{}'.format(clr.yellow,'shd__ ',clr.end,'Executing genotyping workflow..'))
					gtype_report = predict.GenotypePrediction(sequencepair_data, predict_path, self.training_data, self.instance_params).get_report()
					log.info('{}{}{}{}'.format(clr.green,'shd__ ',clr.end,'Genotyping workflow complete!'))

				##
				## Collating the required information for this data pair into a summary dictionary
				## Add dictionary to instance parent dictionary (dict of dicts for all data pairs in run...)
				r1_trimming = ''
				r2_trimming = ''
				r1_align = ''
				r2_align = ''

				##
				## If the stage has been specified, get the report from that stage
				## If sliced, that report is a list indicating fw/rv files
				## Genotype stage not required, due to us directly cleaning results
				## as they were generated; just scrape from relevant array indice
				if seq_qc_flag == 'True':
					r1_trimming = self.scrape_summary_data('trim', trim_report[0])
					r2_trimming = self.scrape_summary_data('trim', trim_report[1])
				if alignment_flag == 'True':
					r1_align = self.scrape_summary_data('align', align_report[0])
					r2_align = self.scrape_summary_data('align', align_report[1])

				datapair_summary = {'R1_Trimming':r1_trimming,
									'R1_Alignment':r1_align,
									'R2_Trimming':r2_trimming,
									'R2_Alignment':r2_align,
									'Sample_Genotype':gtype_report}

				self.instance_summary[sequence_label] = datapair_summary

				##
				## Clear the report lists before next iteration
				## iteritems results in next sample_pair files being appended
				## so indexing of lists breaks -- hence wipe
				del trim_report[:]
				del	align_report[:]
				del gtype_report[:]

				##
				## Finished all desired stages for this file pair, inform user if -v
				log.info('{}{}{}{}'.format(clr.green, 'shd__ ', clr.end, 'Sequence pair workflow complete!\n'))

	@staticmethod
	def scrape_summary_data(stage, input_report_file):

		##
		## If the argument input_report_file is from trimming..
		if stage == 'trim':
			with open(input_report_file, 'r') as trpf:
				trim_lines = trpf.readlines()
				##
				## Determine buffer size to slice from above array
				scraping_buffer = 8
				if '-q' in trim_lines[1]:
					scraping_buffer += 1
				##
				## Get Anchor
				summary_start = 0
				for i in range(0, len(trim_lines)):
					if '== Summary ==' in trim_lines[i]:
						summary_start = i
				##
				## Slice and close
				summary_data = trim_lines[summary_start:summary_start+scraping_buffer]
				trpf.close()
			return summary_data[2:]

		##
		## If the argument input_report_file is from alignment..
		if stage == 'align':
			with open(input_report_file,'r') as alnrpf:
				align_lines = alnrpf.readlines()
				alnrpf.close()
			##
			## No ranges required, only skip first line
			return align_lines[1:]

		##
		## No need to tidy up report for genotyping
		## since we already have the data from our own objects
		if stage == 'gtype':
			pass

	def process_report(self):

		"""
		A very large and ugly method which scrapes required strings from the instance dictionary
		Within each scraped string, some backend methods further scrape specific values
		This data is then slotted into the respective location in the master summary file
		"""

		##
		## Taking summary results dictionary and polishing
		## into a suitably readable CSV master table
		## Huge ass headers string for columns
		master_summary_file = os.path.join(self.instance_rundir, 'InstanceReport.txt')
		report_headers = '{}\t{}\t\t\t\t{}\t\t\t\t{}\t\t\t\t{}\t\t\t\t{}\n'.format('Sample Name','Forward Trimming','Forward Alignment','Reverse Trimming','Reverse Alignment','Genotype')
		report_subheaders = '{}\t{}\t{}\t{}\t{}\t' \
							'{}\t{}\t{}\t{}\t{}\t' \
							'{}\t{}\t{}\t{}\t{}\t' \
							'{}\t{}\t{}\t{}\t{}\t' \
							'{}\t{}\t{}\t{}\t{}\t' \
							'{}\t{}\t{}\n'.format('','Total reads processed','Quality-trimmed','Adapter-trimmed','Total written',
												  'Aligned 0 times','Aligned 1 time', 'Aligned >1 time','Overall alignment rate',
												  'Total reads processed','Quality-trimmed','Adapter-trimmed','Total written',
												  'Aligned 0 times','Aligned 1 time', 'Aligned >1 time','Overall alignment rate',
												  'Allele one','Allele two','CCG Zygosity Disconnect','CCG Expansion Skew',
												  'CCG Peak Ambiguity', 'CCG Density Ambiguity', 'CCG Recall Warning',
												  'CCG Peak Out Of Bounds', 'CAG Recall Warning', 'CAG Consensus Warning',
												  'FP/SP Disconnect')
		##
		## Writing to the summary file all of the information gathered previously
		## Some stages missing from instance/n.a. so check & act accordingly
		with open(master_summary_file, 'w') as msfile:
			msfile.write(report_headers)
			msfile.write(report_subheaders)
			sorted_instance = iter(sorted(self.instance_summary.iteritems()))

			for key, child_dictionary in sorted_instance:

				##
				## Forward Trimming information
				forward_trimming = child_dictionary['R1_Trimming']
				fw_total_reads = sanitise_trimming_output(seek_target(forward_trimming, 'Total reads processed'), forward_trimming)
				fw_quality_trimmed = sanitise_trimming_output(seek_target(forward_trimming, 'Quality-trimmed'), forward_trimming)
				fw_adapter_trimmed = sanitise_trimming_output(seek_target(forward_trimming, 'Reads with adapters'), forward_trimming)
				fw_total_written = sanitise_trimming_output(seek_target(forward_trimming, 'Total written (filtered)'), forward_trimming)

				##
				## Forward Alignment information
				forward_alignment = child_dictionary['R1_Alignment']
				fw_zero_align = sanitise_alignment_output(seek_target(forward_alignment, 'aligned 0 times'), forward_alignment, 0)
				fw_one_align = sanitise_alignment_output(seek_target(forward_alignment, 'aligned exactly 1 time'), forward_alignment, 1)
				fw_onepl_align = sanitise_alignment_output(seek_target(forward_alignment, 'aligned >1 times'), forward_alignment, 2)
				fw_overall_align = sanitise_alignment_output(seek_target(forward_alignment, 'overall alignment rate'), forward_alignment, 3)

				##
				## Reverse Trimming information
				reverse_trimming = child_dictionary['R2_Trimming']
				rv_total_reads = sanitise_trimming_output(seek_target(reverse_trimming, 'Total reads processed'), reverse_trimming)
				rv_quality_trimmed = sanitise_trimming_output(seek_target(reverse_trimming, 'Quality-trimmed'), reverse_trimming)
				rv_adapter_trimmed = sanitise_trimming_output(seek_target(reverse_trimming, 'Reads with adapters'), reverse_trimming)
				rv_total_written = sanitise_trimming_output(seek_target(reverse_trimming, 'Total written (filtered)'), reverse_trimming)

				##
				## Reverse Alignment information
				reverse_alignment = child_dictionary['R2_Alignment']
				rv_zero_align = sanitise_alignment_output(seek_target(reverse_alignment, 'aligned 0 times'), reverse_alignment, 0)
				rv_one_align = sanitise_alignment_output(seek_target(reverse_alignment, 'aligned exactly 1 time'), reverse_alignment, 1)
				rv_onepl_align = sanitise_alignment_output(seek_target(reverse_alignment, 'aligned >1 times'), reverse_alignment, 2)
				rv_overall_align = sanitise_alignment_output(seek_target(reverse_alignment, 'overall alignment rate'), reverse_alignment, 3)

				##
				## Genotype information
				## Work in progress, so not filled out..
				genotype_results = child_dictionary['Sample_Genotype']

				##
				## Generate the string
				## NOTE: Genotype information replaced with WiP for now
				datasample_string = '{}\t' \
									'{}\t{}\t{}\t{}\t' \
									'{}\t{}\t{}\t{}\t' \
									'{}\t{}\t{}\t{}\t' \
									'{}\t{}\t{}\t{}\t' \
									'{}\t{}\t{}\t{}\t' \
									'{}\t{}\t{}\t{}\t' \
									'{}\t{}\t{}'.format(key,
													  fw_total_reads, fw_quality_trimmed, fw_adapter_trimmed, fw_total_written,
													  fw_zero_align, fw_one_align, fw_onepl_align, fw_overall_align,
													  rv_total_reads, rv_quality_trimmed, rv_adapter_trimmed, rv_total_written,
													  rv_zero_align, rv_one_align, rv_onepl_align, rv_overall_align,
													  genotype_results[0],genotype_results[1],genotype_results[2],genotype_results[3],
													  genotype_results[4],genotype_results[5],genotype_results[6],genotype_results[7],
													  genotype_results[8],genotype_results[9],genotype_results[10])

				##
				## Write line to file
				msfile.write(datasample_string + '\n')

			msfile.close()

def main():
	try:
		ScaleHD()
	except KeyboardInterrupt:
		log.error('{}{}{}{}'.format(clr.red,'shd__ ',clr.end,'Fatal: Keyboard Interrupt detected. Exiting.'))
		sys.exit(2)