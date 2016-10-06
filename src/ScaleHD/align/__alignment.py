#/usr/bin/python
__version__ = 0.01
__author__ = 'alastair.maxwell@glasgow.ac.uk'

##
## Generals
import os
import sys
import subprocess
import logging as log

##
## Backend junk
from ..__backend import Colour as clr
from ..__backend import replace_fqfile
from ..seq_qc.__quality_control import THREADS

ALN_REPORT = []

class SeqAlign:

	def __init__(self, sequence_label=None, sequencepair_data=None, target_output=None, reference_indexes=None, instance_params=None):

		##
		## Flag
		self.alignment_errors = False

		##
		## Instance data and workflow
		self.sample_root = sequence_label
		self.sequencepair_data = sequencepair_data
		self.target_output = target_output
		self.reference_indexes = reference_indexes
		self.instance_params = instance_params
		self.alignment_reports = []
		self.alignment_workflow()

	def alignment_workflow(self):

		##
		## Get forward reference, resultant indexes, and reads
		forward_index = self.reference_indexes[0]
		forward_reads = self.sequencepair_data[0]

		##
		## Get reverse reference, resultant indexes, and reads
		reverse_index = self.reference_indexes[1]
		reverse_reads = self.sequencepair_data[1]

		##
		## Align the two FastQ files in the pair
		forward_distribution, forward_report = self.execute_alignment(forward_index,forward_reads,'Aligning forward reads..','R1')
		reverse_distribution, reverse_report = self.execute_alignment(reverse_index,reverse_reads,'Aligning reverse reads..','R2')
		ALN_REPORT.append(forward_report); ALN_REPORT.append(reverse_report)

		##
		## Update sequence pair:: replace file that was to be aligned with the distribution resulting from that file
		## list[0] was foward read FASTQ >> will be >> forward read's repeat distribution
		## list[1] was reverse read FASTQ >> will be >> reverse read's repeat distribution
		self.sequencepair_data = replace_fqfile(self.sequencepair_data, forward_reads, forward_distribution)
		self.sequencepair_data = replace_fqfile(self.sequencepair_data, reverse_reads, reverse_distribution)

	def execute_alignment(self, reference_index, target_fqfile, feedback_string, io_index):

		##
		## So. Many. Flags.
		extension_threshold = self.instance_params.config_dict['alignment_flags']['@extension_threshold']
		seed_size = self.instance_params.config_dict['alignment_flags']['@seed_size']
		align_mismatch = self.instance_params.config_dict['alignment_flags']['@align_mismatch']
		substring_length = self.instance_params.config_dict['alignment_flags']['@substr_length']
		substring_interval_start = self.instance_params.config_dict['alignment_flags']['@substr_interval_start']
		substring_interval_end = self.instance_params.config_dict['alignment_flags']['@substr_interval_end']
		substring_interval = '{}{}{}{}'.format('S,',substring_interval_start,',',substring_interval_end)
		read_gap_open = self.instance_params.config_dict['alignment_flags']['@read_gap_open']
		read_gap_extend = self.instance_params.config_dict['alignment_flags']['@read_gap_extend']
		read_penalties = '{}{}{}'.format(read_gap_open,',',read_gap_extend)
		ref_gap_open = self.instance_params.config_dict['alignment_flags']['@ref_gap_open']
		ref_gap_extend = self.instance_params.config_dict['alignment_flags']['@ref_gap_extend']
		reference_penalties = '{}{}{}'.format(ref_gap_open,',',ref_gap_extend)
		max_mismatch_penalties = self.instance_params.config_dict['alignment_flags']['@max_mismatch_pen']
		min_mismatch_penalties = self.instance_params.config_dict['alignment_flags']['@min_mismatch_pen']
		mismatch_penalties = '{}{}{}'.format(max_mismatch_penalties,',',min_mismatch_penalties)

		##
		##User feedback on alignment progres.. maybe improve later
		##if you're reading this and want better feedback, you probably know 'top' exists
		log.info('{}{}{}{}'.format(clr.bold,'shd__ ',clr.end,feedback_string))
		alignment_outdir = os.path.join(self.target_output, self.sample_root+'_'+io_index)
		if not os.path.exists(alignment_outdir): os.makedirs(alignment_outdir)

		##THREADS             :: -P             :: CPU threads to utilise
		##extension_threshold :: -D             :: give up extending after <int> failed extends in a row
		##seed_size           :: -R             :: for reads w/ repetitive seeds, try <int> sets of seeds
		##align_mismatch      :: -N             :: max # mismatches in seed alignment; 1/0
		##substring_length    :: -L             :: length of seed substring; 3 < x < 32
		##substring_interval  :: -i S,start,end :: interval between seed substrings w/r/t read length
		##read_gap_open/ext   :: --rdg          :: Read gap opening/extension penalty
		##ref_gap_open/ext    :: --rfg          :: Reference gap opening/extension penalty
		##max/min_mismatch..  :: --mp           :: Max/min mismatch score

		alignment_outfile = '{}/{}'.format(alignment_outdir, 'assembly.sam')
		metrics_outfile = os.path.join(alignment_outdir,'performance_metrics.txt')
		bowtie_process = subprocess.Popen(['bowtie2', '-p', THREADS,
										   '-x', reference_index, target_fqfile,
										   '-D', extension_threshold,
										   '-R', seed_size,
										   '-N', align_mismatch,
										   '-L', substring_length,
										   '-i', substring_interval,
										   '--rdg', read_penalties,
										   '--rfg', reference_penalties,
										   '--mp', mismatch_penalties,
										   '--met-file', metrics_outfile,
										   '--end-to-end',
										   '-S', alignment_outfile], stdout=subprocess.PIPE, stderr=subprocess.PIPE)

		##TODO explore getting the CPU usage and updating the above user feedback string periodically...?

		bowtie_rawoutput = bowtie_process.communicate()
		bowtie_stderr = bowtie_rawoutput[1]
		bowtie_process.wait()

		##
		## There's no reason bowtie should die given how secure everything up to this point is
		## but check the output anyway, just incase.
		if '(ERR): bowtie2-align exited with value' in bowtie_stderr:
			log.critical('{}{}{}{}'.format(clr.red, 'shd__ ', clr.end,'Bowtie2 fatal error.'))
			self.alignment_errors = True

		alignment_report = os.path.join(alignment_outdir, 'AlignmentReport.txt')
		report_file = open(alignment_report, 'w')
		report_file.write(bowtie_stderr)
		report_file.close()

		csv_path = self.extract_repeat_distributions(self.sample_root, alignment_outdir, alignment_outfile)
		sys.stdout.flush()

		if self.alignment_errors:
			log.error('{}{}{}{}'.format(clr.red,'shd__ ',clr.end,'Error during alignment. Exiting.'))
			sys.exit(2)

		return csv_path, alignment_report

	@staticmethod
	def extract_repeat_distributions(sample_root, alignment_outdir, alignment_outfile):

		##
		## Scrapes repeat distribution from alignment
		sorted_assembly = '{}{}'.format(alignment_outdir, '/assembly_sorted.bam')
		view_subprocess = subprocess.Popen(['samtools', 'view', '-bS', '-@', THREADS, alignment_outfile], stdout=subprocess.PIPE)
		sort_subprocess = subprocess.Popen(['samtools', 'sort', '-@', THREADS, '-', '-o', sorted_assembly], stdin=view_subprocess.stdout)
		view_subprocess.wait(); sort_subprocess.wait()

		index_subprocess = subprocess.Popen(['samtools', 'index', sorted_assembly], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
		index_subprocess.wait()

		raw_repeat_distribution = os.path.join(alignment_outdir, 'Raw_RepeatDistribution.txt')
		rrd_file = open(raw_repeat_distribution, 'w')
		idxstats_subprocess = subprocess.Popen(['samtools', 'idxstats', sorted_assembly], stdout=rrd_file)
		idxstats_subprocess.wait()
		rrd_file.close()

		##
		## Text to CSV, clean up text distribution
		with open(raw_repeat_distribution) as text_distribution:
			repeat_values = []
			data_string = ''
			for line in text_distribution.readlines()[:-1]:
				values = line.split('\t')
				data_string += values[0] + ',' + values[1] + ',' + values[2] + ',0\n'

		filestring = sample_root + '\n'
		filestring += data_string
		csv_path = os.path.join(alignment_outdir, sample_root+'_RepeatDistribution.csv')
		csv_file = open(csv_path, 'w')
		csv_file.write(filestring)
		csv_file.close()
		os.remove(raw_repeat_distribution)
		##
		## We return this single csv for when the function is called from shd/prediction
		## That call loops through a -i/-b sam input file individually, doesn't need a list
		## -c input utilises the distribution_files list and the below getter function
		return csv_path

def get_alignreport():
	return ALN_REPORT


class ReferenceIndex:

	def __init__(self, reference_file, target_output):

		self.reference = reference_file
		self.target_output = target_output
		self.reference = self.index_reference()

	def index_reference(self):

		##
		## Be paranoid, check existence/validity of reference.. again
		reference_root = self.reference.split('/')[-1].split('.')[0]
		if os.path.isfile(self.reference):
			if not (self.reference.endswith('.fa') or self.reference.endswith('.fas')):
				log.critical('{}{}{}{}'.format(clr.red,'shd__ ',clr.end,'Specified reference does not exist/is not fasta.'))
		##
		## Path to store indexes for this reference
		reference_index = os.path.join(self.target_output, reference_root)
		if not os.path.exists(reference_index): os.makedirs(reference_index)

		##
		## Indexing reference with bowtie2-build
		output_root = os.path.join(reference_index, reference_root)
		build_subprocess = subprocess.Popen(['bowtie2-build', self.reference, '-output', output_root], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
		build_rawoutput = build_subprocess.communicate()
		build_stdout = build_rawoutput[0]
		build_subprocess.wait()

		build_report = os.path.join(reference_index, 'IndexBuildReport.txt')
		report_file = open(build_report, 'w')
		report_file.write(build_stdout)
		report_file.close()

		return output_root

	def getIndexPath(self):

		return self.reference