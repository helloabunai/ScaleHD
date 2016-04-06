##
## Generals
import os
import sys
import subprocess
import logging as log

##
## Backend junk
from ..backpack import Colour as clr
from ..seq_qc.quality_control import THREADS

distribution_files = []

class SeqAlign:

	def __init__(self, args, processed_files, instance_rundir, instance_params):

		self.args = args
		self.input_data = processed_files
		self.instance_rundir = instance_rundir
		self.instance_params = instance_params
		self.alignment_errors = False
		reference_indexes = self.index_reference()
		self.execute_alignment(reference_indexes)
		self.render_distributions()

	def index_reference(self):

		##
		## Be paranoid, check existence/validity of reference.. again
		reference_file = self.instance_params.config_dict['@reference_file']
		reference_root = reference_file.split('/')[-1].split('.')[0]
		if os.path.isfile(reference_file):
			if not (reference_file.endswith('.fa') or reference_file.endswith('.fas')):
				log.critical('{}{}{}{}'.format(clr.red,'shd__ ',clr.end,'Specified reference does not exist/is not fasta.'))
				self.alignment_errors = True
		##
		## Path to store indexes for this reference
		index_dir = os.path.join(self.instance_rundir,'Indexes')
		if not os.path.exists(index_dir): os.makedirs(index_dir)
		reference_index = '{}{}{}'.format(index_dir, '/', reference_root)

		##
		## Indexing reference with bowtie2-build
		build_subprocess = subprocess.Popen(['bowtie2-build', reference_file, '-output', reference_index], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
		build_rawoutput = build_subprocess.communicate()
		build_stdout = build_rawoutput[0]
		build_subprocess.wait()

		build_report = os.path.join(index_dir,'IndexBuildReport.txt')
		report_file = open(build_report, 'w')
		report_file.write(build_stdout)
		report_file.close()

		return reference_index

	def execute_alignment(self, reference_index):

		##
		## Determine if QC was executed (naming consistency purpose)
		quality_control = self.instance_params.config_dict['instance_flags']['@quality_control']

		extension_threshold = self.instance_params.config_dict['alignment_flags']['@extension_threshold']
		seed_size = self.instance_params.config_dict['alignment_flags']['@seed_size']
		align_mismatch = self.instance_params.config_dict['alignment_flags']['@align_mismatch']
		substring_length = self.instance_params.config_dict['alignment_flags']['@substr_length']
		substring_interval_start = self.instance_params.config_dict['alignment_flags']['@substr_interval_start']
		substring_interval_end = self.instance_params.config_dict['alignment_flags']['@substr_interval_end']
		substring_interval = '{}{}{}{}'.format('S,',substring_interval_start,',',substring_interval_end)

		fqdatalen = len(self.input_data)
		for i in range(0, len(self.input_data)):

			##
			##User feedback on alignment progres.. maybe improve later
			##if you're reading this and want better feedback, you probably know 'top' exists
			log.info('{}{}{}{}{}{}{}'.format(clr.bold,'shd__ ',clr.end,'Processing FastQ sequence file.. ',str(i+1),'/',str(fqdatalen)))

			##
			## If no QC, files will not be appended with "_trimmed" or "_dempx"
			## Ensure sample_root is not going to throw an exception from invalid splits
			if quality_control == 'False':
				sample_root = self.input_data[i].split('/')[-1].split('.')[0] ##absolutely_disgusting.jpg
			else:
				sample_root = '_'.join((self.input_data[i].split('/')[-1]).split('_')[:-1]) ##absolutely_disgusting.jpg
			alignment_outdir = os.path.join(self.instance_rundir, sample_root, 'Alignment')
			if not os.path.exists(alignment_outdir): os.makedirs(alignment_outdir)

			##THREADS             :: -P             :: CPU threads to utilise
			##extension_threshold :: -D             :: give up extending after <int> failed extends in a row
			##seed_size           :: -R             :: for reads w/ repetitive seeds, try <int> sets of seeds
			##align_mismatch      :: -N             :: max # mismatches in seed alignment; 1/0
			##substring_length    :: -L             :: length of seed substring; 3 < x < 32
			##substring_interval  :: -i S,start,end :: interval between seed substrings w/r/t read length

			alignment_outfile = '{}{}{}{}'.format(alignment_outdir, '/', sample_root, '_assembly.sam')
			bowtie_process = subprocess.Popen(['bowtie2', '-p', THREADS,
											   '-x', reference_index, self.input_data[i],
											   '-D', extension_threshold,
											   '-R', seed_size,
											   '-N', align_mismatch,
											   '-L', substring_length,
											   '-i', substring_interval,
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

			extract_repeat_distributions(sample_root, alignment_outdir, alignment_outfile)
			sys.stdout.flush()

		if self.alignment_errors:
			log.error('{}{}{}{}'.format(clr.red,'shd__ ',clr.end,'Error during alignment. Exiting.'))
			sys.exit(2)

	def render_distributions(self):
		##TODO graph the distributions (individual graph per ccg contig probably)
		pass

def extract_repeat_distributions(sample_root, alignment_outdir, alignment_outfile):

	##
	## Scrapes repeat distribution from alignment
	sorted_assembly = '{}{}'.format(alignment_outdir, '/assembly_sorted')
	view_subprocess = subprocess.Popen(['samtools', 'view', '-bS', '-@', THREADS, alignment_outfile], stdout=subprocess.PIPE)
	sort_subprocess = subprocess.Popen(['samtools', 'sort', '-@', THREADS, '-', sorted_assembly], stdin=view_subprocess.stdout)
	view_subprocess.wait(); sort_subprocess.wait()

	index_subprocess = subprocess.Popen(['samtools', 'index', sorted_assembly + '.bam'], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
	index_subprocess.wait()

	raw_repeat_distribution = os.path.join(alignment_outdir, 'Raw_RepeatDistribution.txt')
	rrd_file = open(raw_repeat_distribution, 'w')
	idxstats_subprocess = subprocess.Popen(['samtools', 'idxstats', sorted_assembly + '.bam'], stdout=rrd_file)
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

	filestring = str(len(repeat_values))+',3,' + sample_root + '\n'
	filestring += data_string
	csv_path = os.path.join(alignment_outdir,'RepeatDistribution.csv')
	csv_file = open(csv_path, 'w')
	csv_file.write(filestring)
	csv_file.close()
	distribution_files.append(csv_path)
	os.remove(raw_repeat_distribution)
	##
	## We return this single csv for when the function is called from shd/prediction
	## That call loops through a -i/-b sam input file individually, doesn't need a list
	## -c input utilises the distribution_files list
	return csv_path

def get_repeat_distributions():
	return distribution_files


