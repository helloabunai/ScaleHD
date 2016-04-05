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

class SeqAlign:

	def __init__(self, processed_files, instance_rundir, instance_params):

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

		extension_threshold = self.instance_params.config_dict['alignment_flags']['@extension_threshold']
		seed_size = self.instance_params.config_dict['alignment_flags']['@seed_size']
		align_mismatch = self.instance_params.config_dict['alignment_flags']['@align_mismatch']
		substring_length = self.instance_params.config_dict['alignment_flags']['@substr_length']
		substring_interval_start = self.instance_params.config_dict['alignment_flags']['@substr_interval_start']
		substring_interval_end = self.instance_params.config_dict['alignment_flags']['@substr_interval_end']
		substring_interval = '{}{}{}{}'.format('S,',substring_interval_start,',',substring_interval_end)

		for fqfile in self.input_data:
			sample_root = '_'.join((fqfile.split('/')[-1]).split('_')[:-1]) ##absolutely_disgusting.jpg
			alignment_outdir = os.path.join(self.instance_rundir, sample_root, 'Alignment')
			if not os.path.exists(alignment_outdir): os.makedirs(alignment_outdir)

			##threads = -P
			##extension = -D
			##seedsize = -R
			##alignmismatch = -N
			##substringlength = -L
			##interval = -i S,start,end

			alignment_outfile = '{}{}{}{}'.format(alignment_outdir, '/', sample_root, '_assembly.sam')
			bowtie_process = subprocess.Popen(['bowtie2', '-p', THREADS,
											   '-x', reference_index, fqfile,
											   '-D', extension_threshold,
											   '-R', seed_size,
											   '-N', align_mismatch,
											   '-L', substring_length,
											   '-i', substring_interval,
											   '-S', alignment_outfile], stdout=subprocess.PIPE, stderr=subprocess.PIPE)

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

			self.extract_repeat_distributions(alignment_outdir, alignment_outfile)

		if self.alignment_errors:
			log.error('{}{}{}{}'.format(clr.red,'shd__ ',clr.end,'Error during alignment. Exiting.'))
			sys.exit(2)

	@staticmethod
	def extract_repeat_distributions(alignment_outdir, alignment_outfile):

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

		##TODO cleanup results file here, delete raw_RepeatDistribution.txt when cleaned into csv, csv binned into CCG contig columns

	def render_distributions(self):
		##TODO graph the distributions (individual graph per ccg probably)
		pass
