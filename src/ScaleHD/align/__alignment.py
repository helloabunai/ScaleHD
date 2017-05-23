#/usr/bin/python
__version__ = 0.01
__author__ = 'alastair.maxwell@glasgow.ac.uk'

##
## Generals
import os
import sys
import errno
import subprocess
import shutil
import logging as log

##
## Backend junk
from ..__backend import Colour as clr
from ..seq_qc.__quality_control import THREADS

def purge_alignment_map(alignment_outdir, alignment_outfile):
	purged_assembly = '{}{}'.format(alignment_outdir, '/assembly_unique.bam')
	purged_file = open(purged_assembly, 'w')
	view_subprocess = subprocess.Popen(['samtools', 'view', '-bq', '1', '-@', str(THREADS), alignment_outfile], stdout=purged_file)
	view_subprocess.wait()
	purged_file.close()
	os.remove(alignment_outfile)
	return purged_assembly

def extract_repeat_distributions(sample_root, alignment_outdir, alignment_outfile):

	##
	## Scrapes repeat distribution from alignment
	sorted_assembly = '{}{}'.format(alignment_outdir, '/assembly_sorted.bam')
	view_subprocess = subprocess.Popen(['samtools', 'view', '-bS', '-@', str(THREADS), alignment_outfile], stdout=subprocess.PIPE)
	sort_subprocess = subprocess.Popen(['samtools', 'sort', '-@', str(THREADS), '-', '-o', sorted_assembly], stdin=view_subprocess.stdout)
	view_subprocess.wait(); sort_subprocess.wait()

	index_subprocess = subprocess.Popen(['samtools', 'index', sorted_assembly], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
	index_subprocess.wait()

	raw_repeat_distribution = os.path.join(alignment_outdir, 'RawRepeatDistribution.txt')
	rrd_file = open(raw_repeat_distribution, 'w')
	idxstats_subprocess = subprocess.Popen(['samtools', 'idxstats', sorted_assembly], stdout=rrd_file)
	idxstats_subprocess.wait()
	rrd_file.close()

	##
	## Text to CSV, clean up text distribution
	with open(raw_repeat_distribution) as text_distribution:
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
	os.remove(alignment_outfile)
	return csv_path, sorted_assembly

class SeqAlign:

	def __init__(self, sequencepair_object = None, instance_params=None, individual_allele=None):

		##
		## Instance data and workflow
		self.sequencepair_object = sequencepair_object
		self.individual_allele = individual_allele
		self.sample_root = sequencepair_object.get_label()
		self.target_output = sequencepair_object.get_alignpath()
		self.alignment_suffix = ''
		if individual_allele is not None:
			self.reference_indexes = [individual_allele.get_fwidx(), individual_allele.get_rvidx()]
		else:
			self.reference_indexes = [sequencepair_object.get_fwidx(), sequencepair_object.get_rvidx()]
		self.instance_params = instance_params
		self.purge_flag = sequencepair_object.get_purgeflag()
		self.subsample_flag = sequencepair_object.get_subsampleflag()
		self.align_report = []
		self.alignment_workflow()

	def subsample_input(self, target_file, suffix):

		if self.individual_allele is None:
			target_sample = '{}_SUB_{}.fastq'.format(self.sequencepair_object.get_label(), suffix)
			target_output = os.path.join(self.sequencepair_object.get_alignpath(),target_sample)
			target_outfi = open(target_output, 'w')
			seqtk_process = subprocess.Popen(['seqtk', 'sample', '-s100', target_file, str(self.subsample_flag)], stdout=target_outfi)
			seqtk_process.wait(); target_outfi.close()
			os.remove(target_file)
			return target_output
		else:
			return target_file

	def alignment_workflow(self):

		##
		## Get forward/reverse references/indexes
		forward_index = self.reference_indexes[0]
		reverse_index = self.reference_indexes[1]
		forward_reads = ''
		reverse_reads = ''

		##
		## Subsample input files if requested by the user
		if self.subsample_flag:
			forward_reads = self.subsample_input(self.sequencepair_object.get_fwreads(), 'R1')
			reverse_reads = self.subsample_input(self.sequencepair_object.get_rvreads(), 'R2')
			self.sequencepair_object.set_fwreads(forward_reads)
			self.sequencepair_object.set_rvreads(reverse_reads)
		if not self.subsample_flag:
			forward_reads = self.sequencepair_object.get_fwreads()
			reverse_reads = self.sequencepair_object.get_rvreads()
			self.sequencepair_object.set_fwreads(forward_reads)
			self.sequencepair_object.set_rvreads(reverse_reads)

		##
		## Align the two FastQ files in the pair
		if self.individual_allele is not None: typical_flag = 'atypical'
		else: typical_flag = 'typical'
		forward_distribution, forward_report, forward_assembly, fwmapped_pcnt = self.execute_alignment(forward_index,forward_reads,'Aligning forward reads..','R1',typical_flag)
		reverse_distribution, reverse_report, reverse_assembly, rvmapped_pcnt = self.execute_alignment(reverse_index,reverse_reads,'Aligning reverse reads..','R2',typical_flag)
		self.align_report.append(forward_report); self.align_report.append(reverse_report)

		##
		## If the user requested to group alignment output.. do so..
		if self.sequencepair_object.get_groupflag():
			try:
				os.makedirs(self.sequencepair_object.get_instancepath())
			except OSError as exc:
					if exc.errno == errno.EEXIST and os.path.isdir(self.sequencepair_object.get_instancepath()):pass
					else: raise

			forward_samfi = self.alignment_suffix+'_R1.bam'
			forward_idxfi = self.alignment_suffix+'_R1.bam.bai'
			reverse_samfi = self.alignment_suffix+'_R2.bam'
			reverse_idxfi = self.alignment_suffix+'_R2.bam.bai'

			fw_assem = os.path.join(self.sequencepair_object.get_instancepath(), forward_samfi)
			target_fwidxfi = os.path.join(self.sequencepair_object.get_instancepath(), forward_idxfi)
			rv_assem = os.path.join(self.sequencepair_object.get_instancepath(), reverse_samfi)
			target_rvidxfi = os.path.join(self.sequencepair_object.get_instancepath(), reverse_idxfi)

			os.rename(forward_assembly, fw_assem)
			os.rename(forward_assembly+'.bai', target_fwidxfi)
			os.rename(reverse_assembly, rv_assem)
			os.rename(reverse_assembly+'.bai', target_rvidxfi)

			forward_assembly = fw_assem
			reverse_assembly = rv_assem

		##
		## Update object parameters with appropriate distribution/assembly
		if not self.individual_allele:
			self.sequencepair_object.set_fwdist(forward_distribution)
			self.sequencepair_object.set_rvdist(reverse_distribution)
			self.sequencepair_object.set_fwassembly(forward_assembly)
			self.sequencepair_object.set_rvassembly(reverse_assembly)
			self.sequencepair_object.set_fwalnpcnt(fwmapped_pcnt)
			self.sequencepair_object.set_rvalnpcnt(rvmapped_pcnt)
		else:
			self.individual_allele.set_fwdist(forward_distribution)
			self.individual_allele.set_rvdist(reverse_distribution)
			self.individual_allele.set_fwassembly(forward_assembly)
			self.individual_allele.set_rvassembly(reverse_assembly)
			self.individual_allele.set_fwalnpcnt(fwmapped_pcnt)
			self.individual_allele.set_rvalnpcnt(rvmapped_pcnt)

	def execute_alignment(self, reference_index, target_fqfile, feedback_string, io_index, typical_flag):

		##
		## So. Many. Flags.
		min_seed_length = self.instance_params.config_dict['alignment_flags']['@min_seed_length']
		band_width = self.instance_params.config_dict['alignment_flags']['@band_width']
		seed_length_extension = self.instance_params.config_dict['alignment_flags']['@seed_length_extension']
		skip_seed_with_occurrence = self.instance_params.config_dict['alignment_flags']['@skip_seed_with_occurrence']
		chain_drop = self.instance_params.config_dict['alignment_flags']['@chain_drop']
		seeded_chain_drop = self.instance_params.config_dict['alignment_flags']['@seeded_chain_drop']
		seq_match_score = self.instance_params.config_dict['alignment_flags']['@seq_match_score']
		mismatch_penalty = self.instance_params.config_dict['alignment_flags']['@mismatch_penalty']
		indel_penalty = self.instance_params.config_dict['alignment_flags']['@indel_penalty']
		gap_extend_penalty = self.instance_params.config_dict['alignment_flags']['@gap_extend_penalty']
		prime_clipping_penalty = self.instance_params.config_dict['alignment_flags']['@prime_clipping_penalty']
		unpaired_pairing_penalty = self.instance_params.config_dict['alignment_flags']['@unpaired_pairing_penalty']

		##
		##User feedback on alignment progres.. maybe improve later
		##if you're reading this and want better feedback, you probably know 'htop' exists
		log.info('{}{}{}{}'.format(clr.bold,'shd__ ',clr.end,feedback_string))
		sample_string = '{}_{}_{}'.format(self.sample_root, io_index, typical_flag)
		alignment_outdir = os.path.join(self.target_output, sample_string)
		if os.path.exists(alignment_outdir):
			alignment_outdir = os.path.join(self.target_output, '{}_{}'.format(sample_string, 'alternate'))
			os.makedirs(alignment_outdir)
		if not os.path.exists(alignment_outdir):
			os.makedirs(alignment_outdir)
		self.alignment_suffix = alignment_outdir.split('/')[-1]
		aln_outpath = '{}/{}'.format(alignment_outdir, 'assembly.sam')
		aln_outfi = open(aln_outpath, 'w')

		"""
		THREADS                     :: -t <INT>      :: CPU threads to utilise [1]
		min_seed_length             :: -k <INT>      :: minimum seed length [19]
		band_width                  :: -w <INT>      :: band width for banded alignment [100]
		seed_length_extension       :: -r <FLOAT>    :: look for internal seeds inside a seed longer than <val> [1.5]
		skip_seed_with_occurrence   :: -c <INT>      :: skip seeds with more than <val> occurrences [500]
		chain_drop                  :: -D <FLOAT>    :: drop chains shorter than <val> fraction of the overlapping chain [0.50]
		seeded_chain_drop           :: -W <INT>      :: discard chain if seeded bases shorter than <val>
		seq_match_score             :: -A <INT>      :: score for sequence match [1]
		mismatch_penalty            :: -B <INT>      :: penalty for mismatch [4]
		indel_penalty               :: -O [INT, INT] :: gap open penalites for ins/del [6,6]
		gap_extend_penalty          :: -E [INT, INT] :: penalty for extending gaps [1,1]
		prime_clipping_penalty      :: -L [INT, INT] :: 5' & 3' clipping penalty [5,5]
		unpaired_pairing_penalty    :: -U <INT>      :: penalty for unpaired read pair [17]
		"""

		bwa_process = subprocess.Popen(['bwa', 'mem', '-t', str(THREADS), '-k', min_seed_length,
										'-w', band_width, '-r', seed_length_extension,
										'-c', skip_seed_with_occurrence, '-D', chain_drop, '-W', seeded_chain_drop,
										'-A', seq_match_score, '-B', mismatch_penalty, '-O', indel_penalty,
										'-E', gap_extend_penalty, '-L', prime_clipping_penalty,
										'-U', unpaired_pairing_penalty, reference_index, target_fqfile],
									    stdout=aln_outfi, stderr=subprocess.PIPE)
		bwa_error = bwa_process.communicate()[1]
		if 'illegal' in bwa_error: raise Exception('Illegal BWA behaviour: {}'.format(bwa_error))
		bwa_process.wait()
		aln_outfi.close()

		##
		## Generate an alignment report (i.e. console output to file)
		alignment_report = os.path.join(alignment_outdir, 'AlignmentReport.txt')
		report_file = open(alignment_report, 'w')
		report_file.write(bwa_error)
		report_file.close()

		##
		## If the user wants to purge reads which are not uniquely mapped
		## Then we execute that here..
		## the respective files are sent for read count extraction as normal
		if self.purge_flag:
			purged_sam = purge_alignment_map(alignment_outdir, aln_outpath)
			csv_path, sorted_assembly = extract_repeat_distributions(self.sample_root, alignment_outdir, purged_sam)
			sys.stdout.flush()
		else:
			csv_path, sorted_assembly = extract_repeat_distributions(self.sample_root, alignment_outdir, aln_outpath)
			sys.stdout.flush()

		##
		## Run samtools flagstat on alignment file
		## Set allele object's flagstat file variable..
		flagstat_path = '{}/{}'.format(alignment_outdir, 'AlignmentStats.txt')
		#flagstat_file = open(flagstat_path, 'w')
		flagstat_process = subprocess.Popen(['samtools', 'flagstat', sorted_assembly],
											stdout=subprocess.PIPE, stderr=subprocess.PIPE)
		flagstat_output = flagstat_process.communicate(); flagstat_process.wait()

		##
		## Write output to file..
		## Determine % mapped for this assembly
		with open(flagstat_path, 'w') as outfi:
			outfi.write(flagstat_output[0])
			outfi.close()
		mapped_pcnt = [x for x in (flagstat_output[0].split('\n')) if '%' in x]
		aln_pcnt = str(mapped_pcnt[0]).split('(')[1].rsplit('%')[0]

		return csv_path, alignment_report, sorted_assembly, aln_pcnt

	def get_alignreport(self):
		return self.align_report

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
			if not (self.reference.endswith('.fa') or self.reference.endswith('.fas') or self.reference.endswith('.fasta')):
				log.critical('{}{}{}{}'.format(clr.red,'shd__ ',clr.end,'Specified reference does not exist/is not fasta.'))
		##
		## Path to store indexes for this reference
		reference_index = os.path.join(self.target_output, reference_root)
		index_copy = os.path.join(reference_index, self.reference.split('/')[-1])
		if not os.path.exists(reference_index): os.makedirs(reference_index)
		shutil.copy(self.reference, os.path.join(reference_index, self.reference.split('/')[-1]))

		##
		## Indexing reference with bowtie2-build
		build_subprocess = subprocess.Popen(['bwa', 'index', index_copy], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
		build_rawoutput = build_subprocess.communicate()
		build_stderr = build_rawoutput[1]
		build_subprocess.wait()

		build_report = os.path.join(reference_index, 'IndexBuildReport.txt')
		report_file = open(build_report, 'w')
		report_file.write(build_stderr)
		report_file.close()

		return index_copy

	def get_index_path(self):

		return self.reference