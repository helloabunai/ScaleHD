import os
import errno

class SequenceSample:
	def __init__(self):
		self.sample_label = ''
		self.sample_qcpath = ''
		self.sample_alignpath = ''
		self.sample_predictpath = ''
		self.sample_bayespath = ''
		self.purge_flag = False
		self.forward_index = ''
		self.reverse_index = ''
		self.forward_reads = ''
		self.reverse_reads = ''
		self.forward_assembly = ''
		self.reverse_assembly = ''
		self.forward_distribution = []
		self.reverse_distribution = []
		self.trim_report = []
		self.align_report = []
		self.atypical_report = ''
		self.genotype_report = ''
		self.primary_allele = None
		self.secondary_allele = None
		self.atypical_count = 0

	##
	## Setters
	def set_label(self, label): self.sample_label = label
	def set_qcpath(self, qcpath): self.sample_qcpath = qcpath
	def set_alignpath(self, alignpath):	self.sample_alignpath = alignpath
	def set_predictpath(self, predictpath):	self.sample_predictpath = predictpath
	def set_bayespath(self, bayespath):	self.sample_bayespath = bayespath
	def set_purgeflag(self, flag): self.purge_flag = flag
	def set_fwidx(self, idx): self.forward_index = idx
	def set_rvidx(self, idx): self.reverse_index = idx
	def set_fwreads(self, reads): self.forward_reads = reads
	def set_rvreads(self, reads): self.reverse_reads = reads
	def set_fwassembly(self, assembly): self.forward_assembly = assembly
	def set_rvassembly(self, assembly): self.reverse_assembly = assembly
	def set_fwdist(self, dist): self.forward_distribution = dist
	def set_rvdist(self, dist): self.reverse_distribution = dist
	def set_trimreport(self, report): self.trim_report = report
	def set_alignreport(self, report): self.align_report = report
	def set_atypicalreport(self, report): self.atypical_report = report
	def set_genotypereport(self, report): self.genotype_report = report
	def set_primary_allele(self, alleleobj): self.primary_allele = alleleobj
	def set_secondary_allele(self, alleleobj): self.secondary_allele = alleleobj
	def set_atypical_count(self, count): self.atypical_count = count

	##
	## Getters
	def get_label(self): return self.sample_label
	def get_qcpath(self): return self.sample_qcpath
	def get_alignpath(self): return self.sample_alignpath
	def get_predictpath(self): return self.sample_alignpath
	def get_bayespath(self): return self.sample_bayespath
	def get_purgeflag(self): return self.purge_flag
	def get_fwidx(self): return self.forward_index
	def get_rvidx(self): return self.reverse_index
	def get_fwreads(self): return self.forward_reads
	def get_rvreads(self): return self.reverse_reads
	def get_fwassembly(self): return self.forward_assembly
	def get_rvassembly(self): return self.reverse_assembly
	def get_fwdist(self): return self.forward_distribution
	def get_rvdist(self): return self.reverse_distribution
	def get_trimreport(self): return self.trim_report
	def get_alignreport(self): return self.align_report
	def get_atypicalreport(self): return self.atypical_report
	def get_genotypereport(self): return self.genotype_report
	def get_primaryallele(self): return self.primary_allele
	def get_secondaryallele(self): return self.secondary_allele
	def get_atypicalcount(self): return self.atypical_count

	##
	## Functions
	def generate_sampletree(self):
		for path in [self.sample_qcpath, self.sample_alignpath, self.sample_predictpath, self.sample_bayespath]:
			try:
				os.makedirs(path)
			except OSError as exc:
					if exc.errno == errno.EEXIST and os.path.isdir(path):pass
					else: raise

class IndividualAllele:
	def __init__(self):
		self.header = ''
		self.five_prime = ''
		self.cag_value = 0
		self.caacag_value = 0
		self.ccgcca_value = 0
		self.intervening_sequence = ''
		self.ccg_value = 0
		self.cct_value = 0
		self.three_prime = ''

		self.allele_status = ''
		self.reference_label = ''
		self.original_reference = ''
		self.total_reads = 0
		self.typical_count = 0
		self.typical_pcnt = 0.0
		self.atypical_count = 0
		self.atypical_pcnt = 0.0

		self.forward_index = ''
		self.reverse_index = ''
		self.forward_assembly = ''
		self.reverse_assembly = ''
		self.forward_distribution = []
		self.reverse_distribution = []
		self.genotype_status = False

	##
	## Setters
	def set_header(self, header): self.header = header
	def set_fiveprime(self, fp): self.five_prime = fp
	def set_cagval(self, cag): self.cag_value = cag
	def set_caacagval(self, intv1): self.caacag_value = intv1
	def set_ccgccaval(self, intv2): self.ccgcca_value = intv2
	def set_intervening(self, intv): self.intervening_sequence = intv
	def set_ccgval(self, ccg): self.ccg_value = ccg
	def set_cctval(self, cct): self.cct_value = cct
	def set_threeprime(self, tp): self.three_prime = tp

	def set_allelestatus(self, status): self.allele_status = status
	def set_referencelabel(self, label): self.reference_label = label
	def set_originalreference(self, label): self.original_reference = label
	def set_totalreads(self, count): self.total_reads = count
	def set_typicalreads(self, count): self.typical_count = count
	def set_typicalpcnt(self, pcnt): self.typical_pcnt = pcnt
	def set_atypicalreads(self, count): self.atypical_count = count
	def set_atypicalpcnt(self, pcnt): self.atypical_pcnt = pcnt

	def set_fwidx(self, idx): self.forward_index = idx
	def set_rvidx(self, idx): self.reverse_index = idx
	def set_fwassembly(self, assembly): self.forward_assembly = assembly
	def set_rvassembly(self, assembly): self.reverse_assembly = assembly
	def set_fwdist(self, dist): self.forward_distribution = dist
	def set_rvdist(self, dist): self.reverse_distribution = dist
	def set_genotypestatus(self, status): self.genotype_status = status

	##
	## Getters
	def get_header(self): return self.header
	def get_fiveprime(self): return self.five_prime
	def get_cag(self): return self.cag_value
	def get_caacag(self): return self.caacag_value
	def get_ccgcca(self): return self.ccgcca_value
	def get_intervening(self): return self.intervening_sequence
	def get_ccg(self): return self.ccg_value
	def get_cct(self): return self.cct_value
	def get_threeprime(self): return self.three_prime

	def get_allelestatus(self): return self.allele_status
	def get_reflabel(self): return self.reference_label
	def get_originalreference(self): return self.original_reference
	def get_totalreads(self): return self.total_reads
	def get_typicalreads(self): return self.typical_count
	def get_typicalpcnt(self): return self.typical_pcnt
	def get_atypicalreads(self): return self.atypical_count
	def get_atypicalpcnt(self): return self.atypical_pcnt

	def get_fwidx(self): return self.forward_index
	def get_rvidx(self): return self.reverse_index
	def get_fwassembly(self): return self.forward_assembly
	def get_rvassembly(self): return self.reverse_assembly
	def get_fwdist(self): return self.forward_distribution
	def get_rvdist(self): return self.reverse_distribution
	def get_genotypestatus(self): return self.genotype_status