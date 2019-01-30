import os
import errno

class SequenceSample:
	def __init__(self):
		self.sample_label = ''
		self.instance_path = ''
		self.sample_qcpath = ''
		self.sample_alignpath = ''
		self.sample_predictpath = ''
		self.html_path = ''
		self.enshrine_flag = False
		self.subsample_flag = False
		self.snpobservationvalue = 0
		self.snpalgorithm = ''
		self.broad_flag = False
		self.group_flag = False
		self.avoid_furthersubsample = False
		self.total_seqreads = 0
		self.fwalnpcnt = 0.0
		self.rvalnpcnt = 0.0
		self.fwalncount = 0
		self.rvalncount = 0
		self.fwalnrmvd = 0
		self.rvalnrmvd = 0

		self.forward_index = ''
		self.reverse_index = ''
		self.forward_reads = ''
		self.reverse_reads = ''
		self.forward_assembly = ''
		self.reverse_assembly = ''
		self.forward_distribution = []
		self.reverse_distribution = []
		self.forward_trimmed = ''

		self.trim_report = []
		self.fqc_report = []
		self.align_report = []
		self.atypical_report = ''
		self.genotype_report = ''
		self.snp_report = ''

		self.primary_allele = None
		self.secondary_allele = None

		self.exception_raised = ''
		self.atypical_count = 0
		self.recall_count = 0
		self.homozygous_haplotype = False
		self.neighbouringpeaks = False
		self.diminishedpeaks = False
		self.ccguncertainty = False
		self.cctuncertainty = False
		self.svm_failure = False
		self.alignmentwarning = False
		self.atypical_alignmentwarning = False
		self.atypical_ccgrewrite = False
		self.atypical_zygrewrite = False
		self.peakinspection_warning = False
		self.ccgzygstate = ''
		self.fatalreadallele = False
		self.automatic_DSPsubsample = False
		self.distribution_readcount_warning = False
		self.novel_atypical_structure = False
		self.differential_confusion = False
		self.missed_expansion = False
		self.heuristicfilter_fail = False
		self.original_fqcount = 0
		self.subsampled_fqcount = 0

	##
	## Setters
	def set_label(self, label): self.sample_label = label
	def set_instancepath(self, instance_path): self.instance_path = instance_path
	def set_qcpath(self, qcpath): self.sample_qcpath = qcpath
	def set_alignpath(self, alignpath):	self.sample_alignpath = alignpath
	def set_predictpath(self, predictpath):	self.sample_predictpath = predictpath
	def set_htmlpath(self, htmlpath): self.html_path = htmlpath
	def set_enshrineflag(self, flag): self.enshrine_flag = flag
	def set_subsampleflag(self, flag): self.subsample_flag = flag
	def set_snpobservationvalue(self, val): self.snpobservationvalue = val
	def set_snpalgorithm(self, algo): self.snpalgorithm = algo
	def set_broadflag(self, flag): self.broad_flag = flag
	def set_groupflag(self, flag): self.group_flag = flag
	def set_avoidfurthersubsample(self, flag): self.avoid_furthersubsample = flag
	def set_totalseqreads(self, count): self.total_seqreads = count
	def set_fwalnpcnt(self, pcnt): self.fwalnpcnt = pcnt
	def set_rvalnpcnt(self, pcnt): self.rvalnpcnt = pcnt
	def set_fwalncount(self, count): self.fwalncount = count
	def set_rvalncount(self, count): self.rvalncount = count
	def set_fwalnrmvd(self, count): self.fwalnrmvd = count
	def set_rvalnrmvd(self, count): self.rvalnrmvd = count

	def set_fwidx(self, idx): self.forward_index = idx
	def set_rvidx(self, idx): self.reverse_index = idx
	def set_fwreads(self, reads): self.forward_reads = reads
	def set_rvreads(self, reads): self.reverse_reads = reads
	def set_fwassembly(self, assembly): self.forward_assembly = assembly
	def set_rvassembly(self, assembly): self.reverse_assembly = assembly
	def set_fwdist(self, dist): self.forward_distribution = dist
	def set_rvdist(self, dist): self.reverse_distribution = dist
	def set_fwtrimmed(self, reads): self.forward_trimmed = reads

	def set_trimreport(self, report): self.trim_report = report
	def set_fqcreport(self, report): self.fqc_report = report
	def set_alignreport(self, report): self.align_report = report
	def set_atypicalreport(self, report): self.atypical_report = report
	def set_genotypereport(self, report): self.genotype_report = report
	def set_snpreport(self, report): self.snp_report = report

	def set_primary_allele(self, alleleobj): self.primary_allele = alleleobj
	def set_secondary_allele(self, alleleobj): self.secondary_allele = alleleobj

	def set_exceptionraised(self, stage_string): self.exception_raised = stage_string
	def set_atypical_count(self, count): self.atypical_count = count
	def set_recallcount(self, count): self.recall_count = count
	def set_homozygoushaplotype(self, state): self.homozygous_haplotype = state
	def set_neighbouringpeaks(self, state): self.neighbouringpeaks = state
	def set_diminishedpeaks(self, state): self.diminishedpeaks = state
	def set_ccguncertainty(self, state): self.ccguncertainty = state
	def set_cctuncertainty(self, state): self.cctuncertainty = state
	def set_svm_failure(self, state): self.svm_failure = state
	def set_alignmentwarning(self, state): self.alignmentwarning = state
	def set_atypical_alignmentwarning(self, state): self.atypical_alignmentwarning = state
	def set_ccgzygstate(self, state): self.ccgzygstate = state
	def set_atypical_ccgrewrite(self, state): self.atypical_ccgrewrite = state
	def set_atypical_zygrewrite(self, state): self.atypical_zygrewrite = state
	def set_peakinspection_warning(self, state): self.peakinspection_warning = state
	def set_fatalreadallele(self, state): self.fatalreadallele = state
	def set_automatic_DSPsubsample(self, state): self.automatic_DSPsubsample = state
	def set_distribution_readcount_warning(self, state): self.distribution_readcount_warning = state
	def set_novel_atypical_structure(self, state): self.novel_atypical_structure = state
	def set_differential_confusion(self, state): self.differential_confusion = state
	def set_missed_expansion(self, state): self.missed_expansion = state
	def set_heuristicfilter(self, state): self.heuristicfilter_fail = state
	def set_original_fqcount(self, count): self.original_fqcount = count
	def set_subsampled_fqcount(self, count): self.subsampled_fqcount = count

	##
	## Getters
	def get_label(self): return self.sample_label
	def get_instancepath(self): return self.instance_path
	def get_qcpath(self): return self.sample_qcpath
	def get_alignpath(self): return self.sample_alignpath
	def get_predictpath(self): return self.sample_predictpath
	def get_htmlpath(self): return self.html_path
	def get_enshrineflag(self): return self.enshrine_flag
	def get_subsampleflag(self): return self.subsample_flag
	def get_snpobservationvalue(self): return self.snpobservationvalue
	def get_snpalgorithm(self): return self.snpalgorithm
	def get_broadflag(self): return self.broad_flag
	def get_groupflag(self): return self.group_flag
	def get_avoidfurthersubsample(self): return self.avoid_furthersubsample
	def get_totalseqreads(self): return self.total_seqreads
	def get_fwalnpcnt(self): return self.fwalnpcnt
	def get_rvalnpcnt(self): return self.rvalnpcnt
	def get_fwalncount(self): return self.fwalncount
	def get_rvalncount(self): return self.rvalncount
	def get_fwalnrmvd(self): return self.fwalnrmvd
	def get_rvalnrmvd(self): return self.rvalnrmvd

	def get_fwidx(self): return self.forward_index
	def get_rvidx(self): return self.reverse_index
	def get_fwreads(self): return self.forward_reads
	def get_rvreads(self): return self.reverse_reads
	def get_fwassembly(self): return self.forward_assembly
	def get_rvassembly(self): return self.reverse_assembly
	def get_fwdist(self): return self.forward_distribution
	def get_rvdist(self): return self.reverse_distribution
	def get_fwtrimmed(self): return self.forward_trimmed

	def get_trimreport(self): return self.trim_report
	def get_fqcreport(self): return self.fqc_report
	def get_alignreport(self): return self.align_report
	def get_atypicalreport(self): return self.atypical_report
	def get_genotypereport(self): return self.genotype_report
	def get_snpreport(self): return self.snp_report

	def get_primaryallele(self): return self.primary_allele
	def get_secondaryallele(self): return self.secondary_allele

	def get_exceptionraised(self): return self.exception_raised
	def get_atypicalcount(self): return self.atypical_count
	def get_recallcount(self): return self.recall_count
	def get_homozygoushaplotype(self): return self.homozygous_haplotype
	def get_neighbouringpeaks(self): return self.neighbouringpeaks
	def get_diminishedpeaks(self): return self.diminishedpeaks
	def get_ccguncertainty(self): return self.ccguncertainty
	def get_cctuncertainty(self): return self.cctuncertainty
	def get_alignmentwarning(self): return self.alignmentwarning
	def get_atypical_alignmentwarning(self): return self.atypical_alignmentwarning
	def get_ccgzygstate(self): return self.ccgzygstate
	def get_svm_failure(self): return self.svm_failure
	def get_atypical_ccgrewrite(self): return self.atypical_ccgrewrite
	def get_atypical_zygrewrite(self): return self.atypical_zygrewrite
	def get_peakinspection_warning(self): return self.peakinspection_warning
	def get_fatalreadallele(self): return self.fatalreadallele
	def get_automatic_DSPsubsample(self): return self.automatic_DSPsubsample
	def get_distribution_readcount_warning(self): return self.distribution_readcount_warning
	def get_novel_atypical_structure(self): return self.novel_atypical_structure
	def get_differential_confusion(self): return self.differential_confusion
	def get_missed_expansion(self): return self.missed_expansion
	def get_heuristicfilter(self): return self.heuristicfilter_fail
	def get_original_fqcount(self): return self.original_fqcount
	def get_subsampled_fqcount(self): return self.subsampled_fqcount

	##
	## Functions
	def generate_sampletree(self):
		for path in [self.sample_qcpath, self.sample_alignpath, self.sample_predictpath]:
			try:
				os.makedirs(path)
			except OSError as exc:
					if exc.errno == errno.EEXIST and os.path.isdir(path):pass
					else: raise

class IndividualAllele:
	def __init__(self):
		self.header = ''
		self.validation = False
		self.allele_genotype = ''
		self.allele_confidence = 0
		self.five_prime = ''
		self.cag_value = 0
		self.caacag_value = 0
		self.ccgcca_value = 0
		self.intervening_sequence = ''
		self.ccg_value = 0
		self.rewritten_ccg = 0
		self.unrewritten_ccg = 0
		self.cct_value = 0
		self.variant_call = ''
		self.variant_score = 0
		self.three_prime = ''

		self.allele_status = ''
		self.reference_label = ''
		self.original_reference = ''
		self.total_reads = 0
		self.peak_reads = 0
		self.typical_count = 0
		self.typical_pcnt = 0.0
		self.atypical_count = 0
		self.atypical_pcnt = 0.0
		self.fwalnpcnt = 0.0
		self.rvalnpcnt = 0.0
		self.fwalncount = 0
		self.rvalncount = 0
		self.fwalnrmvd = 0
		self.rvalnrmvd = 0

		self.forward_index = ''
		self.reverse_index = ''
		self.forward_assembly = ''
		self.reverse_assembly = ''
		self.forward_distribution = []
		self.reverse_distribution = []
		self.forward_array = []
		self.reverse_array = []
		self.forward_array_original = []
		self.reverse_array_original = []
		self.gatk_file = ''
		self.freebayes_file = ''
		self.ccg_peak_threshold = 0.0
		self.cag_peak_threshold = 0.0

		self.genotype_status = False
		self.fod_ccg = []
		self.fod_cag = []
		self.ccg_valid = False
		self.cag_valid = False
		self.interp_distance = 0.0
		self.vicinity_reads = 0
		self.immediate_dropoff = []
		self.allele_report = ''
		self.allele_graphs = []
		self.confinterval = ''

		self.interpolation_warning = False
		self.nminus_warninglevel = 0
		self.nplus_warninglevel = 0
		self.somaticmosaicism = 0.0
		self.backwards_slippage = 0.0
		self.unexpected_peaks = False
		self.fod_overwrite = False
		self.slippage_overwrite = False
		self.fatalalignmentwarning = False
		self.distribution_readcount_warning = False
		self.differential_confusion = False
		self.neighbouring_candidate = False
		self.ccg_uncertain = False

	##
	## Setters
	def set_header(self, header): self.header = header
	def set_validation(self, validate): self.validation = validate
	def set_allelegenotype(self, genotype): self.allele_genotype = genotype
	def set_alleleconfidence(self, confidence): self.allele_confidence = confidence
	def set_fiveprime(self, fp): self.five_prime = fp
	def set_cagval(self, cag): self.cag_value = cag
	def set_caacagval(self, intv1): self.caacag_value = intv1
	def set_ccgccaval(self, intv2): self.ccgcca_value = intv2
	def set_intervening(self, intv): self.intervening_sequence = intv
	def set_ccgval(self, ccg): self.ccg_value = ccg
	def set_rewrittenccg(self, ccg): self.rewritten_ccg = ccg
	def set_unrewrittenccg(self, ccg): self.unrewritten_ccg = ccg
	def set_cctval(self, cct): self.cct_value = cct
	def set_variantcall(self, call): self.variant_call = call
	def set_variantscore(self, score): self.variant_score = score
	def set_threeprime(self, tp): self.three_prime = tp

	def set_allelestatus(self, status): self.allele_status = status
	def set_referencelabel(self, label): self.reference_label = label
	def set_originalreference(self, label): self.original_reference = label
	def set_totalreads(self, count): self.total_reads = count
	def set_peakreads(self, count): self.peak_reads = count
	def set_typicalreads(self, count): self.typical_count = count
	def set_typicalpcnt(self, pcnt): self.typical_pcnt = pcnt
	def set_atypicalreads(self, count): self.atypical_count = count
	def set_atypicalpcnt(self, pcnt): self.atypical_pcnt = pcnt
	def set_fwalnpcnt(self, pcnt): self.fwalnpcnt = pcnt
	def set_rvalnpcnt(self, pcnt): self.rvalnpcnt = pcnt
	def set_fwalncount(self, count): self.fwalncount = count
	def set_rvalncount(self, count): self.rvalncount = count
	def set_fwalnrmvd(self, count): self.fwalnrmvd = count
	def set_rvalnrmvd(self, count): self.rvalnrmvd = count

	def set_fwidx(self, idx): self.forward_index = idx
	def set_rvidx(self, idx): self.reverse_index = idx
	def set_fwassembly(self, assembly): self.forward_assembly = assembly
	def set_rvassembly(self, assembly): self.reverse_assembly = assembly
	def set_fwdist(self, dist): self.forward_distribution = dist
	def set_rvdist(self, dist): self.reverse_distribution = dist
	def set_fwarray(self, array): self.forward_array = array
	def set_rvarray(self, array): self.reverse_array = array
	def set_fwarray_orig(self, array): self.forward_array_original = array
	def set_rvarray_orig(self, array): self.reverse_array_original = array
	def set_gatk_file(self, infile): self.gatk_file = infile
	def set_freebayes_file(self, infile): self.freebayes_file = infile
	def set_ccgthreshold(self, threshold): self.ccg_peak_threshold = threshold
	def set_cagthreshold(self, threshold): self.cag_peak_threshold = threshold

	def set_genotypestatus(self, status): self.genotype_status = status
	def set_fodccg(self, array): self.fod_ccg = array
	def set_fodcag(self, array): self.fod_cag = array
	def set_ccgvalid(self, status): self.ccg_valid = status
	def set_cagvalid(self, status): self.cag_valid = status
	def set_interpdistance(self,distance): self.interp_distance = distance
	def set_vicinityreads(self, value): self.vicinity_reads = value
	def set_immediate_dropoff(self, value_list): self.immediate_dropoff = value_list
	def set_allelereport(self, path): self.allele_report = path
	def set_allelegraphs(self, graph): self.allele_graphs.append(graph)
	def set_alleleconfinterval(self, value): self.confinterval = value

	def raise_interpolation_warning(self, bool): self.interpolation_warning = bool
	def set_nminuswarninglevel(self, amount): self.nminus_warninglevel = amount
	def set_npluswarninglevel(self, amount): self.nplus_warninglevel = amount
	def set_somaticmosaicism(self, amount): self.somaticmosaicism = amount
	def set_backwardsslippage(self, amount): self.backwards_slippage = amount
	def set_unexpectedpeaks(self, bool): self.unexpected_peaks = bool
	def set_fodoverwrite(self, bool): self.fod_overwrite = bool
	def set_slippageoverwrite(self, state): self.slippage_overwrite = state
	def set_fatalalignmentwarning(self, state): self.fatalalignmentwarning = state
	def set_distribution_readcount_warning(self, state): self.distribution_readcount_warning = state
	def set_differential_confusion(self, state): self.differential_confusion = state
	def set_neighbouring_candidate(self, state): self.neighbouring_candidate = state
	def set_ccguncertainty(self, state): self.ccg_uncertain = state

	##
	## Getters
	def get_header(self): return self.header
	def get_validation(self): return self.validation
	def get_allelegenotype(self): return self.allele_genotype
	def get_alleleconfidence(self): return self.allele_confidence
	def get_fiveprime(self): return self.five_prime
	def get_cag(self): return self.cag_value
	def get_caacag(self): return self.caacag_value
	def get_ccgcca(self): return self.ccgcca_value
	def get_intervening(self): return self.intervening_sequence
	def get_ccg(self): return self.ccg_value
	def get_rewrittenccg(self): return self.rewritten_ccg
	def get_unrewrittenccg(self): return self.unrewritten_ccg
	def get_cct(self): return self.cct_value
	def get_variantcall(self): return self.variant_call
	def get_variantscore(self): return self.variant_score
	def get_threeprime(self): return self.three_prime

	def get_allelestatus(self): return self.allele_status
	def get_reflabel(self): return self.reference_label
	def get_originalreference(self): return self.original_reference
	def get_totalreads(self): return self.total_reads
	def get_peakreads(self): return self.peak_reads
	def get_typicalreads(self): return self.typical_count
	def get_typicalpcnt(self): return self.typical_pcnt
	def get_atypicalreads(self): return self.atypical_count
	def get_atypicalpcnt(self): return self.atypical_pcnt
	def get_fwalnpcnt(self): return self.fwalnpcnt
	def get_rvalnpcnt(self): return self.rvalnpcnt
	def get_fwalncount(self): return self.fwalncount
	def get_rvalncount(self): return self.rvalncount
	def get_fwalnrmvd(self): return self.fwalnrmvd
	def get_rvalnrmvd(self): return self.rvalnrmvd

	def get_fwidx(self): return self.forward_index
	def get_rvidx(self): return self.reverse_index
	def get_fwassembly(self): return self.forward_assembly
	def get_rvassembly(self): return self.reverse_assembly
	def get_fwdist(self): return self.forward_distribution
	def get_rvdist(self): return self.reverse_distribution
	def get_fwarray(self): return self.forward_array
	def get_rvarray(self): return self.reverse_array
	def get_fwarray_orig(self): return self.forward_array_original
	def get_rvarray_orig(self): return self.reverse_array_original
	def get_gatk_file(self): return self.gatk_file
	def get_freebayes_file(self): return self.freebayes_file
	def get_ccgthreshold(self): return self.ccg_peak_threshold
	def get_cagthreshold(self): return self.cag_peak_threshold

	def get_genotypestatus(self): return self.genotype_status
	def get_fodccg(self): return self.fod_ccg
	def get_fodcag(self): return self.fod_cag
	def get_ccgvalid(self): return self.ccg_valid
	def get_cagvalid(self): return self.cag_valid
	def get_interpdistance(self): return self.interp_distance
	def get_vicinityreads(self): return self.vicinity_reads
	def get_immediate_dropoff(self): return self.immediate_dropoff
	def get_allelereport(self): return self.allele_report
	def get_allelegraphs(self): return self.allele_graphs
	def get_alleleconfinterval(self): return self.confinterval

	def get_interpolation_warning(self): return self.interpolation_warning
	def get_nminuswarninglevel(self): return self.nminus_warninglevel
	def get_npluswarninglevel(self): return self.nplus_warninglevel
	def get_somaticmosaicism(self): return self.somaticmosaicism
	def get_backwardsslippage(self): return self.backwards_slippage
	def get_unexpectedpeaks(self): return self.unexpected_peaks
	def get_fodoverwrite(self): return self.fod_overwrite
	def get_slippageoverwrite(self): return self.slippage_overwrite
	def get_fatalalignmentwarning(self): return self.fatalalignmentwarning
	def get_distribution_readcount_warning(self): return self.distribution_readcount_warning
	def get_differential_confusion(self): return self.differential_confusion
	def get_neighbouring_candidate(self): return self.neighbouring_candidate
	def get_ccguncertainty(self): return self.ccg_uncertain
