"""Caller settings: priors, candidate search and flag thresholds."""

from __future__ import annotations

from dataclasses import dataclass

from ...calibration import HTT_MISEQ, StutterCurve, logit

# Default reach of the stutter kernel (below, above); see CallerSettings.stutter_window.
_WINDOW = (20, 30)


## todo: make user exposed on web ui


@dataclass(frozen=True, slots=True)
class CallerSettings:
    stutter: StutterCurve = HTT_MISEQ
    # Reads are PCR copies of a limited number of input templates, not independent
    # molecules, so the likelihood is weighted down to at most this many. Without it,
    # tiny misfits in peak shape outweigh every prior once a sample has 10^5 reads.
    effective_molecules: int | None = 3000
    # How far the stutter kernel reaches (below, above) an allele. The peak region is
    # where the information about N is. With a much wider kernel, distant shoulders and
    # junk (common in ScaleHD 1.x alignments) pulled on the tail ratios and tipped N by
    # one against a clear peak/mode. Below reaches 20 because long alleles' contractions
    # fall off slowly. At 8, the molecules further down pulled N one low from about
    # CAG 55, and the exact-N step then overshot it by one.
    stutter_window: tuple[int, int] = _WINDOW
    # Candidate alleles: the most frequent complete structures, plus CAG +/-1 of the top
    # few.
    max_candidates: int = 6
    neighbour_candidates: int = 3
    min_truncated_fraction: float = 0.03
    # Candidate genotypes that get a full fit after a quick screen.
    refine: int = 6
    # Exact N of each separated allele is chosen among N +/- local_shift using only its own
    # structure's molecules within +/- local_radius of the peak (re: _local_n).
    local_radius: int = 6
    local_shift: int = 2
    # Above this share (relative to those read in full) of molecules near the peak that
    # reads didn't finish or confirm, the peak is cut by read length and N is compared
    # on every molecule instead.
    unread_share: float = 0.05
    # Priors as (median, SD) on the logit scale. Balance is the shorter allele's share
    # of molecules, whose median in the ScaleHD 1.x training matrix is 0.48.
    balance_prior: tuple[float, float] = (0.0, 0.6)
    ccg_prior: tuple[float, float] = (logit(0.01), 1.0)
    # CCG read as any value at all. In ScaleHD 1.x alignments a few percent of reads that
    # did not span the CCG tract landed on the wrong CCG entirely. Without this term
    # they distort the other allele's stutter fit.
    ccg_misassigned_prior: tuple[float, float] = (logit(0.001), 1.5)
    misread_prior: tuple[float, float] = (logit(0.002), 1.0)
    # Share of an allele's molecules spread flat over every CAG length in its own
    # structure. Without it a long flat tail distorts the stutter fit and shifts peak/N.
    floor_prior: tuple[float, float] = (logit(0.005), 1.5)
    background_prior: tuple[float, float] = (logit(0.002), 1.5)
    # Flag thresholds.
    min_molecules: int = 500
    min_posterior: float = 0.99
    max_background: float = 0.05
    max_dropped: float = 0.15
    balance_range: tuple[float, float] = (0.2, 0.8)
    unexplained_fraction: float = 0.02
    # Close alleles = one allele's stutter makes up at least this share of the molecules at
    # the other's peak, by the fit.
    close_alleles_share: float = 0.1
    # Chance that a CAG end read too near the read's own end to confirm was influenced by
    # sequencing error/quality etc, so the tract goes on.
    # About 1 in 50 in simulated reads whose error rate rises to 2% at the end.
    # Pending change upon real data reception
    unconfirmed_error: float = 0.03
