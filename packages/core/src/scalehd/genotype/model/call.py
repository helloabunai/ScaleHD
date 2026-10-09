"""``call_genotype``: candidate genotypes in, the most likely one out."""

from __future__ import annotations

import math
from itertools import combinations_with_replacement, permutations, product

import numpy as np
from scipy.special import logsumexp

from ...counts import SampleCounts
from .candidates import _truncation_bound, candidate_alleles
from .fit import _Fit, _refine, _screen
from .kernel import _log_kernel
from .lengths import _estimate_long_allele, _local_n
from .likelihood import _Model, _n_weights, _per_row, _row_loglik
from .results import AlleleCall, Candidate, Flag, GenotypeCall, NoMoleculesError
from .settings import CallerSettings
from .table import _Table

# Rounds of fitting the best genotype's length variants, to find 'best'.
_VARIANT_ROUNDS = 3


def _distribution_metrics(
    histogram: np.ndarray, n: int
) -> tuple[float | None, float | None, float | None, float | None]:
    if n >= histogram.size or histogram[n] <= 0:
        return None, None, None, None
    peak = histogram[n]
    backward = float(histogram[max(n - 2, 0) : n].sum() / peak)
    somatic = float(histogram[n + 1 : n + 11].sum() / peak)
    shifts = np.arange(histogram.size) - n
    mass = histogram.sum()
    expansion = float((np.clip(shifts, 0, None) * histogram).sum() / mass)
    contraction = float((np.clip(shifts, None, 0) * histogram).sum() / mass)
    return backward, somatic, expansion, contraction


def _allele_calls(
    table: _Table, fit: _Fit, settings: CallerSettings
) -> tuple[AlleleCall, AlleleCall]:
    model, params = fit.model, fit.params
    total, per_allele = _row_loglik(table, model, params)
    weights = _n_weights(table, total)
    shares = [params.balance, 1 - params.balance] if model.heterozygous else [1.0]
    cag = table.value[:, 0]
    exact_cag = table.exact[:, 0]
    size = int(cag.max()) + 41 if table.size else 1
    calls = []
    for allele, stutter, share, own in zip(
        model.alleles, params.stutter, shares, per_allele, strict=True
    ):
        responsibility = _per_row(
            np.exp(math.log1p(-params.background) + math.log(share) + own - total), weights
        )
        attributed = table.weight * responsibility
        histogram = np.bincount(cag[exact_cag], weights=attributed[exact_cag], minlength=size)
        if allele.beyond_read_length:
            backward = somatic = expansion = contraction = None
            estimate = _estimate_long_allele(table, fit, allele, settings)
        else:
            backward, somatic, expansion, contraction = _distribution_metrics(
                histogram, allele.structure.cag
            )
            estimate = None
        calls.append(
            AlleleCall(
                allele,
                share,
                float(attributed.sum()),
                stutter,
                backward,
                somatic,
                expansion,
                contraction,
                estimate,
            )
        )
    if len(calls) == 1:
        half = calls[0]
        half = AlleleCall(
            half.allele,
            0.5,
            half.molecules / 2,
            half.stutter,
            half.backward_slippage,
            half.somatic_mosaicism,
            half.expansion_index,
            half.contraction_index,
            half.cag_estimate,
        )
        return half, half
    return calls[0], calls[1]


def _unexplained(table: _Table, fit: _Fit, settings: CallerSettings) -> tuple[tuple[str, int], ...]:
    total, _ = _row_loglik(table, fit.model, fit.params)
    expected = table.total * _per_row(np.exp(total), _n_weights(table, total))
    called = {a.structure for a in fit.model.alleles}
    found = []
    for i in np.flatnonzero(table.complete):
        observed = table.weight[i]
        if observed < settings.unexplained_fraction * table.total:
            continue
        if observed <= 3 * expected[i] + 10:
            continue
        structure = table.observations[i].structure()
        if structure not in called:
            found.append((structure.label, int(observed)))
    return tuple(sorted(found, key=lambda item: -item[1])[:3])


def _same_configuration(a: _Model, b: _Model, shifts: tuple[int, ...]) -> bool:
    """Same alleles as ``b`` apart from CAG differences of at most ``shifts``, one per
    allele of ``b``."""
    if len(a.alleles) != len(b.alleles):
        return False
    for x, y, shift in zip(a.alleles, b.alleles, shifts, strict=True):
        if x.beyond_read_length != y.beyond_read_length:
            return False
        if x.structure.counts[1:] != y.structure.counts[1:]:
            return False
        if abs(x.structure.cag - y.structure.cag) > shift:
            return False
    return True


def _settle(
    table: _Table, fit: _Fit, settings: CallerSettings, bound: int | None
) -> list[dict[int, float] | None]:
    """Each allele's posterior over its exact N from its own peak, or None where it doesn't fit."""
    return [_local_n(table, fit, i, settings, bound) for i in range(len(fit.model.alleles))]


def _length_variants(
    model: _Model, local_posteriors: list[dict[int, float] | None], shift: int
) -> set[_Model]:
    """``model`` with each allele whose exact N its own peak didn't settle moved by up to
    ``shift`` CAG, every combination. Alleles beyond read length stay as they are."""
    choices = []
    for allele, by_n in zip(model.alleles, local_posteriors, strict=True):
        if by_n is not None or allele.beyond_read_length:
            choices.append([allele])
            continue
        cag = allele.structure.cag
        lengths = range(max(1, cag - shift), cag + shift + 1)
        choices.append([Candidate(allele.structure.with_counts(cag=n)) for n in lengths])
    if len(choices) == 1:
        return {_Model((a,)) for (a,) in product(*choices)}
    return {_Model.of(a, b) for a, b in product(*choices)}


def call_genotype(counts: SampleCounts, settings: CallerSettings | None = None) -> GenotypeCall:
    settings = settings or CallerSettings()
    candidates = candidate_alleles(counts, settings)
    table = _Table(
        counts, settings.effective_molecules, settings.stutter_window, settings.unconfirmed_error
    )
    if table.total <= 0:
        raise NoMoleculesError("no usable molecules to genotype")

    models = sorted(
        {_Model.of(a, b) for a, b in combinations_with_replacement(candidates, 2)},
        key=lambda m: m.label,
    )
    screened = sorted(models, key=lambda m: -_screen(table, m, settings))
    keep = screened[: settings.refine]
    for wanted in (True, False):  # always refine the best of each zygosity
        best = next((m for m in screened if m.heterozygous is wanted), None)
        if best is not None and best not in keep:
            keep.append(best)

    fits = [_refine(table, model, settings) for model in keep]
    leader = max(fits, key=lambda f: f.score)
    fits = [f if f is leader else _refine(table, f.model, settings, warm=leader) for f in fits]

    # Decide each allele's exact N from its own peak. Where that can't be done (alleles
    # close together), the genotypes with that allele a little shorter or longer are
    # fitted too, so they compete with the best on the whole sample.
    # Hopefully resolves automated confidence for longer alleles where spread is an issue
    bound = _truncation_bound(counts, settings)
    best_fit = max(fits, key=lambda f: f.score)
    local_posteriors = _settle(table, best_fit, settings, bound)
    for _ in range(_VARIANT_ROUNDS):
        fitted = {f.model for f in fits}
        variants = _length_variants(best_fit.model, local_posteriors, settings.local_shift)
        new = sorted(variants - fitted, key=lambda m: m.label)
        if not new:
            break
        fits += [_refine(table, model, settings, warm=best_fit) for model in new]
        leader = max(fits, key=lambda f: f.score)
        if leader is not best_fit:
            best_fit = leader
            local_posteriors = _settle(table, best_fit, settings, bound)
    scores = np.array([f.score for f in fits])
    log_posterior = scores - logsumexp(scores)
    order = np.argsort(-log_posterior)

    # refit the genotype with each settled allele at its most likely N
    refined = list(best_fit.model.alleles)
    for i, (allele, by_n) in enumerate(zip(best_fit.model.alleles, local_posteriors, strict=True)):
        if by_n:
            refined[i] = Candidate(
                allele.structure.with_counts(cag=max(by_n, key=by_n.__getitem__))
            )
    configuration = best_fit.model
    if tuple(refined) != best_fit.model.alleles:
        best_fit = _refine(table, _Model(tuple(refined)), settings, warm=best_fit)

    # P(call) = P(this configuration of alleles) x P(each called allele peak N). Only a
    # called allele's N may differ between the fits counted as this configuration. The
    # others' N is then decided by fit fitting (english hurr), so a fit with another N for
    # one of them is an alternative genotype suggestion to be considered
    shifts = tuple(settings.local_shift if by_n else 0 for by_n in local_posteriors)
    same = np.array([_same_configuration(f.model, configuration, shifts) for f in fits])
    log_configuration = float(logsumexp(log_posterior[same]))
    log_n = [
        math.log(max(p[n.structure.cag], 1e-300))
        for p, n in zip(local_posteriors, refined, strict=True)
        if p
    ]
    log_call = log_configuration + sum(log_n)
    posterior = min(1.0, float(np.exp(log_call)))
    alternatives_found: list[tuple[str, float]] = [
        (fits[int(i)].model.label, float(np.exp(log_posterior[i]))) for i in order if not same[i]
    ]
    for i, by_n in enumerate(local_posteriors):
        if not by_n:
            continue
        for n, p in by_n.items():
            if n == refined[i].structure.cag:
                continue
            variant = list(refined)
            variant[i] = Candidate(refined[i].structure.with_counts(cag=n))
            share = float(np.exp(log_call)) / max(by_n[refined[i].structure.cag], 1e-300) * p
            alternatives_found.append((_Model(tuple(variant)).label, share))
    alternatives = tuple(sorted(alternatives_found, key=lambda a: -a[1])[:3])
    error = max(1.0 - posterior, sum(p for _, p in alternatives_found), 1e-10)
    quality = min(99.0, -10 * math.log10(error))

    alleles = _allele_calls(table, best_fit, settings)
    unexplained = _unexplained(table, best_fit, settings)
    params = best_fit.params
    flags = _flags(counts, best_fit, alleles, posterior, unexplained, settings)
    return GenotypeCall(
        alleles=alleles,
        posterior=posterior,
        quality=quality,
        alternatives=alternatives,
        flags=flags,
        molecules=int(table.total),
        background=params.background,
        ccg_slippage=params.ccg_slippage,
        misread=params.misread,
        unexplained=unexplained,
    )


def _stutter_overlap(fit: _Fit, settings: CallerSettings) -> float:
    """The largest share of one allele's peak that is the other allele's stutter.
    should be 0 unless the alleles have the same structure (not counting CAG), since only then do
    their molecules mix in the actual sequencing machines."""
    model, params = fit.model, fit.params
    if not model.heterozygous:
        return 0.0
    first, second = model.alleles
    if first.beyond_read_length or second.beyond_read_length:
        return 0.0
    if first.structure.counts[1:] != second.structure.counts[1:]:
        return 0.0
    shares = (params.balance, 1 - params.balance)
    alleles = list(zip(model.alleles, params.stutter, shares, strict=True))
    largest = 0.0
    for (peak, stutter, share), (other, other_stutter, other_share) in permutations(alleles):
        at = np.array([peak.structure.cag])
        own = share * np.exp(_log_kernel(at, peak.structure.cag, stutter, settings.stutter_window))
        spill = other_share * np.exp(
            _log_kernel(at, other.structure.cag, other_stutter, settings.stutter_window)
        )
        largest = max(largest, float(spill[0] / (own[0] + spill[0])))
    return largest


def _flags(
    counts: SampleCounts,
    fit: _Fit,
    alleles: tuple[AlleleCall, AlleleCall],
    posterior: float,
    unexplained: tuple[tuple[str, int], ...],
    settings: CallerSettings,
) -> tuple[Flag, ...]:
    flags = []
    if counts.molecules - counts.unusable - counts.dropped < settings.min_molecules:
        flags.append(Flag.LOW_DEPTH)
    if posterior < settings.min_posterior:
        flags.append(Flag.LOW_CONFIDENCE)
    first, second = (a.allele for a in alleles)
    if first == second:
        flags.append(Flag.HOMOZYGOUS)
    elif (
        not first.beyond_read_length
        and not second.beyond_read_length
        and abs(first.structure.cag - second.structure.cag) == 1
        and first.structure.counts[1:] == second.structure.counts[1:]
    ):
        flags.append(Flag.NEIGHBOURING)
    elif _stutter_overlap(fit, settings) >= settings.close_alleles_share:
        flags.append(Flag.CLOSE_ALLELES)
    if any(not a.allele.structure.is_typical for a in alleles):
        flags.append(Flag.ATYPICAL)
    if any(a.allele.beyond_read_length for a in alleles):
        flags.append(Flag.BEYOND_READ_LENGTH)
    low, high = settings.balance_range
    if fit.model.heterozygous and not low <= fit.params.balance <= high:
        flags.append(Flag.ALLELE_IMBALANCE)
    if fit.params.background > settings.max_background:
        flags.append(Flag.HIGH_BACKGROUND)
    if unexplained:
        flags.append(Flag.UNEXPLAINED_PEAK)
    if counts.molecules and counts.dropped / counts.molecules > settings.max_dropped:
        flags.append(Flag.HIGH_DISCORDANCE)
    return tuple(flags)
