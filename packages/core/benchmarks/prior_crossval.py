"""Does it bias the legacy matrix score that the priors come from the same matrix?

The caller's stutter priors (``scalehd.calibration.HTT_MISEQ``) are per-length medians
of stutter ratios in the ScaleHD 1.x training matrix, and ``legacy_matrix.py`` scores
the caller on that same matrix. This is a bad smell but again, limited data at time of writing.
Here the curve is rebuilt by the same recipe (``curve_from``) from a random half of the samples,
and every sample is genotyped twice:

  - with the curve from the half it isn't in (that curve never saw it)
  - with the curve from its own half

If the two scores are close, the in-sample calibration isn't what 'defines' the score.

Each matrix sample takes about 20 s to call, so one split (the default) is about 1,200
calls, about 16 minutes on 24 cores. Run with
``uv run python packages/core/benchmarks/prior_crossval.py``.
"""

from __future__ import annotations

import argparse
from collections.abc import Sequence
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import numpy as np
from legacy_matrix import DEFAULT_MATRIX, Result, _match, category, genotype, load
from scalehd.calibration import HTT_MISEQ, StutterCurve
from scalehd.genotype import CallerSettings

type MatrixSample = tuple[int, np.ndarray, tuple[tuple[int, int], tuple[int, int]]]

# The recipe HTT_MISEQ was made with. It reproduces the curve to within 0.03 on
# each scale, except where noted.
# Alleles nearer than this to a partner with the same CCG share their stutter with it.
_PARTNER = 8
# Alleles with fewer reads at N give noisy ratios.
_MIN_READS = 500
# Length bins (CAG from, to). each point of the curve is at its bin's mean CAG.
_BINS = ((0, 20), (20, 30), (30, 40), (40, 50), (50, 60), (60, 1000))
_MIN_ALLELES = 3
# Tail ratios are the geometric mean of this many steps past N-2 or N+2.
_TAIL_STEPS = 4
# Set by hand in HTT_MISEQ, not measured. Its point at CAG 9, and an expansion tail of
# logit 0 below CAG 30, where there are too few molecules above short alleles to measure.
_HAND_SET_BELOW = 30


def _ratio(h: np.ndarray, top: int, bottom: int) -> float:
    if 0 <= top < h.size and 0 <= bottom < h.size and h[bottom] > 0 and h[top] > 0:
        return float(h[top] / h[bottom])
    return float("nan")


def _tail(h: np.ndarray, start: int, step: int) -> float:
    ratios = [_ratio(h, start + step * (j + 1), start + step * j) for j in range(_TAIL_STEPS)]
    found = [r for r in ratios if np.isfinite(r)]
    return float(np.exp(np.mean(np.log(found)))) if found else float("nan")


def _median(values: np.ndarray, log: bool) -> float:
    return float(np.median(_scaled(values, log)))


def _allele_ratios(samples: Sequence[MatrixSample]) -> np.ndarray:
    """One row per labelled allele the recipe keeps: its CAG, then its six ratios."""
    rows = []
    for _, matrix, (a, b) in samples:
        for (cag, ccg), partner in ((a, b), (b, a)):
            if partner[1] == ccg and abs(partner[0] - cag) < _PARTNER:
                continue
            h, n = matrix[ccg - 1], cag - 1
            if h[n] < _MIN_READS:
                continue
            rows.append(
                (
                    cag,
                    _ratio(h, n - 1, n),
                    _ratio(h, n - 2, n - 1),
                    _tail(h, n - 2, -1),
                    _ratio(h, n + 1, n),
                    _ratio(h, n + 2, n + 1),
                    _tail(h, n + 2, 1),
                )
            )
    return np.array(rows)


def _scaled(values: np.ndarray, log: bool) -> np.ndarray:
    """Ratios on the curve's own scale. Log, or logit within 0.001 to 0.97."""
    values = values[np.isfinite(values)]
    if log:
        return np.log(values)
    values = np.clip(values, 1e-3, 0.97)
    return np.log(values) - np.log1p(-values)


def spread_from(samples: Sequence[MatrixSample]) -> list[tuple[float, tuple[float, ...]]]:
    """How much each ratio varies between samples, per length bin."""
    data = _allele_ratios(samples)
    found = []
    for low, high in _BINS:
        bin_rows = data[(data[:, 0] >= low) & (data[:, 0] < high)]
        if len(bin_rows) < _MIN_ALLELES:
            continue
        sds = []
        for k in range(6):
            t = _scaled(bin_rows[:, k + 1], log=k in (0, 3))
            sds.append(1.4826 * float(np.median(np.abs(t - np.median(t)))))
        found.append((float(bin_rows[:, 0].mean()), tuple(sds)))
    return found


def curve_from(samples: Sequence[MatrixSample]) -> StutterCurve:
    """The stutter curve, from the labelled alleles of ``samples``."""
    data = _allele_ratios(samples)
    points = [
        (
            HTT_MISEQ.cag_lengths[0],
            HTT_MISEQ.log_contraction[0],
            HTT_MISEQ.logit_contraction_step[0],
            HTT_MISEQ.logit_contraction_tail[0],
            HTT_MISEQ.log_expansion[0],
            HTT_MISEQ.logit_expansion_step[0],
            HTT_MISEQ.logit_expansion_tail[0],
        )
    ]
    for low, high in _BINS:
        bin_rows = data[(data[:, 0] >= low) & (data[:, 0] < high)]
        if len(bin_rows) < _MIN_ALLELES:
            continue
        mean = float(bin_rows[:, 0].mean())
        medians = [_median(bin_rows[:, k + 1], log=k in (0, 3)) for k in range(6)]
        if mean < _HAND_SET_BELOW:
            medians[5] = 0.0
        points.append((mean, *medians))
    columns = list(zip(*points, strict=True))
    return StutterCurve(*(tuple(float(v) for v in c) for c in columns), spread=HTT_MISEQ.spread)


def _show(name: str, curve: StutterCurve) -> None:
    print(f"  {name}")
    for k, n in enumerate(curve.cag_lengths):
        values = (
            curve.log_contraction[k],
            curve.logit_contraction_step[k],
            curve.logit_contraction_tail[k],
            curve.log_expansion[k],
            curve.logit_expansion_step[k],
            curve.logit_expansion_tail[k],
        )
        print(f"    CAG {n:5.1f}  " + "".join(f"{v:7.2f}" for v in values))


def _score(task: tuple[MatrixSample, StutterCurve]) -> Result:
    sample, curve = task
    return genotype(sample, CallerSettings(stutter=curve))


def _rate(results: Sequence[Result]) -> str:
    hits = sum(_match(r, "exact") for r in results)
    return f"{hits}/{len(results)} ({100 * hits / len(results):.1f}%)" if results else "-"


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--matrix", type=Path, default=DEFAULT_MATRIX)
    parser.add_argument("--splits", type=int, default=1)
    parser.add_argument("--workers", type=int)
    parser.add_argument(
        "--spread-only", action="store_true", help="show the curves and spread, call nothing"
    )
    options = parser.parse_args()
    samples = load(options.matrix)

    print("stutter curves: log N-1/N, logit N-2/N-1, logit tail down, log N+1/N, ...")
    _show("HTT_MISEQ (as shipped)", HTT_MISEQ)
    _show("rebuilt from every sample", curve_from(samples))
    print("  between-sample spread (robust SD, same scales)")
    for n, sds in spread_from(samples):
        print(f"    CAG {n:5.1f}  " + "".join(f"{v:7.2f}" for v in sds))
    print("    priors     " + "".join(f"{v:7.2f}" for v in HTT_MISEQ.spread))
    if options.spread_only:
        return

    tasks: list[tuple[MatrixSample, StutterCurve]] = []
    labels: list[tuple[int, str]] = []  # (split, "unseen" or "own half") per task
    for split in range(options.splits):
        order = np.random.default_rng(split).permutation(len(samples))
        halves = [[samples[i] for i in order[::2]], [samples[i] for i in order[1::2]]]
        curves = [curve_from(half) for half in halves]
        for own in (0, 1):
            for sample in halves[own]:
                tasks += [(sample, curves[1 - own]), (sample, curves[own])]
                labels += [(split, "unseen"), (split, "own half")]
    with ProcessPoolExecutor(options.workers) as pool:
        results = list(pool.map(_score, tasks, chunksize=2))

    print(f"\n{'split':<8}{'curve':<12}{'exact, all':>22}{'exact, labels with reads':>28}")
    for split in range(options.splits):
        for which in ("unseen", "own half"):
            rs = [r for r, label in zip(results, labels, strict=True) if label == (split, which)]
            labelled = [r for r in rs if category(r) != "label at an empty cell"]
            print(f"{split:<8}{which:<12}{_rate(rs):>22}{_rate(labelled):>28}")

    changed = [
        (unseen, own)
        for unseen, own in zip(results[::2], results[1::2], strict=True)
        if _match(unseen, "exact") != _match(own, "exact")
    ]
    print(f"\nsamples whose call changed with the curve: {len(changed)}")
    for unseen, own in changed[:20]:
        print(
            f"  row {unseen.row:>3}: unseen {'right' if _match(unseen, 'exact') else 'wrong'}, "
            f"own half {'right' if _match(own, 'exact') else 'wrong'}"
        )


if __name__ == "__main__":
    main()
