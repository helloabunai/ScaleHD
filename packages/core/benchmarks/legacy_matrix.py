"""Genotype the ScaleHD 1.x labelled training matrix and compare with its labels.

Each row of ``legacy/ScaleHD/train/raw_matrix.csv`` is one MiSeq sample: read counts
for CAG 1-100 x CCG 1-20, as aligned by ScaleHD 1.x, and a label such as
``CAG43CCG7-CAG19CCG10`` that was probably curated by hand. The matrix has no
intervening-sequence information, so every structure is taken as n_1_1_m_2. Atypical
alleles were re-aligned to custom references, which we aren't doing.

The stutter priors were calibrated on this same matrix. That is four smooth curves
fitted to thousands of alleles, so the effect on this comparison is small but not
zero.

Run with ``uv run python packages/core/benchmarks/legacy_matrix.py``.
"""

from __future__ import annotations

import argparse
import csv
import json
import re
from collections import Counter
from concurrent.futures import ProcessPoolExecutor
from dataclasses import asdict, dataclass
from functools import partial
from pathlib import Path

import numpy as np
from scalehd.counts import SampleCounts
from scalehd.genotype import CallerSettings, call_genotype
from scalehd.structure import AlleleStructure

DEFAULT_MATRIX = Path(__file__).parents[3] / "legacy/ScaleHD/train/raw_matrix.csv"
LABEL = re.compile(r"^CAG(\d+)CCG(\d+)-CAG(\d+)CCG(\d+)$")


@dataclass(frozen=True)
class Result:
    row: int
    truth: tuple[tuple[int, int], tuple[int, int]]
    called: tuple[tuple[int, int], tuple[int, int]]
    posterior: float
    quality: float
    flags: tuple[str, ...]
    label_reads: tuple[float, float]


def load(path: Path) -> list[tuple[int, np.ndarray, tuple[tuple[int, int], tuple[int, int]]]]:
    samples = []
    with path.open() as handle:
        reader = csv.reader(handle)
        next(reader)
        for i, row in enumerate(reader):
            match = LABEL.match(row[-1])
            if not match:
                continue
            matrix = np.array([float(x) for x in row[1:-1]]).reshape(20, 100)  # [ccg-1, cag-1]
            a, b = (int(match[1]), int(match[2])), (int(match[3]), int(match[4]))
            samples.append((i, matrix, (a, b)))
    return samples


def counts_from_matrix(matrix: np.ndarray) -> SampleCounts:
    counts = SampleCounts()
    for ccg, cag in zip(*np.nonzero(matrix), strict=True):
        n = round(float(matrix[ccg, cag]))
        if n:
            counts.complete[AlleleStructure(int(cag) + 1, 1, 1, int(ccg) + 1, 2)] += n
    counts.molecules = sum(counts.complete.values())
    return counts


def genotype(
    sample: tuple[int, np.ndarray, tuple[tuple[int, int], tuple[int, int]]],
    settings: CallerSettings | None = None,
) -> Result:
    i, matrix, truth = sample
    call = call_genotype(counts_from_matrix(matrix), settings)
    called = tuple(sorted((a.allele.structure.cag, a.allele.structure.ccg) for a in call.alleles))
    reads = tuple(float(matrix[ccg - 1, cag - 1]) for cag, ccg in truth)
    return Result(
        i,
        tuple(sorted(truth)),  # type: ignore[arg-type]
        called,  # type: ignore[arg-type]
        call.posterior,
        call.quality,
        tuple(str(f) for f in call.flags),
        reads,  # type: ignore[arg-type]
    )


def category(r: Result) -> str:
    (a_cag, a_ccg), (b_cag, b_ccg) = r.truth
    if 0 in r.label_reads:
        return "label at an empty cell"
    if r.truth[0] == r.truth[1]:
        return "homozygous"
    if a_ccg == b_ccg and abs(a_cag - b_cag) <= 2:
        return "within 2 CAG, same CCG"
    if max(a_cag, b_cag) >= 36:
        return "expanded (>=36)"
    return "other"


def report(results: list[Result]) -> None:
    def rate(rs: list[Result], ok: str) -> str:
        if not rs:
            return "-"
        hits = sum(_match(r, ok) for r in rs)
        return f"{hits}/{len(rs)} ({100 * hits / len(rs):.1f}%)"

    print(f"samples {len(results)}\n")
    print(f"{'subset':<28}{'exact':>20}{'CAG only':>20}{'within ±1 CAG':>20}")
    groups: dict[str, list[Result]] = {"all": results}
    for r in results:
        groups.setdefault(category(r), []).append(r)
    for name, rs in groups.items():
        print(f"{name:<28}{rate(rs, 'exact'):>20}{rate(rs, 'cag'):>20}{rate(rs, 'near'):>20}")

    print(f"\n{'call quality':<28}{'samples':>10}{'exact':>12}")
    for low, high in ((0, 10), (10, 20), (20, 40), (40, 99.5), (99.5, 100)):
        rs = [r for r in results if low <= r.quality < high or (high == 100 and r.quality >= 99)]
        hits = sum(_match(r, "exact") for r in rs)
        share = f"{100 * hits / len(rs):.1f}%" if rs else "-"
        print(f"{f'{low}-{high}':<28}{len(rs):>10}{share:>12}")

    flags = Counter(f for r in results for f in r.flags)
    wrong_flags = Counter(f for r in results if not _match(r, "exact") for f in r.flags)
    print("\nflags (all / among mismatches):")
    for flag, n in flags.most_common():
        print(f"  {flag:<20}{n:>6}{wrong_flags[flag]:>6}")

    print("\nmismatches outside 'label at an empty cell' (first 25):")
    shown = 0
    for r in results:
        if _match(r, "exact") or category(r) == "label at an empty cell":
            continue
        print(
            f"  row {r.row:>3}: label {_fmt(r.truth):<22} called {_fmt(r.called):<22} "
            f"Q {r.quality:5.1f}  {','.join(r.flags)}"
        )
        shown += 1
        if shown == 25:
            break


def _fmt(pair: tuple[tuple[int, int], tuple[int, int]]) -> str:
    return "/".join(f"{cag}:{ccg}" for cag, ccg in pair)


def _match(r: Result, how: str) -> bool:
    if how == "exact":
        return r.truth == r.called
    if how == "cag":
        return sorted(c for c, _ in r.truth) == sorted(c for c, _ in r.called)
    truth = sorted(c for c, _ in r.truth)
    called = sorted(c for c, _ in r.called)
    return all(abs(t - c) <= 1 for t, c in zip(truth, called, strict=True))


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--matrix", type=Path, default=DEFAULT_MATRIX)
    parser.add_argument("--limit", type=int)
    parser.add_argument("--workers", type=int)
    parser.add_argument("--window", help="stutter window as BELOW,ABOVE, e.g. 8,12")
    parser.add_argument("--dump", type=Path, help="write every result here as JSON")
    options = parser.parse_args()
    settings = CallerSettings()
    if options.window:
        below, above = (int(v) for v in options.window.split(","))
        settings = CallerSettings(stutter_window=(below, above))
    samples = load(options.matrix)[: options.limit]
    with ProcessPoolExecutor(options.workers) as pool:
        results = list(pool.map(partial(genotype, settings=settings), samples, chunksize=4))
    if options.dump:
        options.dump.write_text(json.dumps([asdict(r) for r in results]))
    report(results)


if __name__ == "__main__":
    main()
