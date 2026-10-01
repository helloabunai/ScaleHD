"""Command-line interface: ``scalehd simulate`` and ``scalehd count``."""

from __future__ import annotations

import argparse
import sys
from collections.abc import Sequence
from pathlib import Path

from . import __version__
from .counts import SampleCounts, count_fastq
from .pairs import DiscordancePolicy
from .simulate import SequencingModel, SimAllele, SimulationSpec, simulate
from .structure import AlleleStructure


def _allele(text: str) -> tuple[AlleleStructure, float]:
    label, _, abundance = text.partition(":")
    try:
        return AlleleStructure.from_label(label), float(abundance or 1.0)
    except ValueError as exc:
        raise argparse.ArgumentTypeError(str(exc)) from exc


def _cmd_simulate(args: argparse.Namespace) -> int:
    alleles = tuple(
        SimAllele(
            structure,
            abundance,
            somatic_fraction=args.somatic_fraction
            if structure.cag >= args.somatic_min_cag
            else 0.0,
            somatic_mean=args.somatic_mean,
        )
        for structure, abundance in args.allele
    )
    spec = SimulationSpec(
        alleles,
        pairs=args.pairs,
        off_target=args.off_target,
        sequencing=SequencingModel(read_length=args.read_length),
        seed=args.seed,
    )
    paths = simulate(spec).write(args.output, args.name)
    for path in paths:
        print(path)
    return 0


def _cmd_count(args: argparse.Namespace) -> int:
    counts = count_fastq(args.r1, args.r2, policy=args.discordant)
    if args.output:
        counts.write_json(args.output)
    _print_summary(counts, args.top)
    return 0


def _print_summary(counts: SampleCounts, top: int) -> None:
    complete = sum(counts.complete.values())
    print(
        f"molecules {counts.molecules:,}  complete {complete:,}  "
        f"partial {sum(counts.partial.values()):,}  dropped {counts.dropped:,}  "
        f"unusable {counts.unusable:,}"
    )
    print(f"\n{'structure':<16}{'molecules':>11}{'%':>8}")
    for structure, n in counts.top(top):
        print(f"{structure.label:<16}{n:>11,}{100 * n / max(complete, 1):>8.2f}")
    if counts.partial:
        print(f"\n{'partial':<16}{'molecules':>11}")
        for observation, n in counts.partial.most_common(min(top, 5)):
            print(f"{observation.label:<16}{n:>11,}")
    print("\nreads: " + ", ".join(f"{k} {v:,}" for k, v in sorted(counts.read_outcomes.items())))
    if counts.discordant:
        print("discordant mates: " + ", ".join(f"{k} {v:,}" for k, v in counts.discordant.items()))


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(prog="scalehd", description=__doc__)
    parser.add_argument("--version", action="version", version=f"%(prog)s {__version__}")
    sub = parser.add_subparsers(dest="command", required=True)

    sim = sub.add_parser("simulate", help="write simulated paired FASTQ with a truth file")
    sim.add_argument(
        "-a",
        "--allele",
        type=_allele,
        action="append",
        required=True,
        metavar="LABEL[:ABUNDANCE]",
        help="e.g. 43_1_1_7_2 or 43_1_1_7_2:0.8",
    )
    sim.add_argument("-n", "--pairs", type=int, default=20_000)
    sim.add_argument("--read-length", type=int, default=300)
    sim.add_argument("--somatic-fraction", type=float, default=0.0)
    sim.add_argument("--somatic-mean", type=float, default=3.0)
    sim.add_argument(
        "--somatic-min-cag",
        type=int,
        default=36,
        help="only alleles at least this long get somatic expansion",
    )
    sim.add_argument("--off-target", type=float, default=0.0)
    sim.add_argument("--seed", type=int, default=0)
    sim.add_argument("-o", "--output", type=Path, required=True, help="output directory")
    sim.add_argument("--name", default="sample", help="sample name used in file names")
    sim.set_defaults(func=_cmd_simulate)

    count = sub.add_parser("count", help="tally molecule repeat structures from FASTQ")
    count.add_argument("r1", type=Path)
    count.add_argument("r2", type=Path, nargs="?")
    count.add_argument("-o", "--output", type=Path, help="write counts JSON here")
    count.add_argument("--top", type=int, default=10)
    count.add_argument(
        "--discordant",
        type=DiscordancePolicy,
        default=DiscordancePolicy.DROP,
        choices=list(DiscordancePolicy),
        help="what to do when mates disagree",
    )
    count.set_defaults(func=_cmd_count)
    return parser


def main(argv: Sequence[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    return int(args.func(args))


if __name__ == "__main__":
    sys.exit(main())
