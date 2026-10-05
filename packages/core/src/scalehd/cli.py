"""Command-line interface: ``scalehd simulate``, ``simulate-run``, ``count``, ``call`` and
``genotype``."""

from __future__ import annotations

import argparse
import json
import sys
from collections.abc import Sequence
from pathlib import Path

from . import __version__
from .counts import SampleCounts, count_fastq
from .genotype import GenotypeCall, NoMoleculesError, call_genotype
from .pairs import DiscordancePolicy
from .simulate import SequencingModel, SimAllele, SimulationSpec, simulate
from .simulate_run import placeholder_samples, write_run
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


def _cmd_simulate_run(args: argparse.Namespace) -> int:
    samples = placeholder_samples(random_count=args.random, seed=args.seed)
    try:
        write_run(args.output, samples, pairs=args.pairs, seed=args.seed, workers=args.workers)
    except FileExistsError as exc:
        print(f"scalehd: {exc}", file=sys.stderr)
        return 1
    print(f"{len(samples)} simulated samples written to {args.output}")
    return 0


def _cmd_count(args: argparse.Namespace) -> int:
    counts = count_fastq(args.r1, args.r2, policy=args.discordant)
    if args.output:
        counts.write_json(args.output)
    _print_summary(counts, args.top)
    return 0


def _cmd_call(args: argparse.Namespace) -> int:
    return _call(SampleCounts.read_json(args.counts), args.output)


def _cmd_genotype(args: argparse.Namespace) -> int:
    counts = count_fastq(args.r1, args.r2, policy=args.discordant)
    if args.counts:
        counts.write_json(args.counts)
    _print_reads(counts)
    print()
    return _call(counts, args.output)


def _call(counts: SampleCounts, output: Path | None) -> int:
    try:
        call = call_genotype(counts)
    except NoMoleculesError as exc:
        print(f"scalehd: {exc}", file=sys.stderr)
        return 1
    if output:
        output.write_text(json.dumps(call.to_dict(), indent=2) + "\n")
    _print_call(call)
    return 0


def _print_call(call: GenotypeCall) -> None:
    print(f"genotype  {call.label}")
    print(
        f"posterior {call.posterior:.4f}  quality {call.quality:.1f}  molecules {call.molecules:,}"
    )
    if call.flags:
        print("flags     " + ", ".join(call.flags))
    for allele in call.alleles if not call.homozygous else call.alleles[:1]:
        line = f"  {allele.label:<16} share {allele.fraction:.2f}"
        if allele.cag_estimate:
            best, low, high = allele.cag_estimate
            line += f"  CAG estimate {best} ({low}-{high}, rough)"
        if allele.backward_slippage is not None:
            line += f"  slippage {allele.backward_slippage:.3f}"
        if allele.somatic_mosaicism is not None:
            line += f"  mosaicism {allele.somatic_mosaicism:.3f}"
        print(line)
    for genotype, posterior in call.alternatives[:2]:
        print(f"  alternative {genotype}  {posterior:.2e}")
    for structure, n in call.unexplained:
        print(f"  unexplained peak {structure} ({n:,} molecules)")


def _print_reads(counts: SampleCounts) -> None:
    print(
        f"molecules {counts.molecules:,}  complete {sum(counts.complete.values()):,}  "
        f"partial {sum(counts.partial.values()):,}  dropped {counts.dropped:,}  "
        f"unusable {counts.unusable:,}"
    )
    print("reads: " + ", ".join(f"{k} {v:,}" for k, v in sorted(counts.read_outcomes.items())))
    if counts.discordant:
        print(
            "discordant read base pairings: "
            + ", ".join(f"{k} {v:,}" for k, v in counts.discordant.items())
        )


def _print_summary(counts: SampleCounts, top: int) -> None:
    complete = sum(counts.complete.values())
    _print_reads(counts)
    print(f"\n{'structure':<16}{'molecules':>11}{'%':>8}")
    for structure, n in counts.top(top):
        print(f"{structure.label:<16}{n:>11,}{100 * n / max(complete, 1):>8.2f}")
    if counts.partial:
        print(f"\n{'partial':<16}{'molecules':>11}")
        for observation, n in counts.partial.most_common(min(top, 5)):
            print(f"{observation.label:<16}{n:>11,}")


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

    run = sub.add_parser(
        "simulate-run", help="write a whole simulated MiSeq run folder of placeholder data"
    )
    run.add_argument("output", type=Path, help="run folder to create (new or empty)")
    run.add_argument("-n", "--pairs", type=int, default=20_000, help="read pairs per sample")
    run.add_argument(
        "--random", type=int, default=20, help="random genotypes on top of the fixed samples"
    )
    run.add_argument("--seed", type=int, default=1)
    run.add_argument("--workers", type=int, help="samples simulated at once (default: every core)")
    run.set_defaults(func=_cmd_simulate_run)

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
        help="what to do when read base pairings disagree",
    )
    count.set_defaults(func=_cmd_count)

    call = sub.add_parser("call", help="call a genotype from a counts JSON")
    call.add_argument("counts", type=Path)
    call.add_argument("-o", "--output", type=Path, help="write the call as JSON here")
    call.set_defaults(func=_cmd_call)

    genotype = sub.add_parser("genotype", help="count and call in one step")
    genotype.add_argument("r1", type=Path)
    genotype.add_argument("r2", type=Path, nargs="?")
    genotype.add_argument("-o", "--output", type=Path, help="write the call as JSON here")
    genotype.add_argument("--counts", type=Path, help="also write the counts JSON here")
    genotype.add_argument(
        "--discordant",
        type=DiscordancePolicy,
        default=DiscordancePolicy.DROP,
        choices=list(DiscordancePolicy),
        help="what to do when read base pairings disagree",
    )
    genotype.set_defaults(func=_cmd_genotype)
    return parser


def main(argv: Sequence[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    return int(args.func(args))


if __name__ == "__main__":
    sys.exit(main())
