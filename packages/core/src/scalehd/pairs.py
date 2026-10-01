"""Join R1 and R2 observations of the same molecule.

Attempting different logic here compared to ScaleHD 1.x

Both strands read the same template, PCR stutter is a thing that exists so
single substitutions almost never coincide at the same base pair position.

However single sub/ins/del can move a repeat tract start/end boundary with the
flanks (e.g. G>A sub in last CAG reads as extra CAACAG intervening seq), so default
a moluecule whose pair/mates disagree is dropped naively.

Where only one read saw a field to the end-point (e.g. very long CAG alleles), assume
truth. R1 reads CAG and intervening sequence where quality is high. R2 does the same
for CCG and CCT.
"""

from __future__ import annotations

from dataclasses import dataclass
from enum import StrEnum

from .structure import FIELDS, FieldStatus, Observation

# Per field (cag, caacag, ccgcca, ccg, cct): True to prefer R1, False to prefer R2.
DEFAULT_PREFER_R1: tuple[bool, ...] = (True, True, True, False, False)


class DiscordancePolicy(StrEnum):
    DROP = "drop"  # discard molecules whose mates disagree
    # Keep the preferred mate's value per field. Can stitch together combinations
    # that neither mate saw. intended for development comparison, not routine use.
    PREFER = "prefer"


@dataclass(frozen=True, slots=True)
class PairResult:
    observation: Observation | None
    discordant: tuple[bool, ...]

    @property
    def is_discordant(self) -> bool:
        return any(self.discordant)


def join_mates(
    r1: Observation | None,
    r2: Observation | None,
    *,
    policy: DiscordancePolicy = DiscordancePolicy.DROP,
    prefer_r1: tuple[bool, ...] = DEFAULT_PREFER_R1,
) -> PairResult:
    if r1 is None and r2 is None:
        return PairResult(None, (False,) * len(FIELDS))

    counts: list[int] = []
    status: list[FieldStatus] = []
    discordant: list[bool] = []
    for k in range(len(FIELDS)):
        ordered = (r1, r2) if prefer_r1[k] else (r2, r1)
        seen = [(o.counts[k], o.status[k]) for o in ordered if o is not None]
        exact = [c for c, s in seen if s is FieldStatus.EXACT]
        bounds = [c for c, s in seen if s is FieldStatus.LOWER_BOUND]
        if exact:
            counts.append(exact[0])
            status.append(FieldStatus.EXACT)
            discordant.append(len(set(exact)) > 1 or any(b > exact[0] for b in bounds))
        elif bounds:
            counts.append(max(bounds))
            status.append(FieldStatus.LOWER_BOUND)
            discordant.append(False)
        else:
            counts.append(0)
            status.append(FieldStatus.UNOBSERVED)
            discordant.append(False)

    flags = tuple(discordant)
    if any(flags) and policy is DiscordancePolicy.DROP:
        return PairResult(None, flags)
    return PairResult(Observation(tuple(counts), tuple(status)), flags)  # type: ignore[arg-type]
