"""Legacy (ScaleHD 1.x) genotyping: a placeholder until it is brought over."""

import pytest
from scalehd.counts import SampleCounts
from scalehd.genotype import legacy


def test_legacy_genotyping_is_not_available_yet() -> None:
    with pytest.raises(NotImplementedError, match="not available yet"):
        legacy.call_genotype(SampleCounts())
