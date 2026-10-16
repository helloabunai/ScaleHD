"""Legacy (ScaleHD 1.x) genotyping!! not ported to this library, ran directly in docker."""

import pytest
from scalehd.counts import SampleCounts
from scalehd.genotype import legacy


def test_legacy_genotyping_isnt_in_this_library() -> None:
    with pytest.raises(NotImplementedError, match=r"runs 1\.x itself"):
        legacy.call_genotype(SampleCounts())
