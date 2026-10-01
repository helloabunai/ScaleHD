from pathlib import Path

import pytest
from scalehd.cli import main
from scalehd.seqio import (
    FastqFormatError,
    FastqRecord,
    open_text,
    read_fastq,
    read_pairs,
    reverse_complement,
    write_fastq,
)


def test_reverse_complement() -> None:
    assert reverse_complement("CAGCCGN") == "NCGGCTG"


def test_pair_id_strips_mate_markers() -> None:
    assert FastqRecord("read1/1", "A", "I").pair_id == "read1"
    assert FastqRecord("read1 2:N:0:ACGT", "A", "I").pair_id == "read1"


@pytest.mark.parametrize("suffix", [".fastq", ".fastq.gz"])
def test_fastq_round_trip(tmp_path: Path, suffix: str) -> None:
    records = [FastqRecord("a", "ACGT", "IIII"), FastqRecord("b", "GG", "##")]
    path = tmp_path / f"x{suffix}"
    with open_text(path, "wt") as handle:
        write_fastq(handle, records)
    assert list(read_fastq(path)) == records


def test_unsynced_mates_rejected(tmp_path: Path) -> None:
    for name, ids in (("r1.fq", ["a", "b"]), ("r2.fq", ["a", "c"])):
        with open_text(tmp_path / name, "wt") as handle:
            write_fastq(handle, [FastqRecord(i, "A", "I") for i in ids])
    with pytest.raises(FastqFormatError, match="out of sync"):
        list(read_pairs(tmp_path / "r1.fq", tmp_path / "r2.fq"))


def test_malformed_fastq(tmp_path: Path) -> None:
    path = tmp_path / "bad.fq"
    path.write_text("@a\nACGT\n+\nII\n")
    with pytest.raises(FastqFormatError, match="length"):
        list(read_fastq(path))


def test_cli_simulate_then_count(tmp_path: Path, capsys: pytest.CaptureFixture[str]) -> None:
    assert (
        main(
            [
                "simulate",
                "-a",
                "17_1_1_7_2",
                "-a",
                "43_1_1_7_2:0.8",
                "-n",
                "400",
                "-o",
                str(tmp_path),
                "--name",
                "s",
            ]
        )
        == 0
    )
    out = tmp_path / "s.counts.json"
    assert (
        main(
            [
                "count",
                str(tmp_path / "s_R1.fastq.gz"),
                str(tmp_path / "s_R2.fastq.gz"),
                "-o",
                str(out),
            ]
        )
        == 0
    )
    printed = capsys.readouterr().out
    assert "17_1_1_7_2" in printed
    assert "43_1_1_7_2" in printed
    assert out.exists()
