"""Self tests to ensure the legacy wrapper isn't changing anything"""

from __future__ import annotations

import gzip
import importlib.util
import io
import json
import logging
import re
import xml.etree.ElementTree as ET
from pathlib import Path
from types import ModuleType
from typing import Any

import numpy as np
import pytest

RUNNER_FILE = Path(__file__).resolve().parents[1] / "scalehd_legacy_run.py"
LEGACY = Path(__file__).resolve().parents[3] / "legacy" / "ScaleHD"


def _load() -> ModuleType:
    spec = importlib.util.spec_from_file_location("scalehd_legacy_run", RUNNER_FILE)
    assert spec is not None
    assert spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


runner = _load()

# How 1.x writes its log lines (ScaleHD/sherpa.py), colours included.
BOLD, YELLOW, RED, END = "\033[1m", "\033[93m", "\033[91m", "\033[0m"


def _config() -> dict[str, dict[str, str]]:
    return {
        "instance_flags": {"demultiplex": "False", "quality_control": "True"},
        "demultiplex_flags": {},
        "trim_flags": {"trim_type": "Adapter", "error_tolerance": "0.39"},
        "alignment_flags": {"min_seed_length": "19"},
        "prediction_flags": {"snp_observation_threshold": "2", "quality_cutoff": "0"},
    }


@pytest.fixture
def request_data(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> dict[str, Any]:
    root = tmp_path / "opt-legacy"
    (root / "references").mkdir(parents=True)
    for name in (runner.FORWARD_REFERENCE, runner.REVERSE_REFERENCE):
        (root / "references" / name).write_text(">ref\nACGT\n")
    monkeypatch.setattr(runner, "LEGACY_ROOT", root)
    reads = tmp_path / "reads"
    reads.mkdir()
    (reads / "a_1.fq").write_text("@r\nACGT\n+\nIIII\n")
    (reads / "a_2.fq").write_text("@r\nACGT\n+\nIIII\n")
    return {
        "label": "sample-1",
        "r1": str(reads / "a_1.fq"),
        "r2": str(reads / "a_2.fq"),
        "threads": 4,
        "config": _config(),
    }


def test_a_request_gets_the_references_1x_ships(request_data: dict[str, Any]) -> None:
    request = runner.read_request(request_data)

    assert request["label"] == "sample-1"
    assert request["threads"] == 4
    assert request["enshrine"] is False
    assert request["broadscope"] is False
    assert Path(request["forward_reference"]).name == "4k-HD-INTER.fasta"
    assert Path(request["reverse_reference"]).name == "20-TypicalReverse.fasta"
    assert request["config"] == _config()


@pytest.mark.parametrize(
    ("change", "complaint"),
    [
        ({"label": "a.b"}, "label"),
        ({"label": "a b"}, "label"),
        ({"label": "_a"}, "label"),
        ({"label": ""}, "label"),
        ({"r1": "/no/such/file"}, "r1"),
        ({"threads": 0}, "threads"),
        ({"threads": True}, "threads"),
        ({"threads": "4"}, "threads"),
        ({"enshrine": "yes"}, "enshrine"),
        ({"forward_reference": "/no/such.fasta"}, "forward_reference"),
    ],
)
def test_a_bad_request_is_refused(
    request_data: dict[str, Any], change: dict[str, Any], complaint: str
) -> None:
    with pytest.raises(runner.RequestError, match=complaint):
        runner.read_request({**request_data, **change})


def test_a_reference_must_be_named_so_1x_can_use_it(
    request_data: dict[str, Any], tmp_path: Path
) -> None:
    reference = tmp_path / "my reference.fasta"
    reference.write_text(">ref\nACGT\n")
    with pytest.raises(runner.RequestError, match="characters"):
        runner.read_request({**request_data, "forward_reference": str(reference)})


@pytest.mark.parametrize(
    "config",
    [
        {k: v for k, v in _config().items() if k != "trim_flags"},
        {**_config(), "extra_flags": {}},
        {**_config(), "trim_flags": {"quality_threshold": 5}},
    ],
)
def test_a_bad_config_is_refused(request_data: dict[str, Any], config: dict[str, Any]) -> None:
    with pytest.raises(runner.RequestError, match="config"):
        runner.read_request({**request_data, "config": config})


def test_1x_can_demultiplex_the_sample_itself(request_data: dict[str, Any]) -> None:
    config = _config()
    config["instance_flags"]["demultiplex"] = "True"

    request = runner.read_request({**request_data, "config": config})

    assert request["config"]["instance_flags"]["demultiplex"] == "True"


def test_the_config_is_written_as_given_in_the_dtds_order() -> None:
    config = _config()
    config["trim_flags"]["forward_adapter"] = 'A"<&>'

    xml = runner.config_xml(config, "/data/in", "/refs/fw.fasta", "/refs/rv.fasta")
    root = ET.fromstring(xml)

    assert root.tag == "config"
    assert root.attrib == {
        "data_dir": "/data/in",
        "forward_reference": "/refs/fw.fasta",
        "reverse_reference": "/refs/rv.fasta",
    }
    assert [child.tag for child in root] == list(runner.ELEMENTS)
    assert [child.attrib for child in root] == [config[element] for element in runner.ELEMENTS]
    assert "\t<demultiplex_flags/>\n" in xml


def test_the_elements_are_in_the_order_1x_dtd_has_them() -> None:
    dtd = (LEGACY / "config" / "config.dtd").read_text()
    declared = re.search(r"<!ELEMENT config\s+\(([^)]*)\)>", dtd)

    assert declared is not None
    assert tuple(name.strip() for name in declared.group(1).split(",")) == runner.ELEMENTS


def test_reads_are_copied_under_the_names_1x_pairs(tmp_path: Path) -> None:
    r1 = tmp_path / "reads.1.FQ.GZ"
    with gzip.open(r1, "wt") as file:
        file.write("@r\nACGT\n+\nIIII\n")
    r2 = tmp_path / "whatever"
    r2.write_text("@r\nTTTT\n+\nIIII\n")

    staged = runner.stage_reads("s1", str(r1), str(r2), tmp_path / "input")

    assert [path.name for path in staged] == ["s1_R1.fastq.gz", "s1_R2.fastq"]
    assert staged[0].read_bytes() == r1.read_bytes()
    assert staged[1].read_bytes() == r2.read_bytes()
    assert r1.is_file()
    assert r2.is_file()


def test_the_output_folder_must_be_empty(tmp_path: Path) -> None:
    assert runner.prepare_out_dir(tmp_path / "new" / "out").is_dir()
    assert runner.prepare_out_dir(tmp_path / "new" / "out").is_dir()
    (tmp_path / "new" / "out" / "file").write_text("")

    with pytest.raises(runner.RequestError, match="empty"):
        runner.prepare_out_dir(tmp_path / "new" / "out")


def test_the_output_folder_must_be_safe_for_1xs_shell_commands(tmp_path: Path) -> None:
    with pytest.raises(runner.RequestError, match="characters"):
        runner.prepare_out_dir(tmp_path / "my output")


def test_1xs_programs_come_first_and_keep_to_the_threads() -> None:
    environ = {"PATH": "/usr/bin", "PYTHONPATH": "/elsewhere", "HOME": "/home/x"}

    env = runner.legacy_environment(environ, 3, Path("/opt/legacy"))

    assert env["PATH"].split(":") == [
        "/opt/legacy/bin",
        "/opt/legacy/python/bin",
        "/opt/legacy/tools/bin",
        "/usr/bin",
    ]
    assert env["OMP_NUM_THREADS"] == env["OPENBLAS_NUM_THREADS"] == env["MKL_NUM_THREADS"] == "3"
    assert "PYTHONPATH" not in env
    assert env["HOME"] == "/home/x"
    assert environ["PYTHONPATH"] == "/elsewhere"


def _log(colour: str, message: str) -> str:
    return f"{colour}shd__ {END}{message}"


def test_the_listener_reports_stages_as_1x_logs_them() -> None:
    events = io.StringIO()
    listener = runner.Listener(events)
    for colour, message in [
        (BOLD, "Indexing reference(s) before initialising sample pair cycle.."),
        (BOLD, "Demultiplexing reads.. (this may take some time)"),
        (BOLD, "Processing sequence pair: 1/1 (s1)"),
        (YELLOW, "Executing sequence quality control workflow.."),
        (YELLOW, "Executing alignment workflow.."),
        (BOLD, "Scanning for atypical alleles.."),
        (YELLOW, "User specified sequence re-alignment. Generating custom reference.."),
        (YELLOW, "Genotyping alleles.. "),
        (YELLOW, "Calling SNPs.. "),
    ]:
        listener.heard(logging.INFO, _log(colour, message))
    listener.heard(logging.DEBUG, "Took 0.000s to load font")

    lines = [json.loads(line) for line in events.getvalue().splitlines()]
    assert [line["stage"] for line in lines] == [
        "indexing",
        "demultiplexing",
        "quality_control",
        "alignment",
        "atypical_scanning",
        "realignment",
        "genotyping",
        "snp_calling",
    ]
    assert listener.stages == [line["stage"] for line in lines]
    assert listener.failure is None


def test_the_listener_keeps_why_the_sample_failed() -> None:
    listener = runner.Listener(None)
    error = "<200 aligned reads in Allele #1. Data un-usable."
    listener.heard(logging.INFO, _log(RED, f"Atypical scanning failure on s1: {error}\n"))
    listener.heard(logging.INFO, _log(RED, "Genotyping failure on s1: later"))

    assert listener.failure == error


def test_the_listener_knows_a_demultiplexed_sample() -> None:
    listener = runner.Listener(None)
    listener.heard(logging.INFO, _log(RED, "Alignment failure on DMPX_s1: no reads\n"))

    assert listener.failure == "no reads"


def test_the_listener_keeps_1xs_errors() -> None:
    listener = runner.Listener(None)
    dtd = "DTD validation failure /out/config.xml: No declaration for attribute adapter"
    listener.heard(logging.ERROR, dtd)
    listener.heard(logging.ERROR, _log(RED, "XML Config: Failure, exiting."))
    listener.heard(logging.WARNING, _log(RED, "Alignment contains too few reads."))

    assert listener.errors == [dtd, "XML Config: Failure, exiting."]


def test_the_listener_hears_1xs_logging() -> None:
    events = io.StringIO()
    listener = runner.Listener(events)
    logger = logging.getLogger("test_legacy_runner")
    logger.addHandler(listener)
    try:
        logger.warning("%s", _log(YELLOW, "Calling SNPs.. "))
    finally:
        logger.removeHandler(listener)

    assert json.loads(events.getvalue()) == {"stage": "snp_calling"}


class Allele:
    def get_cag(self) -> Any:
        return np.int64(17)

    def get_fwarray(self) -> Any:
        return np.array([0, 3, 1])

    def get_validation(self) -> Any:
        return np.bool_(True)

    def get_somaticmosaicism(self) -> float:
        return float("nan")

    def get_backwardsslippage(self) -> Any:
        return np.float64("-inf")

    def set_cag(self, value: int) -> None:
        pass


class Sample:
    def get_label(self) -> str:
        return "s1"

    def get_primaryallele(self) -> Allele:
        return Allele()

    def get_secondaryallele(self) -> None:
        return None

    def get_trimreport(self) -> list[str]:
        return ["/out/run/s1/SeqQC/s1_R1_TrimmingReport.txt"]

    def get_exceptionraised(self) -> str:
        return "N/A"

    def get_flag(self, name: str) -> bool:
        return True


# As if they were 1.x's own classes.
Allele.__module__ = Sample.__module__ = "ScaleHD.__allelecontainer"


def test_every_getter_value_is_read_as_json() -> None:
    values = runner.to_json(Sample())

    assert values == {
        "exceptionraised": "N/A",
        "label": "s1",
        "primaryallele": {
            "backwardsslippage": "-inf",
            "cag": 17,
            "fwarray": [0, 3, 1],
            "somaticmosaicism": "nan",
            "validation": True,
        },
        "secondaryallele": None,
        "trimreport": ["/out/run/s1/SeqQC/s1_R1_TrimmingReport.txt"],
    }
    json.dumps(values, allow_nan=False)


def test_anything_else_is_kept_by_its_repr() -> None:
    assert runner.to_json({1: (2, Path("/x"))}) == {
        "1": [2, {"type": "PosixPath", "repr": "PosixPath('/x')"}]
    }


class Instance:
    def __init__(self, samples: list[Any]) -> None:
        self.instance_objects = samples


class FailedSample(Sample):
    def get_exceptionraised(self) -> str:
        return "DSP"


def test_a_complete_sample() -> None:
    found = runner.outcome(Instance([Sample()]), runner.Listener(None))

    assert found["status"] == "complete"
    assert found["failed_stage"] is None
    assert found["error"] is None


def test_a_sample_1x_could_not_finish() -> None:
    listener = runner.Listener(None)
    listener.heard(logging.INFO, _log(RED, "Atypical scanning failure on s1: too few reads"))

    found = runner.outcome(Instance([FailedSample()]), listener)

    assert found["status"] == "failed"
    assert found["failed_stage"] == "DSP"
    assert found["error"] == "too few reads"


def test_1x_stopping_before_the_sample() -> None:
    listener = runner.Listener(None)
    listener.heard(logging.CRITICAL, _log(RED, "Missing binary: freebayes!"))
    listener.heard(logging.ERROR, _log(RED, "Detected missing library from system/$PATH."))

    for instance in (None, Instance([])):
        found = runner.outcome(instance, listener)
        assert found["status"] == "failed"
        assert found["sample"] is None
        assert found["error"] == (
            "Missing binary: freebayes!; Detected missing library from system/$PATH."
        )


def test_demultiplexing_that_wrote_nothing_is_said_so(
    request_data: dict[str, Any], tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    config = _config()
    config["instance_flags"]["demultiplex"] = "True"
    request = runner.read_request({**request_data, "config": config})
    out = runner.prepare_out_dir(tmp_path / "out")

    def run_scalehd(request: dict[str, Any], out_dir: Path, listener: Any) -> None:
        # As 1.x does when batchadapt fails i.e. an empty folder, and no samples.
        (out_dir / "input_demultiplexed").mkdir()

    monkeypatch.setattr(runner, "run_scalehd", run_scalehd)
    found = runner.run(request, out, runner.Listener(None))

    assert found["status"] == "failed"
    assert found["error"] == "demultiplexing (batchadapt) wrote no reads"
    assert sorted(path.name for path in out.iterdir()) == ["config.xml"]


def test_versions_say_where_1x_came_from(tmp_path: Path) -> None:
    (tmp_path / "UPSTREAM").write_text("url=https://example.org/x\ncommit=abc\ndate=2024-11-27\n")

    found = runner.versions(tmp_path)

    assert found["upstream"] == {
        "url": "https://example.org/x",
        "commit": "abc",
        "date": "2024-11-27",
    }


def test_a_bad_request_stops_before_anything_runs(
    request_data: dict[str, Any], tmp_path: Path, capsys: pytest.CaptureFixture[str]
) -> None:
    request = tmp_path / "request.json"
    request.write_text(json.dumps({**request_data, "threads": 0}))

    assert runner.main([str(request), str(tmp_path / "out")]) == 2
    assert "threads" in capsys.readouterr().err
    assert not (tmp_path / "out").exists()


# What the runner relies on in 1.x, checked against legacy/ so an update can't break the
# wrapper quietly

SHERPA = (LEGACY / "sherpa.py").read_text()


@pytest.mark.parametrize(("start", "stage"), runner.STAGES)
def test_1x_still_logs_each_stage_starting(start: str, stage: str) -> None:
    assert f"'{start}" in SHERPA


def test_1x_still_says_why_a_sample_failed_the_same_way() -> None:
    failures = [line for line in SHERPA.splitlines() if "failure on " in line]

    assert len(failures) == 7
    for line in failures:
        assert ".format(clr.red" in line
        assert re.search(r"failure on ',\s*seqpair_lbl,\s*str\(e\)\)", line), line


def test_threads_are_only_set_where_the_runner_sets_them() -> None:
    modules = {
        "ScaleHD." + ".".join(path.relative_to(LEGACY).with_suffix("").parts)
        for path in LEGACY.rglob("*.py")
        if re.search(r"\bTHREADS\b|cpu_count", path.read_text())
    }

    assert modules == {*runner.THREAD_MODULES, "ScaleHD.sherpa"}
    # sherpa's own THREADS is only the range of its -t option
    assert not any("args.threads" in path.read_text() for path in LEGACY.rglob("*.py"))


def test_1x_still_has_the_methods_the_runner_wraps() -> None:
    assert "\tdef sequence_workflow(self):" in SHERPA
    assert "\tdef html_workflow(self):" in SHERPA
    assert "self.sequence_workflow()" in SHERPA
    assert "self.html_workflow()" in SHERPA
    assert "self.instance_objects.append(current_seqpair)" in SHERPA
