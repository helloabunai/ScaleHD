#!/opt/legacy/python/bin/python -I
"""Run ScaleHD 1.x on one sample, as 1.x itself would, and save output for parsing by new-app.

    scalehd-legacy-run REQUEST.json OUT_DIR

REQUEST.json names the sample and how to run it (i.e. old XML config but in JSON)

    {
      "label": "sample-1",
      "r1": "/data/.../reads_R1.fastq.gz",
      "r2": "/data/.../reads_R2.fastq.gz",
      "threads": 4,
      "enshrine": false,
      "broadscope": false,
      "config": {
        "instance_flags": {"demultiplex": "False", "quality_control": "True", ...},
        "demultiplex_flags": {...},
        "trim_flags": {...},
        "alignment_flags": {...},
        "prediction_flags": {...}
      }
    }

With demultiplex "True", 1.x demultiplexes the sample itself first (batchadapt, which runs
cutadapt on each read file on its own). Same behaviour as CLI old version.

OUT_DIR must be empty or not exist yet. It ends up containing:

    config.xml    the configuration 1.x was given
    run/          1.x's own output = InstanceReport.csv, the sample's folder, ScaleHDLog.txt, etc
    result.json   distilled 1.x results for use in reporting to web-app (todo)

While it runs, a line of JSON goes to stdout as each 1.x stage starts, {"stage": "alignment"}.
1.x's log and the output of the programs it runs go to stderr.

1.x itself is not changed for scientific reproducability for those who still use it.
Around the code, this only:
- sets the THREADS that three 1.x modules copy at import to the request's threads, where
  1.x would use every core on the machine
- skips the HTML report 1.x writes at the end (cos got web frontend now)
- listens to 1.x's log for each stage starting and for failure messages (for stage reporting
  on frontend)
- reads every get_* value from the sample and its two alleles afterwards

Runs on 1.x's Python (3.8) as that is what was used in MVLS-RSE repository.
"""

from __future__ import annotations

import argparse
import importlib
import inspect
import json
import logging
import math
import os
import platform
import re
import shutil
import sys
import time
import traceback
from pathlib import Path
from typing import Any, TextIO
from xml.sax.saxutils import quoteattr

SCHEMA = "scalehd.legacy-result/1"
# /opt/legacy in the image. bin/ (this), python/, tools/, references/.
LEGACY_ROOT = Path(__file__).resolve().parent.parent
FORWARD_REFERENCE = "4k-HD-INTER.fasta"
REVERSE_REFERENCE = "20-TypicalReverse.fasta"
# 1.x's -j: its output goes to OUT_DIR/run.
JOB_NAME = "run"

# XML config elements
ELEMENTS = (
    "instance_flags",
    "demultiplex_flags",
    "trim_flags",
    "alignment_flags",
    "prediction_flags",
)
# Sample labels and paths but safe
LABEL = re.compile(r"[A-Za-z0-9][A-Za-z0-9_-]*")
SHELL_SAFE = re.compile(r"[A-Za-z0-9_./+:=@%,-]+")
THREAD_MODULES = (
    "ScaleHD.seq_qc.__quality_control",
    "ScaleHD.align.__alignment",
    "ScaleHD.align.__atypical",
)
# Stage equivalents for reporting back to web-app (eventually).
STAGES = (
    ("Indexing reference(s)", "indexing"),
    ("Demultiplexing reads", "demultiplexing"),
    ("Executing sequence quality control workflow", "quality_control"),
    ("Executing alignment workflow", "alignment"),
    ("Scanning for atypical alleles", "atypical_scanning"),
    ("User specified sequence re-alignment", "realignment"),
    ("Genotyping alleles", "genotyping"),
    ("Calling SNPs", "snp_calling"),
)
_COLOUR = re.compile(r"\x1b\[[0-9;]*m")
_PREFIX = "shd__"
# 1.x's "<stage> failure on <sample>: <error>". Its sample is DMPX_<label> once demultiplexed.
_FAILURE = re.compile(r"(.+?) failure on \S+: (.*)", re.DOTALL)


class RequestError(ValueError):
    pass


def read_request(data: Any) -> dict[str, Any]:
    """The request, checked, with the defaults filled in."""
    if not isinstance(data, dict):
        raise RequestError("the request must be a JSON object")
    label = data.get("label")
    if not isinstance(label, str) or not LABEL.fullmatch(label):
        raise RequestError(
            "label must be letters, digits, '_' and '-', starting with a letter or digit"
        )
    references = {
        "forward_reference": str(LEGACY_ROOT / "references" / FORWARD_REFERENCE),
        "reverse_reference": str(LEGACY_ROOT / "references" / REVERSE_REFERENCE),
    }
    files = {}
    for name in ("r1", "r2", *references):
        path = data.get(name, references.get(name))
        if not isinstance(path, str) or not Path(path).is_file():
            raise RequestError(f"{name} must be the path of a file: {path!r}")
        files[name] = str(Path(path).resolve())
    for name in references:
        if not SHELL_SAFE.fullmatch(Path(files[name]).name):
            raise RequestError(f"{name}'s file name has characters 1.x can't use: {files[name]}")
    threads = data.get("threads")
    if isinstance(threads, bool) or not isinstance(threads, int) or threads < 1:
        raise RequestError("threads must be a whole number, at least 1")
    flags = {}
    for name in ("enshrine", "broadscope"):
        flags[name] = data.get(name, False)
        if not isinstance(flags[name], bool):
            raise RequestError(f"{name} must be true or false")
    return {
        "label": label,
        **files,
        "threads": threads,
        **flags,
        "config": _read_config(data.get("config")),
    }


def _read_config(config: Any) -> dict[str, dict[str, str]]:
    if not isinstance(config, dict) or set(config) != set(ELEMENTS):
        raise RequestError("config must have exactly these elements: " + ", ".join(ELEMENTS))
    for element, attributes in config.items():
        if not isinstance(attributes, dict) or not all(
            isinstance(name, str) and isinstance(value, str) for name, value in attributes.items()
        ):
            raise RequestError(f"config's {element} must map attribute names to strings")
    return {element: dict(config[element]) for element in ELEMENTS}


def prepare_out_dir(path: Path) -> Path:
    out_dir = path.resolve()
    if not SHELL_SAFE.fullmatch(str(out_dir)):
        raise RequestError(f"the output folder's path has characters 1.x can't use: {out_dir}")
    if out_dir.exists() and (not out_dir.is_dir() or any(out_dir.iterdir())):
        raise RequestError(f"the output folder must be empty: {out_dir}")
    out_dir.mkdir(parents=True, exist_ok=True)
    return out_dir


def stage_reads(label: str, r1: str, r2: str, data_dir: Path) -> list[Path]:
    """Copy the two read files to data_dir, named as 1.x pairs them (<label>_R1, <label>_R2).

    A copy, because 1.x decompresses everything in its data folder where it is sourced.
    """
    data_dir.mkdir(parents=True)
    staged = []
    for read, source in (("R1", r1), ("R2", r2)):
        suffix = ".fastq.gz" if _is_gzip(source) else ".fastq"
        target = data_dir / f"{label}_{read}{suffix}"
        shutil.copyfile(source, str(target))
        staged.append(target)
    return staged


def _is_gzip(path: str) -> bool:
    with open(path, "rb") as file:
        return file.read(2) == b"\x1f\x8b"


def config_xml(config: dict[str, dict[str, str]], data_dir: str, forward: str, reverse: str) -> str:
    paths = {"data_dir": data_dir, "forward_reference": forward, "reverse_reference": reverse}
    lines = [f"<config{_attributes(paths)}>"]
    lines += [f"\t<{element}{_attributes(config[element])}/>" for element in ELEMENTS]
    lines.append("</config>")
    return "\n".join(lines) + "\n"


def _attributes(values: dict[str, str]) -> str:
    return "".join(f" {name}={quoteattr(value)}" for name, value in values.items())


def legacy_environment(
    environ: dict[str, str], threads: int, root: Path = LEGACY_ROOT
) -> dict[str, str]:
    """The environment 1.x runs in, binaries on PATH, etc"""
    env = dict(environ)
    paths = [str(root / "bin"), str(root / "python" / "bin"), str(root / "tools" / "bin")]
    if environ.get("PATH"):
        paths.append(environ["PATH"])
    env["PATH"] = os.pathsep.join(paths)
    for name in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
        env[name] = str(threads)
    env["MPLBACKEND"] = "Agg"
    # The 1.x programs written in Python (cutadapt, generatr) must not see another Python's.
    for name in ("PYTHONPATH", "PYTHONHOME"):
        env.pop(name, None)
    return env


class Listener(logging.Handler):
    """Follows 1.x's log to stderr, reports each stage as it starts, and keeps failure messages."""

    def __init__(self, events: TextIO | None) -> None:
        super().__init__(logging.DEBUG)
        self.events = events
        self.stages: list[str] = []
        self.errors: list[str] = []
        self.failure: str | None = None

    def emit(self, record: logging.LogRecord) -> None:
        try:
            self.heard(record.levelno, record.getMessage())
        except Exception:
            self.handleError(record)

    def heard(self, level: int, text: str) -> None:
        text = _COLOUR.sub("", text)
        ours = _PREFIX in text
        message = text.split(_PREFIX, 1)[1].strip() if ours else text.strip()
        if level >= logging.ERROR:
            # Including 1.x's DTD check, which logs without the prefix.
            self.errors.append(message)
        if not ours:
            # Not from 1.x itself (matplotlib and fontTools log to same stream).
            return
        for start, stage in STAGES:
            if message.startswith(start):
                self.stage(stage)
        failed = _FAILURE.fullmatch(message)
        if failed and self.failure is None:
            self.failure = failed.group(2).strip()

    def stage(self, stage: str) -> None:
        self.stages.append(stage)
        if self.events is not None:
            self.events.write(json.dumps({"stage": stage}) + "\n")
            self.events.flush()


def run_scalehd(request: dict[str, Any], out_dir: Path, listener: Listener) -> Any:
    """Run 1.x on out_dir/config.xml. Returns the 1.x instance, or None if 1.x stopped
    before it got to the sample."""
    original_basic_config = logging.basicConfig

    def basic_config(*args: Any, **kwargs: Any) -> None:
        # 1.x sets up logging first thing, listen in once it has.
        original_basic_config(*args, **kwargs)
        logging.getLogger().addHandler(listener)

    logging.basicConfig = basic_config

    sherpa: Any = importlib.import_module("ScaleHD.sherpa")
    for name in THREAD_MODULES:
        module: Any = importlib.import_module(name)
        if not hasattr(module, "THREADS"):
            raise RuntimeError(f"{name} has no THREADS to set")
        module.THREADS = str(request["threads"])

    instances = []
    sequence_workflow = sherpa.ScaleHD.sequence_workflow

    def keep_instance(self: Any) -> None:
        instances.append(self)
        sequence_workflow(self)

    sherpa.ScaleHD.sequence_workflow = keep_instance
    sherpa.ScaleHD.html_workflow = lambda self: None

    sys.argv = ["ScaleHD", "-v", "-c", str(out_dir / "config.xml"), "-o", str(out_dir)]
    sys.argv += ["-j", JOB_NAME]
    if request["enshrine"]:
        sys.argv.append("-e")
    if request["broadscope"]:
        sys.argv.append("-b")
    try:
        sherpa.ScaleHD()
    except SystemExit:
        # 1.x logged why before it exited
        pass
    except Exception as error:
        traceback.print_exc()
        listener.errors.append(f"1.x stopped: {type(error).__name__}: {error}")
    return instances[0] if instances else None


def getter_values(thing: Any) -> dict[str, Any]:
    """Every get_* value of a 1.x object, keyed by the getter's name without "get_"."""
    values = {}
    for name in sorted(dir(thing)):
        getter = getattr(thing, name)
        if name.startswith("get_") and callable(getter) and _takes_nothing(getter):
            values[name[len("get_") :]] = to_json(getter())
    return values


def _takes_nothing(function: Any) -> bool:
    try:
        return not inspect.signature(function).parameters
    except (TypeError, ValueError):
        return False


def to_json(value: Any) -> Any:
    if value is None or isinstance(value, (bool, str, int)):
        return value
    if isinstance(value, float):
        return value if math.isfinite(value) else repr(float(value))
    if hasattr(value, "tolist"):
        return to_json(value.tolist())
    if isinstance(value, dict):
        return {str(key): to_json(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [to_json(item) for item in value]
    if type(value).__module__.startswith("ScaleHD."):
        return getter_values(value)
    return {"type": type(value).__name__, "repr": repr(value)}


def outcome(instance: Any, listener: Listener) -> dict[str, Any]:
    """The sample's status, the stage it failed at (1.x's name for it) and any given reason"""
    samples = getattr(instance, "instance_objects", None) or []
    if not samples:
        error = "; ".join(listener.errors) or "1.x stopped before the sample"
        return {"status": "failed", "failed_stage": None, "error": error, "sample": None}
    sample = samples[0]
    stage = sample.get_exceptionraised()
    if stage == "N/A":
        return {"status": "complete", "failed_stage": None, "error": None, "sample": sample}
    return {"status": "failed", "failed_stage": stage, "error": listener.failure, "sample": sample}


def run(request: dict[str, Any], out_dir: Path, listener: Listener) -> dict[str, Any]:
    listener.stage("preparing")
    data_dir = out_dir / "input"
    stage_reads(request["label"], request["r1"], request["r2"], data_dir)
    config = config_xml(
        request["config"],
        str(data_dir),
        request["forward_reference"],
        request["reverse_reference"],
    )
    (out_dir / "config.xml").write_text(config)
    # Where 1.x's demultiplexing puts the sample's reads.
    demultiplexed = out_dir / "input_demultiplexed"
    try:
        instance = run_scalehd(request, out_dir, listener)
        if request["config"]["instance_flags"].get("demultiplex") == "True" and not any(
            demultiplexed.glob("*")
        ):
            listener.errors.append("demultiplexing (batchadapt) wrote no reads")
    finally:
        # Only the copies made above and by demultiplexing (1.x keeps its own trimmed reads).
        shutil.rmtree(str(data_dir), ignore_errors=True)
        shutil.rmtree(str(demultiplexed), ignore_errors=True)
    found = outcome(instance, listener)
    found["sample"] = to_json(found["sample"])
    found["instance_report"] = _text(out_dir / JOB_NAME / "InstanceReport.csv")
    found["aligned_distributions"] = _text(out_dir / JOB_NAME / "AlignedDistributions.csv")
    return found


def _text(path: Path) -> str | None:
    return path.read_text() if path.is_file() else None


def versions(root: Path = LEGACY_ROOT) -> dict[str, Any]:
    from importlib.metadata import PackageNotFoundError, version

    try:
        scalehd: str | None = version("ScaleHD")
    except PackageNotFoundError:
        scalehd = None
    upstream = {}
    record = root / "UPSTREAM"
    if record.is_file():
        for line in record.read_text().splitlines():
            name, _, value = line.partition("=")
            upstream[name] = value
    return {"scalehd": scalehd, "upstream": upstream or None, "python": platform.python_version()}


def write_json(path: Path, data: dict[str, Any]) -> None:
    aside = path.with_name(f".{path.name}.part")
    aside.write_text(json.dumps(data, indent=2, allow_nan=False) + "\n")
    os.replace(str(aside), str(path))


def _take_stdout() -> TextIO:
    """Keep stdout for stage lines, and send everything else written there to stderr."""
    sys.stdout.flush()
    events = os.fdopen(os.dup(1), "w")
    os.dup2(2, 1)
    devnull = os.open(os.devnull, os.O_RDONLY)
    os.dup2(devnull, 0)
    os.close(devnull)
    return events


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(
        prog="scalehd-legacy-run", description="Run ScaleHD 1.x on one sample."
    )
    parser.add_argument("request", type=Path, help="the request (JSON)")
    parser.add_argument("out_dir", type=Path, help="an empty folder for the output")
    args = parser.parse_args(argv)
    try:
        request = read_request(json.loads(args.request.read_text()))
        out_dir = prepare_out_dir(args.out_dir)
    except (OSError, ValueError) as error:
        print(f"scalehd-legacy-run: {error}", file=sys.stderr)
        return 2

    environment = legacy_environment(dict(os.environ), request["threads"])
    os.environ.clear()
    os.environ.update(environment)
    listener = Listener(_take_stdout())
    started = time.monotonic()
    try:
        found = run(request, out_dir, listener)
    except Exception as error:
        traceback.print_exc()
        found = {
            "status": "failed",
            "failed_stage": None,
            "error": f"the runner failed: {type(error).__name__}: {error}",
            "sample": None,
            "instance_report": None,
            "aligned_distributions": None,
        }
    result = {
        "schema": SCHEMA,
        "label": request["label"],
        "status": found["status"],
        "failed_stage": found["failed_stage"],
        "error": found["error"],
        "threads": request["threads"],
        "seconds": round(time.monotonic() - started, 1),
        "versions": versions(),
        "stages": listener.stages,
        "sample": found["sample"],
        "instance_report": found["instance_report"],
        "aligned_distributions": found["aligned_distributions"],
    }
    write_json(out_dir / "result.json", result)
    return 0 if result["status"] == "complete" else 1


if __name__ == "__main__":
    sys.exit(main())
