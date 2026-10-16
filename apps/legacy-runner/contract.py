"""Checks to ensure legacy ScaleHD is what is expected.

Not duplicating 1.x code in this repo for maintenance ease reasons.
The image fetches the commit pinned in apps/legacy-runner/UPSTREAM into /opt/legacy/upstream.
tools/legacy-parity.sh runs these checks on that location within docker.

If you want to manually run a check outside the image, point it at a checkout of 1.x:

    python apps/legacy-runner/contract.py [UPSTREAM_DIR]
"""

from __future__ import annotations

import re
import sys
import xml.etree.ElementTree as ET
from collections.abc import Callable
from pathlib import Path

import scalehd_legacy_run as runner
from scalehd_server import legacy
from scalehd_server.schemas import LegacySettings

UPSTREAM = Path("/opt/legacy/upstream")


class Source:
    """dirs/files to read"""

    def __init__(self, root: Path) -> None:
        self.package = root / "src" / "ScaleHD"
        self.docs = root / "src" / "docs"
        self.sherpa = (self.package / "sherpa.py").read_text()
        self.dtd = (self.package / "config" / "config.dtd").read_text()


def stages_logged(source: Source) -> list[str]:
    return [
        f"no longer logs {start!r}"
        for start, _ in runner.STAGES
        if f"'{start}" not in source.sherpa
    ]


def failure_messages(source: Source) -> list[str]:
    lines = [line for line in source.sherpa.splitlines() if "failure on " in line]
    problems = [] if len(lines) == 7 else [f"{len(lines)} failure messages, not 7"]
    pattern = re.compile(r"failure on ',\s*seqpair_lbl,\s*str\(e\)\)")
    return problems + [
        f"a failure message changed: {line.strip()}"
        for line in lines
        if ".format(clr.red" not in line or not pattern.search(line)
    ]


def threads(source: Source) -> list[str]:
    files = sorted(source.package.rglob("*.py"))
    modules = {
        "ScaleHD." + ".".join(path.relative_to(source.package).with_suffix("").parts)
        for path in files
        if re.search(r"\bTHREADS\b|cpu_count", path.read_text())
    }
    # sherpa's own THREADS is only the range of its -t option
    expected = {*runner.THREAD_MODULES, "ScaleHD.sherpa"}
    problems = [] if modules == expected else [f"THREADS or cpu_count in {sorted(modules)}"]
    return problems + [
        f"{path.name} reads args.threads" for path in files if "args.threads" in path.read_text()
    ]


def wrapped_methods(source: Source) -> list[str]:
    """1.x has the methods the runner wraps, and keeps its samples where the runner reads them"""
    needed = (
        "\tdef sequence_workflow(self):",
        "\tdef html_workflow(self):",
        "self.sequence_workflow()",
        "self.html_workflow()",
        "self.instance_objects.append(current_seqpair)",
    )
    return [
        f"sherpa.py no longer has {text.strip()!r}" for text in needed if text not in source.sherpa
    ]


def element_order(source: Source) -> list[str]:
    """The runner writes the XML's elements in the DTD's order"""
    declared = re.search(r"<!ELEMENT config\s+\(([^)]*)\)>", source.dtd)
    found = tuple(name.strip() for name in declared.group(1).split(",")) if declared else ()
    return [] if found == runner.ELEMENTS else [f"the DTD's elements are {found}"]


def required_attributes(source: Source) -> list[str]:
    """The server's settings have every attribute the DTD requires"""
    config = LegacySettings().config()
    required = re.findall(r"<!ATTLIST (\w+)\s+(\w+) CDATA #REQUIRED>", source.dtd)
    return [
        f"no {element} {name}"
        for element, name in required
        if element != "config" and name not in config.get(element, {})
    ]


def documented_defaults(source: Source) -> list[str]:
    """The server's defaults are 1.x's documented configuration"""
    text = (source.docs / "ProvidingInput.rst").read_text()
    root = ET.fromstring(text[text.index("<config") : text.index("</config>") + len("</config>")])
    documented = {child.tag: dict(child.attrib) for child in root}
    config = LegacySettings().config()
    problems = [] if list(config) == list(documented) else [f"elements {list(documented)}"]
    for element, attributes in documented.items():
        ours = config.get(element, {})
        if ours.keys() != attributes.keys():
            problems.append(f"{element}'s attributes are {list(attributes)}")
        problems += [
            f"{element} {name} is {ours.get(name)!r}, documented {value!r}"
            for name, value in attributes.items()
            if not _same(ours.get(name, ""), value)
        ]
    return problems


def _same(written: str, documented: str) -> bool:
    try:
        return float(written) == float(documented)
    except ValueError:
        return written == documented


def flag_columns(source: Source) -> list[str]:
    """The server's flags are 1.x's InstanceReport columns, in order"""
    header = source.sherpa[
        source.sherpa.index("'Exception Raised'") : source.sherpa.index("'Low Peak Reads'")
    ]
    columns = [*re.findall(r"'([^']+)'", header)[1:], "Low Peak Reads"]
    ours = [column for column, _ in legacy.FLAGS]
    return [] if ours == columns else [f"1.x's flag columns are {columns}"]


def flag_values(source: Source) -> list[str]:
    """The server reads each flag from the value 1.x's report does"""
    listed = source.sherpa[
        source.sherpa.index("'get_exceptionraised'") : source.sherpa.index("report_string =")
    ]
    getters = re.findall(r"'get_(\w+)'", listed)[1:]
    ours = [key for _, key in legacy.FLAGS]
    return [] if ours == getters else [f"1.x's report reads {getters}"]


def failed_stages(source: Source) -> list[str]:
    """Every stage 1.x can fail in has words on the server"""
    stages = set(re.findall(r"set_exceptionraised\('([^']+)'\)", source.sherpa)) - {"N/A"}
    return [] if stages == set(legacy.FAILED_STAGES) else [f"1.x's stages are {sorted(stages)}"]


CHECKS: tuple[Callable[[Source], list[str]], ...] = (
    stages_logged,
    failure_messages,
    threads,
    wrapped_methods,
    element_order,
    required_attributes,
    documented_defaults,
    flag_columns,
    flag_values,
    failed_stages,
)


def check(root: Path = UPSTREAM) -> list[str]:
    """Run every check, printing how each went. Returns what was wrong."""
    source = Source(root)
    failures = []
    for each in CHECKS:
        try:
            found = each(source)
        except (OSError, ValueError) as error:
            found = [f"couldn't check: {error}"]
        print(f"{each.__doc__:<92} {'ok' if not found else 'DIFFERS'}")
        failures += [f"{each.__name__}: {problem}" for problem in found]
    return failures


if __name__ == "__main__":
    failures = check(Path(sys.argv[1]) if len(sys.argv) > 1 else UPSTREAM)
    for failure in failures:
        print(failure)
    sys.exit(1 if failures else 0)
