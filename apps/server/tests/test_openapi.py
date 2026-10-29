"""automated api generation checks"""

import json
from typing import Any

import pytest
from scalehd.genotype import Flag
from scalehd_server.app import create_app
from scalehd_server.config import ServerSettings


def _openapi() -> dict[str, Any]:
    return create_app(ServerSettings(web_dir=None)).openapi()


def _ref(name: str) -> dict[str, str]:
    return {"$ref": f"#/components/schemas/{name}"}


def test_sample_detail_describes_the_genotype_call() -> None:
    schemas = _openapi()["components"]["schemas"]
    assert _ref("GenotypeCallJson") in schemas["SampleDetail"]["properties"]["call"]["anyOf"]
    call = schemas["GenotypeCallJson"]["properties"]
    assert call["alleles"]["items"] == _ref("AlleleCallJson")
    assert call["flags"]["items"] == _ref("Flag")
    assert schemas["AlleleCallJson"]["properties"]["stutter"] == _ref("StutterJson")


def test_sample_flags_are_the_callers_flags() -> None:
    schemas = _openapi()["components"]["schemas"]
    assert schemas["SampleOut"]["properties"]["flags"]["items"] == _ref("Flag")
    assert schemas["Flag"]["enum"] == [flag.value for flag in Flag]


def test_sample_files_are_the_downloadable_kinds() -> None:
    files = _openapi()["components"]["schemas"]["SampleDetail"]["properties"]["files"]
    assert files["items"]["enum"] == ["call", "counts", "r1", "r2"]


def test_openapi_command_prints_the_apps_description(capsys: pytest.CaptureFixture[str]) -> None:
    from scalehd_server.openapi import main

    main()
    assert json.loads(capsys.readouterr().out) == _openapi()
