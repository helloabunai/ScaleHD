import pytest
import subprocess


@pytest.mark.integration
def test_can_call():
    subprocess.call([
        "ScaleHD"
    ])
