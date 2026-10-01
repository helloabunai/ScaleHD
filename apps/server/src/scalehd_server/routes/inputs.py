"""FASTQ files on the server that jobs can use. Not built yet.

Files are read from the input directory (``SCALEHD_INPUT_DIR``) rather than uploaded,
because a MiSeq run is gigabytes and the server is usually the machine the run was
copied to. Uploads for small one-off samples could come later.
"""

from fastapi import APIRouter

from ..auth import CurrentUser
from ..config import ServerConfig
from ..errors import not_implemented
from ..schemas import InputPair

router = APIRouter(prefix="/inputs", tags=["inputs"])


@router.get("")
def list_inputs(user: CurrentUser, config: ServerConfig, folder: str = "") -> list[InputPair]:
    """FASTQ files in one folder of the input directory, paired into samples.

    TODO: pair ``<sample>_R1*.fastq[.gz]`` with ``<sample>_R2*``, list unpaired files
    as single-end, and refuse paths that leave the input directory.
    """
    raise not_implemented("inputs")
