"""Request and response bodies of the HTTP API.

The web frontend mirrors these in ``apps/web/src/api.ts``.
"""

from __future__ import annotations

from dataclasses import replace
from datetime import datetime
from enum import StrEnum
from typing import Annotated, Any, Literal

from pydantic import (
    AfterValidator,
    BaseModel,
    BeforeValidator,
    ConfigDict,
    Field,
    StringConstraints,
    field_validator,
)
from scalehd.genotype import CallerSettings
from scalehd.pairs import DiscordancePolicy

from .models import Job, JobStatus, SampleStatus, Theme
from .progress import SampleProgress, read_progress
from .workspace import sample_folders


class Health(BaseModel):
    status: Literal["ok"]
    version: str
    core_version: str
    python: str
    platform: str
    sqlite: str
    libraries: dict[str, str]


# Usernames are case-insensitive, so they are stored and compared in lower case.
# Scope of users is like 4 people in a lab so these are not as strict requirements as
# a "real" webserver
_Username = Annotated[
    str, StringConstraints(strip_whitespace=True, max_length=64), AfterValidator(str.lower)
]
_NewPassword = Annotated[str, Field(min_length=8, max_length=256)]


class Login(BaseModel):
    username: _Username
    password: str = Field(max_length=256)


class NewAccount(BaseModel):
    # Letters, digits, dots, dashes and underscores, starting with a letter or digit.
    username: Annotated[_Username, StringConstraints(pattern=r"^[A-Za-z0-9][A-Za-z0-9._-]*$")]
    password: _NewPassword


class PasswordChange(BaseModel):
    current_password: str = Field(max_length=256)
    new_password: _NewPassword


class Folders(BaseModel):
    """Data folders from host machine i.e. sequencing data to pick from, and where results go."""

    data_root: str | None
    workspace: str
    # <workspace>/<username>: this user's jobs are saved here.
    your_folder: str


class Registration(BaseModel):
    """Whether accounts can be created, and whether the next one is the first (the admin)."""

    open: bool
    first_account: bool


class UserOut(BaseModel):
    model_config = ConfigDict(from_attributes=True)

    id: int
    username: str
    is_admin: bool
    created_at: datetime
    theme: Theme


class ThemeChange(BaseModel):
    theme: Theme


class GenotypeMethod(StrEnum):
    """How a job genotypes its samples. The web interface holds the names people see."""

    # ScaleHD 1.x, unchanged. Align reads to the reference library, then the 1.x
    # genotyper. Runs in the Docker image (apps/legacy-runner).
    LEGACY = "legacy"
    # Read the repeat structure straight from each read, then the model-based caller
    # (scalehd.genotype). Work in progress.
    MODEL = "model"


# ScaleHD 1.x's settings: the attributes of its configuration XML
# (its config/config.dtd), element by element, under the XML's own names.
# Defaults are scraped from 1.x's documented configuration.


def _tidy_bases(value: Any) -> Any:
    return value.strip().upper() if isinstance(value, str) else value


Bases = Annotated[str, BeforeValidator(_tidy_bases), StringConstraints(pattern=r"^[ACGT]+$")]
# A bwa mem option taking two whole numbers, e.g. gap open penalties for insertions and
# deletions.
NumberPair = tuple[Annotated[int, Field(ge=0)], Annotated[int, Field(ge=0)]]
# Where a demultiplexing adapter sits: 5' or 3' end. 1.x also takes "AP" (either end), but
# then hands batchadapt no flag for it, and batchadapt writes nothing.
AdapterPosition = Literal["5P", "3P"]
# Error tolerances 1.x refuses although they're between 0 and 1: it checks the value is in
# numpy.arange(0, 1.1, 0.01), and these aren't, exactly.
REFUSED_ERROR_TOLERANCES = frozenset({0.35, 0.41, 0.47, 0.57, 0.69, 0.70, 0.82, 0.83, 0.94, 0.95})


class LegacyInstanceFlags(BaseModel):
    """Which 1.x stages run."""

    demultiplex: bool = False
    quality_control: bool = True
    # Always on. No alignment = no fun and games
    sequence_alignment: Literal[True] = True
    atypical_realignment: bool = True
    genotype_prediction: Literal[True] = True
    snp_calling: bool = True


class LegacyDemultiplexFlags(BaseModel):
    """batchadapt: keeps the reads with these adapters (cutadapt on each read file)."""

    forward_adapter: Bases = "GCGACCCTGG"
    forward_position: AdapterPosition = "5P"
    reverse_adapter: Bases = "GCAGCGGCTG"
    reverse_position: AdapterPosition = "5P"
    error_rate: int = Field(default=0, ge=0)
    min_overlap: int = Field(default=10, ge=0)


class LegacyTrimFlags(BaseModel):
    """cutadapt: quality and adapter trimming."""

    trim_type: Literal["Quality", "Adapter", "Both"] = "Adapter"
    quality_threshold: int = Field(default=5, ge=0, le=38)
    adapter_flag: Literal["-a", "-g", "-a$", "-g^", "-b"] = "-a"
    forward_adapter: Bases = "GATCGGAAGAGCACACGTCTGAACTCCAGTCAC"
    reverse_adapter: Bases = "AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT"
    error_tolerance: float = Field(default=0.39, ge=0, le=1)

    @field_validator("error_tolerance")
    @classmethod
    def _tolerance_1x_takes(cls, value: float) -> float:
        if round(value, 2) != value:
            raise ValueError("legacy ScaleHD takes error tolerances in steps of 0.01")
        if value in REFUSED_ERROR_TOLERANCES:
            raise ValueError(
                f"legacy ScaleHD refuses an error tolerance of {value:.2f}: try one either side"
            )
        return value


class LegacyAlignmentFlags(BaseModel):
    """bwa mem: the option each one sets is in its comment."""

    min_seed_length: int = Field(default=19, ge=0)  # -k
    band_width: int = Field(default=100, ge=0)  # -w
    seed_length_extension: float = Field(default=1.5, ge=0)  # -r
    skip_seed_with_occurrence: int = Field(default=500, ge=0)  # -c
    chain_drop: float = Field(default=0.5, ge=0)  # -D
    seeded_chain_drop: int = Field(default=0, ge=0)  # -W
    seq_match_score: int = Field(default=1, ge=0)  # -A
    mismatch_penalty: int = Field(default=4, ge=0)  # -B
    indel_penalty: NumberPair = (6, 6)  # -O
    gap_extend_penalty: NumberPair = (6, 6)  # -E
    prime_clipping_penalty: NumberPair = (5, 5)  # -L
    unpaired_pairing_penalty: int = Field(default=17, ge=0)  # -U


class LegacyPredictionFlags(BaseModel):
    """SNP calling (freebayes)."""

    # Reads that must show a variant (freebayes -C). 1.x's documentation: 1 to 5.
    snp_observation_threshold: int = Field(default=2, ge=1, le=5)
    # Lowest variant quality (VCF QUAL) kept.
    quality_cutoff: int = Field(default=0, ge=0)


class LegacySettings(BaseModel):
    """How a legacy job runs ScaleHD 1.x."""

    instance_flags: LegacyInstanceFlags = Field(default_factory=LegacyInstanceFlags)
    demultiplex_flags: LegacyDemultiplexFlags = Field(default_factory=LegacyDemultiplexFlags)
    trim_flags: LegacyTrimFlags = Field(default_factory=LegacyTrimFlags)
    alignment_flags: LegacyAlignmentFlags = Field(default_factory=LegacyAlignmentFlags)
    prediction_flags: LegacyPredictionFlags = Field(default_factory=LegacyPredictionFlags)
    # 1.x's -e: keep reads that map equally well to more than one reference.
    enshrine: bool = False
    # 1.x's -b: don't subsample samples with many reads.
    broadscope: bool = False
    # Cores each sample gets (1.x's THREADS). The server runs as many samples at once as
    # its cores allow, so more here means fewer at once.
    cores_per_sample: int = Field(default=1, ge=1)

    def config(self) -> dict[str, dict[str, str]]:
        """The configuration XML's attributes, as the legacy runner takes them."""
        elements = {
            "instance_flags": self.instance_flags,
            "demultiplex_flags": self.demultiplex_flags,
            "trim_flags": self.trim_flags,
            "alignment_flags": self.alignment_flags,
            "prediction_flags": self.prediction_flags,
        }
        config = {
            element: {name: _xml_value(value) for name, value in flags}
            for element, flags in elements.items()
        }
        config["demultiplex_flags"] |= {"min_length": "", "max_length": ""}
        return config


def _xml_value(value: bool | int | float | str | tuple[int, int]) -> str:
    """A setting as 1.x reads it from its XML."""
    if isinstance(value, bool):
        return str(value)
    if isinstance(value, tuple):
        return ",".join(str(n) for n in value)
    return str(value)


class JobSettings(BaseModel):
    """What a job runs, and with which settings. Unset thresholds keep the core defaults."""

    method: GenotypeMethod = GenotypeMethod.MODEL
    # Call genotypes, or only count each sample's repeat structures.
    call: bool = True
    discordant: DiscordancePolicy = DiscordancePolicy.DROP
    # Flag thresholds, as in scalehd.genotype.CallerSettings.
    min_molecules: int | None = Field(default=None, ge=0)
    min_posterior: float | None = Field(default=None, ge=0, le=1)
    max_background: float | None = Field(default=None, ge=0, le=1)
    max_dropped: float | None = Field(default=None, ge=0, le=1)
    # only for legacy jobs....
    legacy: LegacySettings | None = None

    def caller_settings(self) -> CallerSettings:
        thresholds = ("min_molecules", "min_posterior", "max_background", "max_dropped")
        overrides: dict[str, Any] = {
            name: value for name in thresholds if (value := getattr(self, name)) is not None
        }
        return replace(CallerSettings(), **overrides)

    def with_legacy(self) -> JobSettings:
        """These settings with 1.x's filled in (defaults, if unset by user)."""
        return self.model_copy(update={"legacy": self.legacy or LegacySettings()})

    def for_job(self) -> JobSettings:
        """The settings a job keeps"""
        if self.method is GenotypeMethod.LEGACY:
            return self.with_legacy()
        return self.model_copy(update={"legacy": None})


class ServerOptions(BaseModel):
    """What this server can run jobs with."""

    # Legacy needs ScaleHD 1.x, which is only in the Docker image.
    methods: list[GenotypeMethod]
    # Cores for running samples (SCALEHD_WORKERS). A legacy job's samples share them.
    cores: int
    # ScaleHD 1.x's documented settings, for defaulting
    legacy_defaults: LegacySettings = Field(default_factory=LegacySettings)


class InputPair(BaseModel):
    """A sample's FASTQ files, as paths relative to the server's data rpot."""

    name: str
    r1: str
    r2: str | None = None


class InputSample(BaseModel):
    """A sample found in the data folder: its FASTQ files, paired by name."""

    name: str
    files: list[str]
    r1: str | None  # None when it can't be run for any reason
    r2: str | None  # None when it can't be run for any reason
    size: int
    undetermined: bool  # couldn't regex the sample name from filename.
    skipped: str | None  # reason for skip e.g. R2 without R1


class InputSubfolder(BaseModel):
    """A folder inside the open one, with what it holds, for the folder tree."""

    name: str
    folders: int
    samples: int  # runnable samples (name regex OK)


class InputFolder(BaseModel):
    """a run/data dir in the server's data root."""

    folder: str  # Relative to the data folder; "" is the data folder itself.
    path: str  # complete dir path
    folders: list[InputSubfolder]
    samples: list[InputSample]
    other_files: list[str]  # non FASTQ files if any


MAX_TAGS = 5
MAX_TAG_LENGTH = 15
TagName = Annotated[
    str, StringConstraints(strip_whitespace=True, min_length=1, max_length=MAX_TAG_LENGTH)
]


class JobTag(BaseModel):
    """self explanatory"""

    model_config = ConfigDict(from_attributes=True)

    id: int
    name: str


class TagOut(JobTag):
    """A tag, with how many jobs use it (from any user)."""

    jobs: int


class TagChange(BaseModel):
    """Making or renaming a tag."""

    name: TagName


class JobTags(BaseModel):
    """Tags sorted by id."""

    tags: list[int] = Field(max_length=MAX_TAGS)


class JobCreate(BaseModel):
    name: str = Field(min_length=1, max_length=200)
    samples: list[InputPair] = Field(min_length=1)
    settings: JobSettings | None = None  # None = unset = default
    tags: list[int] = Field(default_factory=list, max_length=MAX_TAGS)


class AdminChange(BaseModel):
    is_admin: bool


class PasswordSet(BaseModel):
    """admin reset pwrd"""

    password: _NewPassword


class SampleOut(BaseModel):
    model_config = ConfigDict(from_attributes=True)

    id: int
    name: str
    r1: str | None
    r2: str | None
    status: SampleStatus
    genotype: str | None
    confidence: float | None
    flags: list[str]
    truth: str | None
    matches_truth: bool | None
    error: str | None
    # How far a running sample has got, once its worker has said.
    progress: SampleProgress | None = None


class JobSummary(BaseModel):
    id: int
    name: str
    demo: bool
    method: GenotypeMethod
    status: JobStatus
    created_at: datetime
    started_at: datetime | None
    finished_at: datetime | None
    output_dir: str | None
    sample_count: int
    samples_done: int
    tags: list[JobTag]


class JobOut(JobSummary):
    settings: JobSettings
    samples: list[SampleOut]
    error: str | None


def _summary(job: Job) -> dict[str, Any]:
    done = (SampleStatus.FINISHED, SampleStatus.FAILED)
    return {
        "id": job.id,
        "name": job.name,
        "demo": job.demo,
        "method": JobSettings.model_validate(job.settings).method,
        "status": job.status,
        "created_at": job.created_at,
        "started_at": job.started_at,
        "finished_at": job.finished_at,
        "output_dir": job.output_dir,
        "sample_count": len(job.samples),
        "samples_done": sum(sample.status in done for sample in job.samples),
        "tags": [JobTag.model_validate(tag) for tag in job.tags],
    }


def job_summary(job: Job) -> JobSummary:
    """A job for the jobs list. Call inside an open session."""
    return JobSummary(**_summary(job))


def job_out(job: Job) -> JobOut:
    """A job with its samples, and the progress of those. Call inside an
    open session."""
    samples = [SampleOut.model_validate(sample) for sample in job.samples]
    running = [sample for sample in samples if sample.status is SampleStatus.RUNNING]
    if running:
        folders = sample_folders(job)
        for sample in running:
            sample.progress = read_progress(folders[sample.id])
    return JobOut(
        **_summary(job),
        settings=JobSettings.model_validate(job.settings),
        samples=samples,
        error=job.error,
    )


class CagBar(BaseModel):
    cag: int
    # Molecules at this length, including those read to the tract's end too near the
    # read's own end to confirm it.
    molecules: int
    # Molecules whose CAG is only known to be at least this (reads ended in the tract).
    lower_bound: int


class CagChart(BaseModel):
    """Molecules at each CAG length within one structure of the called alleles."""

    caacag: int
    ccgcca: int
    ccg: int
    cct: int
    # The called alleles on this structure, and their CAG lengths (to highlight).
    alleles: list[str]
    called: list[int]
    bars: list[CagBar]


class CcgBar(BaseModel):
    ccg: int
    molecules: int


class Cell(BaseModel):
    cag: int
    ccg: int
    molecules: int


class Reads(BaseModel):
    """Sample read information"""

    molecules: int
    complete: int
    partial: int
    dropped: int
    unusable: int
    read_outcomes: dict[str, int]
    discordant: dict[str, int]


class LegacyCagChart(BaseModel):
    """Aligned reads at each CAG length (1-200) at one CCG: what ScaleHD 1.x charted."""

    # The called alleles it shows (two alleles with the same CCG share one), and their CAG.
    alleles: list[str]
    called: list[int]
    ccg: int
    reads: list[int]


class LegacyCcgChart(BaseModel):
    """Aligned reverse reads at each CCG length (1-20), as from 1.x."""

    alleles: list[str]
    reads: list[int]


class LegacySample(BaseModel):
    """What ScaleHD 1.x made of a legacy job's sample (its result.json)."""

    status: Literal["complete", "failed"]
    # 1.x's name for the stage it failed in (e.g. "DSP"), and why?
    failed_stage: str | None
    error: str | None
    # The 1.x stages it went through.
    stages: list[str]
    seconds: float | None
    threads: int | None
    versions: dict[str, Any] | None
    # The call (scalehd.legacy-call/1), if 1.x got as far as genotyping.
    call: dict[str, Any] | None
    cag_charts: list[LegacyCagChart]
    ccg_charts: list[LegacyCcgChart]


class SampleDetail(BaseModel):
    """All data needed in an indv sample result page"""

    sample: SampleOut
    job_id: int
    job_name: str
    demo: bool
    tags: list[JobTag]
    folder: str | None
    # The full call (scalehd.call/2), once the sample has been called.
    call: dict[str, Any] | None
    cag_charts: list[CagChart]
    ccg: list[CcgBar]
    cells: list[Cell]
    reads: Reads | None
    legacy: LegacySample | None = None
    # Which files can be downloaded (results.FILE_KINDS).
    files: list[str]
    previous_id: int | None
    next_id: int | None
