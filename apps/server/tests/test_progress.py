"""test sample generated progress json dumps"""

from pathlib import Path

from scalehd_server.progress import (
    ProgressReporter,
    SampleProgress,
    SampleStage,
    clear_progress,
    read_progress,
)


def test_a_new_stage_is_written_at_once_and_shares_at_most_every_so_often(tmp_path: Path) -> None:
    folder = tmp_path / "s1"
    progress = ProgressReporter(folder, every=60)
    progress.stage(SampleStage.COUNTING, 0.0)
    progress.counted(0.5)
    assert read_progress(folder) == SampleProgress(stage=SampleStage.COUNTING, share=0.0)
    progress.stage(SampleStage.CALLING)
    assert read_progress(folder) == SampleProgress(stage=SampleStage.CALLING)
    assert [path.name for path in folder.iterdir()] == ["progress.json"]


def test_shares_are_written_once_enough_time_has_passed(tmp_path: Path) -> None:
    progress = ProgressReporter(tmp_path, every=0)
    progress.counted(0.1)
    progress.counted(0.6)
    assert read_progress(tmp_path) == SampleProgress(stage=SampleStage.COUNTING, share=0.6)


def test_no_progress_or_an_unreadable_file_reads_as_none(tmp_path: Path) -> None:
    assert read_progress(tmp_path) is None
    (tmp_path / "progress.json").write_text('{"stage": "count')
    assert read_progress(tmp_path) is None
    clear_progress(tmp_path)
    clear_progress(tmp_path)
    assert list(tmp_path.iterdir()) == []
