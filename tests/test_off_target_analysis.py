from __future__ import annotations

import importlib.util
from pathlib import Path

import pytest


PROJECT_ROOT = Path(__file__).resolve().parents[1]
SCRIPT_PATH = PROJECT_ROOT / "scripts" / "python" / "off_target_analysis.py"

SPEC = importlib.util.spec_from_file_location(
    "off_target_analysis",
    SCRIPT_PATH,
)

if SPEC is None or SPEC.loader is None:
    raise RuntimeError(f"Could not load {SCRIPT_PATH}")

off_target = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(off_target)


def test_base_summary_is_consistent() -> None:
    result = off_target.summarise_base_counts(
        total_retained_aligned_bases=1_000,
        off_target_aligned_bases=750,
    )

    assert result["total_retained_aligned_bases"] == 1_000
    assert result["on_target_aligned_bases"] == 250
    assert result["off_target_aligned_bases"] == 750
    assert result["on_target_fraction_percent"] == pytest.approx(25.0)
    assert result["off_target_fraction_percent"] == pytest.approx(75.0)

    assert (
        result["on_target_aligned_bases"]
        + result["off_target_aligned_bases"]
        == result["total_retained_aligned_bases"]
    )

    assert (
        result["on_target_fraction_percent"]
        + result["off_target_fraction_percent"]
        == pytest.approx(100.0)
    )


def test_zero_total_returns_zero_percentages() -> None:
    result = off_target.summarise_base_counts(0, 0)

    assert result["on_target_aligned_bases"] == 0
    assert result["on_target_fraction_percent"] == 0.0
    assert result["off_target_fraction_percent"] == 0.0


def test_off_target_bases_cannot_exceed_total() -> None:
    with pytest.raises(
        ValueError,
        match="exceed total retained",
    ):
        off_target.summarise_base_counts(100, 101)


def test_module_does_not_create_filtered_bam() -> None:
    source = SCRIPT_PATH.read_text(encoding="utf-8")

    assert ".primary.bam" not in source
    assert "samtools view -u" not in source
    assert "run_pipeline_to_file" in source
    assert '"-u"' in source
