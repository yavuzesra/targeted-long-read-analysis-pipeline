from __future__ import annotations

import importlib.util
from pathlib import Path

import pytest


PROJECT_ROOT = Path(__file__).resolve().parents[1]
MODULE_PATH = PROJECT_ROOT / "scripts/python/target_metrics.py"

SPEC = importlib.util.spec_from_file_location(
    "target_metrics",
    MODULE_PATH,
)

if SPEC is None or SPEC.loader is None:
    raise RuntimeError(
        f"Unable to load target_metrics module: {MODULE_PATH}"
    )

target_metrics = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(target_metrics)


def test_merge_overlapping_intervals() -> None:
    intervals = [
        ("chr1", 100, 200),
        ("chr1", 150, 250),
        ("chr1", 300, 350),
    ]

    merged = target_metrics.merge_intervals(intervals)

    assert merged == [
        ("chr1", 100, 250),
        ("chr1", 300, 350),
    ]


def test_merge_bookended_intervals() -> None:
    intervals = [
        ("chr1", 100, 200),
        ("chr1", 200, 300),
    ]

    merged = target_metrics.merge_intervals(intervals)

    assert merged == [
        ("chr1", 100, 300),
    ]


def test_intervals_on_different_chromosomes_are_not_merged() -> None:
    intervals = [
        ("chr1", 100, 200),
        ("chr2", 100, 200),
    ]

    merged = target_metrics.merge_intervals(intervals)

    assert merged == [
        ("chr1", 100, 200),
        ("chr2", 100, 200),
    ]


def test_bed_half_open_length_calculation() -> None:
    interval = ("chr1", 100, 200)

    _, start, end = interval
    length = end - start

    assert length == 100


def test_reference_target_fraction() -> None:
    panel_bases = 1_203_874
    reference_bases = 3_099_922_541

    fraction = panel_bases / reference_bases

    assert fraction == pytest.approx(
        0.000388355959,
        rel=1e-6,
    )


def test_alignment_record_metrics() -> None:
    values = target_metrics.calculate_record_values(
        eligible_records=40_000,
        on_target_records=8_000,
    )

    assert values["eligible_mapped_primary_records"] == 40_000
    assert values["on_target_records"] == 8_000
    assert values["off_target_records"] == 32_000
    assert values[
        "on_target_record_percentage"
    ] == pytest.approx(20.0)
    assert values[
        "off_target_record_percentage"
    ] == pytest.approx(80.0)


def test_target_base_enrichment_factor() -> None:
    values = target_metrics.calculate_base_values(
        target_depth_sum=92_220_339,
        total_depth_sum=463_312_433,
        target_length_bp=1_203_874,
        reference_length_bp=3_099_922_541,
    )

    assert values[
        "on_target_aligned_base_fraction_percent"
    ] == pytest.approx(
        19.904568,
        rel=1e-6,
    )

    assert values[
        "mean_target_depth"
    ] == pytest.approx(
        76.602983,
        rel=1e-6,
    )

    assert values[
        "mean_genome_wide_depth"
    ] == pytest.approx(
        0.149459358,
        rel=1e-6,
    )

    assert values[
        "target_base_enrichment_factor"
    ] == pytest.approx(
        512.534,
        rel=1e-5,
    )


def test_run_summary_uses_arithmetic_means() -> None:
    rows = [
        {
            "run_id": "trio_chuv",
            "on_target_record_percentage": 20.0,
            "on_target_aligned_base_fraction_percent": 22.0,
            "target_base_enrichment_factor": 500.0,
        },
        {
            "run_id": "trio_chuv",
            "on_target_record_percentage": 24.0,
            "on_target_aligned_base_fraction_percent": 26.0,
            "target_base_enrichment_factor": 700.0,
        },
    ]

    summaries = target_metrics.summarise_runs(rows)

    assert summaries == [
        {
            "run_id": "trio_chuv",
            "n_samples": 2,
            "mean_on_target_record_percentage": 22.0,
            "mean_on_target_aligned_base_fraction_percent": 24.0,
            "mean_target_base_enrichment_factor": 600.0,
        },
    ]


def test_read_bed_rejects_end_equal_to_start(
    tmp_path: Path,
) -> None:
    bed_path = tmp_path / "invalid.bed"
    bed_path.write_text(
        "chr1\t100\t100\n",
        encoding="utf-8",
    )

    with pytest.raises(
        ValueError,
        match="end > start",
    ):
        target_metrics.read_bed_intervals(bed_path)


def test_validate_bed_rejects_unknown_chromosome() -> None:
    intervals = [
        ("chrUnknown", 0, 100),
    ]

    reference_lengths = {
        "chr1": 1_000,
    }

    with pytest.raises(
        ValueError,
        match="absent from the reference FAI",
    ):
        target_metrics.validate_bed_against_reference(
            intervals,
            reference_lengths,
        )


def test_validate_bed_rejects_interval_beyond_chromosome() -> None:
    intervals = [
        ("chr1", 900, 1_100),
    ]

    reference_lengths = {
        "chr1": 1_000,
    }

    with pytest.raises(
        ValueError,
        match="exceeds chromosome length",
    ):
        target_metrics.validate_bed_against_reference(
            intervals,
            reference_lengths,
        )
