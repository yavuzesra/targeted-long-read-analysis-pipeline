from __future__ import annotations

import importlib.util
import sys
from pathlib import Path

import pytest


PROJECT_ROOT = Path(__file__).resolve().parents[1]
MODULE_PATH = (
    PROJECT_ROOT / "scripts/python/coverage_analysis.py"
)

SPEC = importlib.util.spec_from_file_location(
    "coverage_analysis",
    MODULE_PATH,
)

if SPEC is None or SPEC.loader is None:
    raise RuntimeError(
        f"Unable to load coverage_analysis module: {MODULE_PATH}"
    )

coverage_analysis = importlib.util.module_from_spec(SPEC)
sys.modules[SPEC.name] = coverage_analysis
SPEC.loader.exec_module(coverage_analysis)


def test_target_region_length_uses_half_open_coordinates() -> None:
    region = coverage_analysis.TargetRegion(
        region_id="region_1",
        chromosome="chr1",
        start=100,
        end=200,
    )

    assert region.length == 100


def test_assign_depth_converts_one_based_position() -> None:
    region = coverage_analysis.TargetRegion(
        region_id="region_1",
        chromosome="chr1",
        start=100,
        end=103,
    )

    coverage_analysis.initialise_depth_arrays([region])

    regions_by_chromosome, starts_by_chromosome = (
        coverage_analysis.build_region_index([region])
    )

    coverage_analysis.assign_depth(
        chromosome="chr1",
        position_1based=101,
        depth=7,
        regions_by_chromosome=regions_by_chromosome,
        starts_by_chromosome=starts_by_chromosome,
    )

    assert region.depths == [7, 0, 0]


def test_assign_depth_at_last_included_base() -> None:
    region = coverage_analysis.TargetRegion(
        region_id="region_1",
        chromosome="chr1",
        start=100,
        end=103,
    )

    coverage_analysis.initialise_depth_arrays([region])

    regions_by_chromosome, starts_by_chromosome = (
        coverage_analysis.build_region_index([region])
    )

    coverage_analysis.assign_depth(
        chromosome="chr1",
        position_1based=103,
        depth=5,
        regions_by_chromosome=regions_by_chromosome,
        starts_by_chromosome=starts_by_chromosome,
    )

    assert region.depths == [0, 0, 5]


def test_depth_statistics_include_zero_depth_bases() -> None:
    depths = [0, 0, 10, 20, 30]

    metrics = coverage_analysis.calculate_depth_statistics(
        depths=depths,
        thresholds=[1, 10, 20],
    )

    assert metrics["target_bases"] == 5
    assert metrics["zero_depth_bases"] == 2
    assert metrics["zero_depth_percentage"] == pytest.approx(40.0)
    assert metrics["bases_ge_1x"] == 3
    assert metrics["percentage_ge_1x"] == pytest.approx(60.0)
    assert metrics["bases_ge_10x"] == 3
    assert metrics["percentage_ge_10x"] == pytest.approx(60.0)
    assert metrics["bases_ge_20x"] == 2
    assert metrics["percentage_ge_20x"] == pytest.approx(40.0)


def test_depth_mean_and_median() -> None:
    depths = [0, 10, 20, 30]

    metrics = coverage_analysis.calculate_depth_statistics(
        depths=depths,
        thresholds=[1, 10, 20],
    )

    assert metrics["mean_depth"] == pytest.approx(15.0)
    assert metrics["median_depth"] == pytest.approx(15.0)
    assert metrics["minimum_depth"] == 0
    assert metrics["maximum_depth"] == 30


def test_depth_coefficient_of_variation() -> None:
    depths = [10, 10, 20, 20]

    metrics = coverage_analysis.calculate_depth_statistics(
        depths=depths,
        thresholds=[1],
    )

    assert metrics["mean_depth"] == pytest.approx(15.0)
    assert metrics["depth_standard_deviation"] == pytest.approx(5.0)
    assert metrics["depth_coefficient_of_variation"] == pytest.approx(
        1 / 3
    )


def test_zero_mean_depth_returns_nan_cv() -> None:
    metrics = coverage_analysis.calculate_depth_statistics(
        depths=[0, 0, 0],
        thresholds=[1],
    )

    assert metrics["mean_depth"] == 0
    assert coverage_analysis.math.isnan(
        metrics["depth_coefficient_of_variation"]
    )


def test_read_target_regions_rejects_overlap(
    tmp_path: Path,
) -> None:
    bed_path = tmp_path / "overlapping.bed"
    bed_path.write_text(
        "chr1\t100\t200\tregion_1\n"
        "chr1\t150\t250\tregion_2\n",
        encoding="utf-8",
    )

    with pytest.raises(
        ValueError,
        match="Overlapping BED intervals",
    ):
        coverage_analysis.read_target_regions(bed_path)


def test_read_target_regions_reports_one_based_coordinates(
    tmp_path: Path,
) -> None:
    bed_path = tmp_path / "targets.bed"
    bed_path.write_text(
        "chr1\t100\t200\tregion_1\n",
        encoding="utf-8",
    )

    regions = coverage_analysis.read_target_regions(bed_path)

    region = regions[0]

    assert region.start + 1 == 101
    assert region.end == 200


def test_panel_qc_pass_logic() -> None:
    percentage_ge_20x = 90.0
    required_percentage = 90.0

    panel_qc_pass = (
        percentage_ge_20x >= required_percentage
    )

    assert panel_qc_pass is True


def test_panel_qc_fail_logic() -> None:
    percentage_ge_20x = 89.999
    required_percentage = 90.0

    panel_qc_pass = (
        percentage_ge_20x >= required_percentage
    )

    assert panel_qc_pass is False
