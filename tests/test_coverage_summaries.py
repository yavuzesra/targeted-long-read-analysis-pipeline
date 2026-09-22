from __future__ import annotations

import importlib.util
from pathlib import Path

import pytest


PROJECT_ROOT = Path(__file__).resolve().parents[1]


def load_module(name: str, relative_path: str):
    module_path = PROJECT_ROOT / relative_path

    spec = importlib.util.spec_from_file_location(
        name,
        module_path,
    )

    if spec is None or spec.loader is None:
        raise RuntimeError(
            f"Unable to load module: {module_path}"
        )

    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)

    return module


annotation_module = load_module(
    "prepare_target_annotation",
    "scripts/python/prepare_target_annotation.py",
)

summary_module = load_module(
    "coverage_summaries",
    "scripts/python/coverage_summaries.py",
)


def test_split_multivalue() -> None:
    result = annotation_module.split_multivalue(
        "chr1; chr2;chr3"
    )

    assert result == ["chr1", "chr2", "chr3"]


def test_split_multivalue_rejects_na_as_content() -> None:
    assert annotation_module.split_multivalue("NA") == []
    assert annotation_module.split_multivalue(None) == []


def test_target_key_uses_bed_coordinates() -> None:
    result = summary_module.make_target_key(
        chromosome="chr6",
        start=32840716,
        end=32860734,
    )

    assert result == "chr6:32840716-32860734"


def test_annotation_index_allows_shared_target() -> None:
    rows = [
        {
            "gene_symbol": "PSMB8",
            "chromosome": "chr6",
            "bed_start_0based": "32840716",
            "bed_end_0based_exclusive": "32860734",
            "target_key": "chr6:32840716-32860734",
        },
        {
            "gene_symbol": "PSMB9",
            "chromosome": "chr6",
            "bed_start_0based": "32840716",
            "bed_end_0based_exclusive": "32860734",
            "target_key": "chr6:32840716-32860734",
        },
    ]

    index = summary_module.build_annotation_index(rows)

    assert len(index["chr6:32840716-32860734"]) == 2
    assert {
        row["gene_symbol"]
        for row in index["chr6:32840716-32860734"]
    } == {"PSMB8", "PSMB9"}


def test_annotation_rejects_inconsistent_target_key() -> None:
    rows = [
        {
            "gene_symbol": "GENE1",
            "chromosome": "chr1",
            "bed_start_0based": "100",
            "bed_end_0based_exclusive": "200",
            "target_key": "chr1:100-201",
        }
    ]

    with pytest.raises(
        ValueError,
        match="does not match coordinates",
    ):
        summary_module.build_annotation_index(rows)


def test_region_threshold_detection() -> None:
    rows = [
        {
            "sample_id": "Sample01",
            "chromosome": "chr1",
            "bed_start_0based": "100",
            "bed_end_0based_exclusive": "200",
            "region_length_bp": "100",
            "target_bases": "100",
            "mean_depth": "25",
            "median_depth": "24",
            "minimum_depth": "0",
            "maximum_depth": "50",
            "depth_standard_deviation": "5",
            "depth_coefficient_of_variation": "0.2",
            "zero_depth_bases": "1",
            "zero_depth_percentage": "1",
            "bases_ge_1x": "99",
            "percentage_ge_1x": "99",
            "bases_ge_10x": "90",
            "percentage_ge_10x": "90",
            "bases_ge_20x": "80",
            "percentage_ge_20x": "80",
        }
    ]

    thresholds = summary_module.validate_region_columns(rows)

    assert thresholds == [1, 10, 20]


def test_shared_target_produces_two_gene_rows() -> None:
    region_rows = [
        {
            "sample_id": "Sample01",
            "chromosome": "chr6",
            "bed_start_0based": "32840716",
            "bed_end_0based_exclusive": "32860734",
            "region_length_bp": "20018",
            "target_bases": "20018",
            "mean_depth": "23.1",
            "median_depth": "20",
            "minimum_depth": "0",
            "maximum_depth": "70",
            "depth_standard_deviation": "7",
            "depth_coefficient_of_variation": "0.30",
            "zero_depth_bases": "0",
            "zero_depth_percentage": "0",
            "bases_ge_1x": "20018",
            "percentage_ge_1x": "100",
            "bases_ge_10x": "19000",
            "percentage_ge_10x": "94.91",
            "bases_ge_20x": "12000",
            "percentage_ge_20x": "59.95",
        }
    ]

    annotation_index = {
        "chr6:32840716-32860734": [
            {
                "gene_symbol": "PSMB8",
                "match_status": "Exact full-gene match",
                "coverage_percentage": "100",
            },
            {
                "gene_symbol": "PSMB9",
                "match_status": (
                    "Full gene covered with extra custom panel bases"
                ),
                "coverage_percentage": "100",
            },
        ]
    }

    rows = summary_module.create_gene_rows(
        region_rows=region_rows,
        annotation_index=annotation_index,
        thresholds=[1, 10, 20],
    )

    assert len(rows) == 2
    assert {row["gene_symbol"] for row in rows} == {
        "PSMB8",
        "PSMB9",
    }
    assert all(
        row["percentage_ge_20x"] == pytest.approx(59.95)
        for row in rows
    )
