from __future__ import annotations

import csv
import importlib.util
from pathlib import Path

import pytest


ROOT = Path(__file__).parents[1]


def load_module(name: str, path: Path):
    spec = importlib.util.spec_from_file_location(name, path)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"Unable to load module: {path}")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


metadata = load_module(
    "prepare_giab_sensitivity_metadata",
    ROOT / "scripts/python/prepare_giab_sensitivity_metadata.py",
)
comparison = load_module(
    "compare_giab_truthsets",
    ROOT / "scripts/python/compare_giab_truthsets.py",
)


def make_file(path: Path) -> str:
    path.write_text("test\n", encoding="utf-8")
    return str(path)


def source_row(tmp_path: Path) -> dict[str, str]:
    return {
        "sample_id": "HG002",
        "analysis_role": "validation_control",
        "run_id": "run1",
        "bam_path": "sample.bam",
        "bai_path": "sample.bam.bai",
        "small_variant_vcf": make_file(tmp_path / "query.vcf.gz"),
        "small_variant_vcf_index": make_file(tmp_path / "query.vcf.gz.tbi"),
        "clinvar_vcf": "NA",
        "clinvar_vcf_index": "NA",
        "truth_vcf": "primary.vcf.gz",
        "truth_vcf_index": "primary.vcf.gz.tbi",
        "truth_bed": "primary.bed",
        "truth_version": "v5.0q",
        "sex": "male",
        "include": "true",
        "notes": "primary",
    }


def sensitivity_row(tmp_path: Path) -> dict[str, str]:
    return {
        "benchmark_id": "HG002_v4_2_1",
        "source_sample_id": "HG002",
        "truth_vcf": make_file(tmp_path / "truth.vcf.gz"),
        "truth_vcf_index": make_file(tmp_path / "truth.vcf.gz.tbi"),
        "truth_bed": make_file(tmp_path / "truth.bed"),
        "truth_version": "v4.2.1",
        "include": "true",
        "notes": "alternative truth set",
    }


def test_build_rows_reuses_query_and_changes_truth(tmp_path: Path) -> None:
    source = source_row(tmp_path)
    alternative = sensitivity_row(tmp_path)

    rows = metadata.build_rows([source], [alternative])

    assert len(rows) == 1
    assert rows[0]["sample_id"] == "HG002_v4_2_1"
    assert rows[0]["small_variant_vcf"] == source["small_variant_vcf"]
    assert rows[0]["truth_vcf"] == alternative["truth_vcf"]
    assert rows[0]["truth_version"] == "v4.2.1"
    assert rows[0]["sex"] == "male"


def test_build_rows_requires_explicit_enablement(tmp_path: Path) -> None:
    alternative = sensitivity_row(tmp_path)
    alternative["include"] = "false"

    with pytest.raises(ValueError, match="No sensitivity analyses are enabled"):
        metadata.build_rows([source_row(tmp_path)], [alternative])


def test_build_rows_rejects_duplicate_benchmark_id(tmp_path: Path) -> None:
    alternative = sensitivity_row(tmp_path)

    with pytest.raises(ValueError, match="Duplicate benchmark_id"):
        metadata.build_rows(
            [source_row(tmp_path)], [alternative, alternative.copy()]
        )


def metric_rows(sample_id: str, version: str) -> list[dict[str, str]]:
    return [
        {
            "sample_id": sample_id,
            "truth_version": version,
            "variant_type": "SNP",
            "truth_total": "100",
            "truth_tp": "90",
            "truth_fn": "10",
            "query_total": "96",
            "query_tp": "90",
            "query_fp": "5",
            "query_unk": "1",
            "precision": "0.947368",
            "recall": "0.900000",
            "f1_score": "0.923077",
        },
        {
            "sample_id": sample_id,
            "truth_version": version,
            "variant_type": "INDEL",
            "truth_total": "50",
            "truth_tp": "30",
            "truth_fn": "20",
            "query_total": "42",
            "query_tp": "30",
            "query_fp": "10",
            "query_unk": "2",
            "precision": "0.750000",
            "recall": "0.600000",
            "f1_score": "0.666667",
        },
    ]


def test_output_rows_keeps_truth_and_query_counts_separate() -> None:
    rows = comparison.output_rows(
        "primary",
        "HG002",
        "HG002",
        comparison.metrics_by_sample(metric_rows("HG002", "v5.0q"), "HG002"),
        {"benchmarkable_bases": "1192800"},
    )

    assert [row["variant_type"] for row in rows] == ["SNP", "INDEL"]
    assert rows[0]["truth_tp"] == "90"
    assert rows[0]["truth_fn"] == "10"
    assert rows[0]["query_tp"] == "90"
    assert rows[0]["query_fp"] == "5"
    assert rows[0]["query_unk"] == "1"
    assert rows[0]["precision"] == "0.947368"
    assert rows[0]["recall"] == "0.900000"


def test_manifest_row_requires_unique_match() -> None:
    with pytest.raises(ValueError, match="Expected one manifest row"):
        comparison.manifest_row([], "HG002")
