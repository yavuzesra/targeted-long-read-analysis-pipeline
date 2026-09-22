from __future__ import annotations

import csv
import importlib.util
from pathlib import Path

import pytest


MODULE_PATH = (
    Path(__file__).parents[1]
    / "scripts"
    / "python"
    / "run_variant_application.py"
)

SPEC = importlib.util.spec_from_file_location(
    "run_variant_application",
    MODULE_PATH,
)

if SPEC is None or SPEC.loader is None:
    raise RuntimeError(
        f"Unable to load module: {MODULE_PATH}"
    )

variant_application = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(variant_application)


def write_samples_table(
    path: Path,
    rows: list[dict[str, str]],
) -> None:
    fieldnames = [
        "sample_id",
        "analysis_role",
        "include",
        "small_variant_vcf",
        "small_variant_vcf_index",
        "clinvar_vcf",
        "clinvar_vcf_index",
    ]

    with path.open(
        "w",
        encoding="utf-8",
        newline="",
    ) as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=fieldnames,
            delimiter="\t",
            lineterminator="\n",
        )
        writer.writeheader()
        writer.writerows(rows)


def make_row(
    sample_id: str,
    role: str = "workflow_application",
    include: str = "true",
) -> dict[str, str]:
    return {
        "sample_id": sample_id,
        "analysis_role": role,
        "include": include,
        "small_variant_vcf": "small.vcf.gz",
        "small_variant_vcf_index": "small.vcf.gz.tbi",
        "clinvar_vcf": "clinvar.vcf.gz",
        "clinvar_vcf_index": "clinvar.vcf.gz.tbi",
    }


def test_read_samples_selects_workflow_application(
    tmp_path: Path,
) -> None:
    samples_path = tmp_path / "samples.tsv"

    write_samples_table(
        samples_path,
        [
            make_row("CLINICAL_SAMPLE_01"),
            make_row(
                "HG002",
                role="validation_control",
            ),
            make_row(
                "CLINICAL_SAMPLE_02",
                include="false",
            ),
        ],
    )

    rows = variant_application.read_samples(samples_path)

    assert len(rows) == 1
    assert rows[0]["sample_id"] == "CLINICAL_SAMPLE_01"


def test_read_samples_rejects_missing_columns(
    tmp_path: Path,
) -> None:
    samples_path = tmp_path / "samples.tsv"

    samples_path.write_text(
        "sample_id\tanalysis_role\tinclude\n"
        "CLINICAL_SAMPLE_01\tworkflow_application\ttrue\n",
        encoding="utf-8",
    )

    with pytest.raises(
        ValueError,
        match="missing columns",
    ):
        variant_application.read_samples(samples_path)


def test_read_samples_rejects_empty_selection(
    tmp_path: Path,
) -> None:
    samples_path = tmp_path / "samples.tsv"

    write_samples_table(
        samples_path,
        [
            make_row(
                "HG002",
                role="validation_control",
            ),
        ],
    )

    with pytest.raises(
        ValueError,
        match="No included workflow_application samples",
    ):
        variant_application.read_samples(samples_path)


def test_require_file_accepts_existing_file(
    tmp_path: Path,
) -> None:
    file_path = tmp_path / "input.vcf.gz"
    file_path.write_text("test", encoding="utf-8")

    observed = variant_application.require_file(
        "CLINICAL_SAMPLE_01",
        "small_variant_vcf",
        str(file_path),
    )

    assert observed == file_path


def test_require_file_rejects_missing_value() -> None:
    with pytest.raises(
        ValueError,
        match="small_variant_vcf is missing",
    ):
        variant_application.require_file(
            "CLINICAL_SAMPLE_01",
            "small_variant_vcf",
            "NA",
        )


def test_require_file_rejects_missing_path(
    tmp_path: Path,
) -> None:
    missing_path = tmp_path / "missing.vcf.gz"

    with pytest.raises(
        ValueError,
        match="not found",
    ):
        variant_application.require_file(
            "CLINICAL_SAMPLE_01",
            "small_variant_vcf",
            str(missing_path),
        )
