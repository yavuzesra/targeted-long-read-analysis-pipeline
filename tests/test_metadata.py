from __future__ import annotations

import subprocess
import sys
from pathlib import Path

import yaml


PROJECT_ROOT = Path(__file__).resolve().parents[1]
VALIDATOR = PROJECT_ROOT / "scripts/python/validate_metadata.py"


def write_project_config(tmp_path: Path) -> Path:
    reference = tmp_path / "reference.fna"
    fai = tmp_path / "reference.fna.fai"
    bed = tmp_path / "targets.bed"
    annotation = tmp_path / "annotation.xlsx"

    reference.write_text(">chr1\nACGT\n", encoding="utf-8")
    fai.write_text("chr1\t4\t6\t4\t5\n", encoding="utf-8")
    bed.write_text("chr1\t0\t4\n", encoding="utf-8")
    annotation.write_text("placeholder\n", encoding="utf-8")

    config = {
        "resources": {
            "reference_fasta": str(reference),
            "reference_fai": str(fai),
            "target_bed": str(bed),
            "annotation_xlsx": str(annotation),
        }
    }

    config_path = tmp_path / "project.yaml"

    with config_path.open("w", encoding="utf-8") as handle:
        yaml.safe_dump(config, handle)

    return config_path


def write_runs(tmp_path: Path) -> Path:
    bam_root = tmp_path / "bam_root"
    bam_root.mkdir()

    runs_path = tmp_path / "runs.tsv"
    runs_path.write_text(
        "run_id\tdataset_label\tworkflow\tworkflow_version\t"
        "bam_root\treference_build\tnotes\n"
        f"run1\ttest\twf-alignment\t1.2.5\t"
        f"{bam_root}\tGRCh38\ttest run\n",
        encoding="utf-8",
    )

    return runs_path


def create_alignment_files(tmp_path: Path) -> tuple[Path, Path]:
    bam = tmp_path / "sample.bam"
    bai = tmp_path / "sample.bam.bai"

    bam.write_bytes(b"placeholder")
    bai.write_bytes(b"placeholder")

    return bam, bai


def write_samples(
    tmp_path: Path,
    rows: list[list[str]],
) -> Path:
    header = [
        "sample_id",
        "analysis_role",
        "run_id",
        "bam_path",
        "bai_path",
        "small_variant_vcf",
        "small_variant_vcf_index",
        "clinvar_vcf",
        "clinvar_vcf_index",
        "truth_vcf",
        "truth_vcf_index",
        "truth_bed",
        "truth_version",
        "sex",
        "include",
        "notes",
    ]

    samples_path = tmp_path / "samples.tsv"

    lines = ["\t".join(header)]
    lines.extend("\t".join(row) for row in rows)

    samples_path.write_text(
        "\n".join(lines) + "\n",
        encoding="utf-8",
    )

    return samples_path


def run_validator(
    project_config: Path,
    samples: Path,
    runs: Path,
) -> subprocess.CompletedProcess[str]:
    return subprocess.run(
        [
            sys.executable,
            str(VALIDATOR),
            "--project-config",
            str(project_config),
            "--samples",
            str(samples),
            "--runs",
            str(runs),
        ],
        capture_output=True,
        text=True,
        check=False,
    )


def base_method_row(
    bam: Path,
    bai: Path,
    sample_id: str = "Sample01",
    include: str = "true",
) -> list[str]:
    return [
        sample_id,
        "method_development",
        "run1",
        str(bam),
        str(bai),
        "NA",
        "NA",
        "NA",
        "NA",
        "NA",
        "NA",
        "NA",
        "NA",
        "unknown",
        include,
        "test sample",
    ]


def test_valid_method_development_metadata(tmp_path: Path) -> None:
    project = write_project_config(tmp_path)
    runs = write_runs(tmp_path)
    bam, bai = create_alignment_files(tmp_path)

    samples = write_samples(
        tmp_path,
        [base_method_row(bam, bai)],
    )

    result = run_validator(project, samples, runs)

    assert result.returncode == 0
    assert "Metadata validation completed successfully." in result.stdout
    assert "method_development: 1" in result.stdout


def test_duplicate_sample_id_fails(tmp_path: Path) -> None:
    project = write_project_config(tmp_path)
    runs = write_runs(tmp_path)
    bam, bai = create_alignment_files(tmp_path)

    row = base_method_row(bam, bai)

    samples = write_samples(
        tmp_path,
        [row, row],
    )

    result = run_validator(project, samples, runs)

    assert result.returncode == 1
    assert "duplicate sample_id" in result.stderr


def test_invalid_analysis_role_fails(tmp_path: Path) -> None:
    project = write_project_config(tmp_path)
    runs = write_runs(tmp_path)
    bam, bai = create_alignment_files(tmp_path)

    row = base_method_row(bam, bai)
    row[1] = "patient"

    samples = write_samples(tmp_path, [row])

    result = run_validator(project, samples, runs)

    assert result.returncode == 1
    assert "invalid analysis_role" in result.stderr


def test_validation_control_requires_truth_files(
    tmp_path: Path,
) -> None:
    project = write_project_config(tmp_path)
    runs = write_runs(tmp_path)
    bam, bai = create_alignment_files(tmp_path)

    row = base_method_row(bam, bai, sample_id="HG002")
    row[1] = "validation_control"
    row[5] = "NA"
    row[6] = "NA"

    samples = write_samples(tmp_path, [row])

    result = run_validator(project, samples, runs)

    assert result.returncode == 1
    assert "truth_vcf is required" in result.stderr
    assert "truth_bed is required" in result.stderr
    assert "small_variant_vcf is required" in result.stderr


def test_workflow_application_rejects_truth_metadata(
    tmp_path: Path,
) -> None:
    project = write_project_config(tmp_path)
    runs = write_runs(tmp_path)
    bam, bai = create_alignment_files(tmp_path)

    truth_vcf = tmp_path / "truth.vcf.gz"
    truth_index = tmp_path / "truth.vcf.gz.tbi"
    truth_bed = tmp_path / "truth.bed"

    truth_vcf.write_bytes(b"placeholder")
    truth_index.write_bytes(b"placeholder")
    truth_bed.write_text("chr1\t0\t4\n", encoding="utf-8")

    row = base_method_row(bam, bai, sample_id="CLINICAL_SAMPLE_01")
    row[1] = "workflow_application"
    row[9] = str(truth_vcf)
    row[10] = str(truth_index)
    row[11] = str(truth_bed)
    row[12] = "v1"

    samples = write_samples(tmp_path, [row])

    result = run_validator(project, samples, runs)

    assert result.returncode == 1
    assert "truth-set fields are only allowed" in result.stderr


def test_excluded_sample_is_not_counted(tmp_path: Path) -> None:
    project = write_project_config(tmp_path)
    runs = write_runs(tmp_path)
    bam, bai = create_alignment_files(tmp_path)

    included = base_method_row(
        bam,
        bai,
        sample_id="Sample01",
        include="true",
    )

    excluded = base_method_row(
        bam,
        bai,
        sample_id="Sample02",
        include="false",
    )

    samples = write_samples(
        tmp_path,
        [included, excluded],
    )

    result = run_validator(project, samples, runs)

    assert result.returncode == 0
    assert "method_development: 1" in result.stdout
    assert "total: 1" in result.stdout


def create_variant_files(
    tmp_path: Path,
) -> tuple[Path, Path, Path, Path]:
    small_vcf = tmp_path / "small.vcf.gz"
    small_index = tmp_path / "small.vcf.gz.tbi"
    clinvar_vcf = tmp_path / "clinvar.vcf.gz"
    clinvar_index = tmp_path / "clinvar.vcf.gz.tbi"

    for path in (
        small_vcf,
        small_index,
        clinvar_vcf,
        clinvar_index,
    ):
        path.write_bytes(b"placeholder")

    return (
        small_vcf,
        small_index,
        clinvar_vcf,
        clinvar_index,
    )


def workflow_application_row(
    bam: Path,
    bai: Path,
    small_vcf: Path,
    small_index: Path,
    clinvar_vcf: Path,
    clinvar_index: Path,
    sample_id: str = "CLINICAL_SAMPLE_01",
) -> list[str]:
    return [
        sample_id,
        "workflow_application",
        "run1",
        str(bam),
        str(bai),
        str(small_vcf),
        str(small_index),
        str(clinvar_vcf),
        str(clinvar_index),
        "NA",
        "NA",
        "NA",
        "NA",
        "unknown",
        "true",
        "workflow-application test sample",
    ]


def test_valid_workflow_application_metadata(
    tmp_path: Path,
) -> None:
    project = write_project_config(tmp_path)
    runs = write_runs(tmp_path)
    bam, bai = create_alignment_files(tmp_path)

    (
        small_vcf,
        small_index,
        clinvar_vcf,
        clinvar_index,
    ) = create_variant_files(tmp_path)

    row = workflow_application_row(
        bam,
        bai,
        small_vcf,
        small_index,
        clinvar_vcf,
        clinvar_index,
    )

    samples = write_samples(tmp_path, [row])
    result = run_validator(project, samples, runs)

    assert result.returncode == 0
    assert (
        "Metadata validation completed successfully."
        in result.stdout
    )
    assert "workflow_application: 1" in result.stdout


def test_workflow_application_requires_clinvar_vcf(
    tmp_path: Path,
) -> None:
    project = write_project_config(tmp_path)
    runs = write_runs(tmp_path)
    bam, bai = create_alignment_files(tmp_path)

    (
        small_vcf,
        small_index,
        clinvar_vcf,
        clinvar_index,
    ) = create_variant_files(tmp_path)

    row = workflow_application_row(
        bam,
        bai,
        small_vcf,
        small_index,
        clinvar_vcf,
        clinvar_index,
    )
    row[7] = "NA"

    samples = write_samples(tmp_path, [row])
    result = run_validator(project, samples, runs)

    assert result.returncode == 1
    assert (
        "clinvar_vcf is required for workflow_application samples"
        in result.stderr
    )


def test_workflow_application_requires_clinvar_vcf_index(
    tmp_path: Path,
) -> None:
    project = write_project_config(tmp_path)
    runs = write_runs(tmp_path)
    bam, bai = create_alignment_files(tmp_path)

    (
        small_vcf,
        small_index,
        clinvar_vcf,
        clinvar_index,
    ) = create_variant_files(tmp_path)

    row = workflow_application_row(
        bam,
        bai,
        small_vcf,
        small_index,
        clinvar_vcf,
        clinvar_index,
    )
    row[8] = "NA"

    samples = write_samples(tmp_path, [row])
    result = run_validator(project, samples, runs)

    assert result.returncode == 1
    assert (
        "clinvar_vcf_index is required "
        "for workflow_application samples"
        in result.stderr
    )


def test_workflow_application_requires_small_variant_vcf(
    tmp_path: Path,
) -> None:
    project = write_project_config(tmp_path)
    runs = write_runs(tmp_path)
    bam, bai = create_alignment_files(tmp_path)

    (
        small_vcf,
        small_index,
        clinvar_vcf,
        clinvar_index,
    ) = create_variant_files(tmp_path)

    row = workflow_application_row(
        bam,
        bai,
        small_vcf,
        small_index,
        clinvar_vcf,
        clinvar_index,
    )
    row[5] = "NA"

    samples = write_samples(tmp_path, [row])
    result = run_validator(project, samples, runs)

    assert result.returncode == 1
    assert (
        "small_variant_vcf is required "
        "for workflow_application samples"
        in result.stderr
    )


def test_workflow_application_requires_small_variant_vcf_index(
    tmp_path: Path,
) -> None:
    project = write_project_config(tmp_path)
    runs = write_runs(tmp_path)
    bam, bai = create_alignment_files(tmp_path)

    (
        small_vcf,
        small_index,
        clinvar_vcf,
        clinvar_index,
    ) = create_variant_files(tmp_path)

    row = workflow_application_row(
        bam,
        bai,
        small_vcf,
        small_index,
        clinvar_vcf,
        clinvar_index,
    )
    row[6] = "NA"

    samples = write_samples(tmp_path, [row])
    result = run_validator(project, samples, runs)

    assert result.returncode == 1
    assert (
        "small_variant_vcf_index is required "
        "for workflow_application samples"
        in result.stderr
    )
