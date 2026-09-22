#!/usr/bin/env python3

"""
Validate project configuration, run metadata, and sample metadata.

The validator does not modify any input file.

Supported analysis roles:
    - method_development
    - validation_control
    - workflow_application

Exit codes:
    0: validation completed successfully
    1: one or more validation errors were detected
"""

from __future__ import annotations

import argparse
import csv
import sys
from pathlib import Path
from typing import Any

import yaml


ALLOWED_ROLES = {
    "method_development",
    "validation_control",
    "workflow_application",
}

ALLOWED_SEX_VALUES = {
    "male",
    "female",
    "unknown",
    "NA",
}

TRUE_VALUES = {"true", "yes", "1"}
FALSE_VALUES = {"false", "no", "0"}

SAMPLE_REQUIRED_COLUMNS = [
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

RUN_REQUIRED_COLUMNS = [
    "run_id",
    "dataset_label",
    "workflow",
    "workflow_version",
    "bam_root",
    "reference_build",
    "notes",
]

EMPTY_VALUES = {"", "NA", "N/A", "none", "None", "."}


def parse_arguments() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Validate QIAseq xHYB pipeline metadata."
    )
    parser.add_argument(
        "--project-config",
        required=True,
        type=Path,
        help="Path to project.yaml",
    )
    parser.add_argument(
        "--samples",
        required=True,
        type=Path,
        help="Path to samples.tsv",
    )
    parser.add_argument(
        "--runs",
        required=True,
        type=Path,
        help="Path to runs.tsv",
    )
    return parser.parse_args()


def is_empty(value: Any) -> bool:
    return str(value).strip() in EMPTY_VALUES


def normalise_path(value: str) -> Path | None:
    if is_empty(value):
        return None
    return Path(value).expanduser()


def read_yaml(path: Path) -> dict[str, Any]:
    if not path.is_file():
        raise ValueError(f"Configuration file not found: {path}")

    with path.open("r", encoding="utf-8") as handle:
        data = yaml.safe_load(handle)

    if not isinstance(data, dict):
        raise ValueError(f"YAML root must be a mapping: {path}")

    return data


def read_tsv(
    path: Path,
    required_columns: list[str],
) -> list[dict[str, str]]:
    if not path.is_file():
        raise ValueError(f"TSV file not found: {path}")

    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")

        if reader.fieldnames is None:
            raise ValueError(f"TSV file has no header: {path}")

        missing = [
            column
            for column in required_columns
            if column not in reader.fieldnames
        ]

        if missing:
            raise ValueError(
                f"Missing required column(s) in {path}: "
                + ", ".join(missing)
            )

        rows = list(reader)

    if not rows:
        raise ValueError(f"TSV file contains no data rows: {path}")

    return rows


def parse_include(value: str) -> bool:
    normalised = value.strip().lower()

    if normalised in TRUE_VALUES:
        return True

    if normalised in FALSE_VALUES:
        return False

    raise ValueError(
        f"Invalid include value '{value}'. "
        "Use true or false."
    )


def validate_project_resources(
    project_config: dict[str, Any],
    errors: list[str],
) -> None:
    resources = project_config.get("resources")

    if not isinstance(resources, dict):
        errors.append(
            "project.yaml must contain a 'resources' mapping."
        )
        return

    required_resources = [
        "reference_fasta",
        "reference_fai",
        "target_bed",
        "annotation_xlsx",
    ]

    for key in required_resources:
        value = resources.get(key)

        if value is None or is_empty(value):
            errors.append(
                f"Missing resource path in project.yaml: resources.{key}"
            )
            continue

        path = Path(str(value)).expanduser()

        if not path.is_file():
            errors.append(
                f"Resource file not found for resources.{key}: {path}"
            )


def validate_runs(
    run_rows: list[dict[str, str]],
    errors: list[str],
) -> set[str]:
    run_ids: set[str] = set()

    for line_number, row in enumerate(run_rows, start=2):
        run_id = row["run_id"].strip()

        if not run_id:
            errors.append(
                f"runs.tsv line {line_number}: run_id is empty."
            )
            continue

        if run_id in run_ids:
            errors.append(
                f"runs.tsv line {line_number}: "
                f"duplicate run_id '{run_id}'."
            )

        run_ids.add(run_id)

        bam_root = normalise_path(row["bam_root"])

        if bam_root is None:
            errors.append(
                f"runs.tsv line {line_number}: "
                f"bam_root is empty for run '{run_id}'."
            )
        elif not bam_root.is_dir():
            errors.append(
                f"runs.tsv line {line_number}: "
                f"bam_root directory not found for run '{run_id}': "
                f"{bam_root}"
            )

    return run_ids


def validate_optional_file_pair(
    sample_id: str,
    first_label: str,
    first_value: str,
    second_label: str,
    second_value: str,
    errors: list[str],
) -> None:
    first_empty = is_empty(first_value)
    second_empty = is_empty(second_value)

    if first_empty and second_empty:
        return

    if first_empty != second_empty:
        errors.append(
            f"Sample '{sample_id}': {first_label} and "
            f"{second_label} must either both be provided or both be NA."
        )
        return

    first_path = normalise_path(first_value)
    second_path = normalise_path(second_value)

    if first_path is not None and not first_path.is_file():
        errors.append(
            f"Sample '{sample_id}': {first_label} not found: "
            f"{first_path}"
        )

    if second_path is not None and not second_path.is_file():
        errors.append(
            f"Sample '{sample_id}': {second_label} not found: "
            f"{second_path}"
        )


def validate_samples(
    sample_rows: list[dict[str, str]],
    valid_run_ids: set[str],
    errors: list[str],
    warnings: list[str],
) -> dict[str, int]:
    sample_ids: set[str] = set()

    role_counts = {
        "method_development": 0,
        "validation_control": 0,
        "workflow_application": 0,
    }

    included_count = 0

    for line_number, row in enumerate(sample_rows, start=2):
        sample_id = row["sample_id"].strip()
        role = row["analysis_role"].strip()
        run_id = row["run_id"].strip()

        if not sample_id:
            errors.append(
                f"samples.tsv line {line_number}: sample_id is empty."
            )
            continue

        if sample_id in sample_ids:
            errors.append(
                f"samples.tsv line {line_number}: "
                f"duplicate sample_id '{sample_id}'."
            )

        sample_ids.add(sample_id)

        try:
            include_sample = parse_include(row["include"])
        except ValueError as exc:
            errors.append(
                f"Sample '{sample_id}': {exc}"
            )
            continue

        if not include_sample:
            continue

        included_count += 1

        if role not in ALLOWED_ROLES:
            errors.append(
                f"Sample '{sample_id}': invalid analysis_role "
                f"'{role}'. Allowed values: "
                + ", ".join(sorted(ALLOWED_ROLES))
            )
            continue

        role_counts[role] += 1

        if run_id not in valid_run_ids:
            errors.append(
                f"Sample '{sample_id}': run_id '{run_id}' "
                "is not defined in runs.tsv."
            )

        sex = row["sex"].strip()

        if sex not in ALLOWED_SEX_VALUES:
            errors.append(
                f"Sample '{sample_id}': invalid sex value '{sex}'. "
                "Use male, female, unknown, or NA."
            )

        bam_path = normalise_path(row["bam_path"])
        bai_path = normalise_path(row["bai_path"])

        if bam_path is None:
            errors.append(
                f"Sample '{sample_id}': bam_path is required."
            )
        elif not bam_path.is_file():
            errors.append(
                f"Sample '{sample_id}': BAM file not found: {bam_path}"
            )

        if bai_path is None:
            errors.append(
                f"Sample '{sample_id}': bai_path is required."
            )
        elif not bai_path.is_file():
            errors.append(
                f"Sample '{sample_id}': BAM index not found: {bai_path}"
            )

        validate_optional_file_pair(
            sample_id=sample_id,
            first_label="small_variant_vcf",
            first_value=row["small_variant_vcf"],
            second_label="small_variant_vcf_index",
            second_value=row["small_variant_vcf_index"],
            errors=errors,
        )

        validate_optional_file_pair(
            sample_id=sample_id,
            first_label="clinvar_vcf",
            first_value=row["clinvar_vcf"],
            second_label="clinvar_vcf_index",
            second_value=row["clinvar_vcf_index"],
            errors=errors,
        )

        truth_fields = [
            row["truth_vcf"],
            row["truth_vcf_index"],
            row["truth_bed"],
            row["truth_version"],
        ]

        if role == "validation_control":
            required_truth_labels = [
                "truth_vcf",
                "truth_vcf_index",
                "truth_bed",
                "truth_version",
            ]

            for label in required_truth_labels:
                if is_empty(row[label]):
                    errors.append(
                        f"Sample '{sample_id}': {label} is required "
                        "for validation_control samples."
                    )

            if is_empty(row["small_variant_vcf"]):
                errors.append(
                    f"Sample '{sample_id}': small_variant_vcf is "
                    "required for validation_control samples."
                )

            for label in [
                "truth_vcf",
                "truth_vcf_index",
                "truth_bed",
            ]:
                path = normalise_path(row[label])

                if path is not None and not path.is_file():
                    errors.append(
                        f"Sample '{sample_id}': {label} not found: "
                        f"{path}"
                    )

        else:
            if any(not is_empty(value) for value in truth_fields):
                errors.append(
                    f"Sample '{sample_id}': truth-set fields are only "
                    "allowed for validation_control samples."
                )

        if role == "workflow_application":
            required_variant_fields = [
                "small_variant_vcf",
                "small_variant_vcf_index",
                "clinvar_vcf",
                "clinvar_vcf_index",
            ]

            for label in required_variant_fields:
                if is_empty(row[label]):
                    errors.append(
                        f"Sample '{sample_id}': {label} is required "
                        "for workflow_application samples."
                    )

    if included_count == 0:
        errors.append(
            "No samples are enabled. At least one sample must have "
            "include=true."
        )

    role_counts["included_total"] = included_count
    return role_counts


def main() -> int:
    args = parse_arguments()

    errors: list[str] = []
    warnings: list[str] = []

    try:
        project_config = read_yaml(args.project_config)
        sample_rows = read_tsv(
            args.samples,
            SAMPLE_REQUIRED_COLUMNS,
        )
        run_rows = read_tsv(
            args.runs,
            RUN_REQUIRED_COLUMNS,
        )
    except ValueError as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 1

    validate_project_resources(project_config, errors)
    valid_run_ids = validate_runs(run_rows, errors)

    role_counts = validate_samples(
        sample_rows=sample_rows,
        valid_run_ids=valid_run_ids,
        errors=errors,
        warnings=warnings,
    )

    if warnings:
        print("Warnings:")
        for warning in warnings:
            print(f"  - {warning}")
        print()

    if errors:
        print("Metadata validation failed:", file=sys.stderr)
        for error in errors:
            print(f"  - {error}", file=sys.stderr)

        print(
            f"\nTotal errors: {len(errors)}",
            file=sys.stderr,
        )
        return 1

    print("Metadata validation completed successfully.")
    print()
    print("Included samples:")
    print(
        f"  method_development: "
        f"{role_counts['method_development']}"
    )
    print(
        f"  validation_control: "
        f"{role_counts['validation_control']}"
    )
    print(
        f"  workflow_application: "
        f"{role_counts['workflow_application']}"
    )
    print(
        f"  total: {role_counts['included_total']}"
    )

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
