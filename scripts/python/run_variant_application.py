#!/usr/bin/env python3

"""
Run the small-variant application workflow for selected samples.

For each included workflow_application sample, the runner creates:
1. A technical small-variant summary.
2. A complete ClinVar-annotated variant table.
3. A review-candidate table.

It then creates cohort-level TSV and Excel summaries.

The workflow reports EPI2ME and ClinVar annotations. It does not perform
independent clinical classification or diagnostic interpretation.
"""

from __future__ import annotations

import argparse
import csv
import subprocess
import sys
from pathlib import Path


TRUE_VALUES = {"true", "yes", "1"}
MISSING_VALUES = {"", ".", "NA", "N/A", "None", "none"}


def parse_arguments() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Run variant summaries for workflow_application samples."
    )
    parser.add_argument(
        "--samples",
        required=True,
        type=Path,
        help="Selected samples.tsv file.",
    )
    parser.add_argument(
        "--scripts-dir",
        required=True,
        type=Path,
        help="Directory containing the clinical Python scripts.",
    )
    parser.add_argument(
        "--output-root",
        required=True,
        type=Path,
        help="Root directory for workflow-application outputs.",
    )
    return parser.parse_args()


def is_included(value: str) -> bool:
    return value.strip().lower() in TRUE_VALUES


def is_missing(value: str | None) -> bool:
    return value is None or value.strip() in MISSING_VALUES


def read_samples(path: Path) -> list[dict[str, str]]:
    if not path.is_file():
        raise ValueError(f"Samples metadata not found: {path}")

    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")

        required_columns = {
            "sample_id",
            "analysis_role",
            "include",
            "small_variant_vcf",
            "small_variant_vcf_index",
            "clinvar_vcf",
            "clinvar_vcf_index",
        }

        if reader.fieldnames is None:
            raise ValueError("Samples metadata has no header.")

        missing_columns = required_columns - set(reader.fieldnames)

        if missing_columns:
            raise ValueError(
                "Samples metadata is missing columns: "
                + ", ".join(sorted(missing_columns))
            )

        selected = [
            row
            for row in reader
            if (
                row["analysis_role"].strip()
                == "workflow_application"
                and is_included(row["include"])
            )
        ]

    if not selected:
        raise ValueError(
            "No included workflow_application samples were found."
        )

    return selected


def require_file(
    sample_id: str,
    label: str,
    value: str | None,
) -> Path:
    if is_missing(value):
        raise ValueError(
            f"Sample '{sample_id}': {label} is missing."
        )

    path = Path(value).expanduser()

    if not path.is_file():
        raise ValueError(
            f"Sample '{sample_id}': {label} not found: {path}"
        )

    return path


def run_command(command: list[str], description: str) -> None:
    print(f"[INFO] {description}")

    try:
        subprocess.run(command, check=True)
    except FileNotFoundError as exc:
        raise RuntimeError(
            f"Required executable was not found: {command[0]}"
        ) from exc
    except subprocess.CalledProcessError as exc:
        raise RuntimeError(
            f"Command failed with exit code {exc.returncode}: "
            + " ".join(command)
        ) from exc


def main() -> int:
    args = parse_arguments()

    samples_path = args.samples.expanduser().resolve()
    scripts_dir = args.scripts_dir.expanduser().resolve()
    output_root = args.output_root.expanduser().resolve()

    technical_script = scripts_dir / "summarize_clinical_variants.py"
    extraction_script = (
        scripts_dir / "extract_epi2me_clinical_variants.py"
    )
    cohort_script = scripts_dir / "summarize_clinical_cohort.py"

    for script_path in (
        technical_script,
        extraction_script,
        cohort_script,
    ):
        if not script_path.is_file():
            print(
                f"[ERROR] Required script not found: {script_path}",
                file=sys.stderr,
            )
            return 1

    try:
        samples = read_samples(samples_path)

        print(
            f"[INFO] Workflow-application samples: {len(samples)}"
        )

        for row in samples:
            sample_id = row["sample_id"].strip()

            small_variant_vcf = require_file(
                sample_id,
                "small_variant_vcf",
                row["small_variant_vcf"],
            )
            require_file(
                sample_id,
                "small_variant_vcf_index",
                row["small_variant_vcf_index"],
            )

            clinvar_vcf = require_file(
                sample_id,
                "clinvar_vcf",
                row["clinvar_vcf"],
            )
            require_file(
                sample_id,
                "clinvar_vcf_index",
                row["clinvar_vcf_index"],
            )

            sample_root = (
                output_root / sample_id / "small_variants"
            )
            summary_dir = sample_root / "summary"
            variants_dir = sample_root / "variants"

            summary_dir.mkdir(parents=True, exist_ok=True)
            variants_dir.mkdir(parents=True, exist_ok=True)

            technical_output = (
                summary_dir / f"{sample_id}.variant_summary.tsv"
            )
            annotated_output = (
                variants_dir
                / f"{sample_id}.all_annotated_variants.tsv"
            )
            review_output = (
                variants_dir
                / f"{sample_id}.review_candidates.tsv"
            )

            print()
            print(f"[INFO] Processing sample: {sample_id}")

            run_command(
                [
                    sys.executable,
                    str(technical_script),
                    "--sample-id",
                    sample_id,
                    "--input-vcf",
                    str(small_variant_vcf),
                    "--output-tsv",
                    str(technical_output),
                ],
                f"Creating technical variant summary: {sample_id}",
            )

            run_command(
                [
                    sys.executable,
                    str(extraction_script),
                    "--sample-id",
                    sample_id,
                    "--input-vcf",
                    str(clinvar_vcf),
                    "--output-all",
                    str(annotated_output),
                    "--output-review",
                    str(review_output),
                ],
                f"Creating ClinVar variant tables: {sample_id}",
            )

        cohort_dir = output_root / "cohort_summary"
        cohort_tsv = (
            cohort_dir / "workflow_application_variant_summary.tsv"
        )
        cohort_xlsx = (
            cohort_dir / "workflow_application_variant_summary.xlsx"
        )

        run_command(
            [
                sys.executable,
                str(cohort_script),
                "--samples",
                str(samples_path),
                "--outputs-root",
                str(output_root),
                "--output-tsv",
                str(cohort_tsv),
                "--output-xlsx",
                str(cohort_xlsx),
            ],
            "Creating workflow-application cohort summaries",
        )

    except (ValueError, RuntimeError) as exc:
        print(f"[ERROR] {exc}", file=sys.stderr)
        return 1

    print()
    print(
        f"[INFO] Variant application completed for "
        f"{len(samples)} sample(s)."
    )
    print(f"[INFO] Cohort TSV: {cohort_tsv}")
    print(f"[INFO] Cohort Excel: {cohort_xlsx}")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
