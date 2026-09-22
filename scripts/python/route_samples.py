#!/usr/bin/env python3

"""
Determine the analytical modules required by the included samples.

This script does not run any bioinformatics analysis. It reads samples.tsv
and reports the execution plan based on sample roles and available inputs.
"""

from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path

import yaml
from typing import Any


TRUE_VALUES = {"true", "yes", "1"}
EMPTY_VALUES = {"", "NA", "N/A", "none", "None", "."}


def parse_arguments() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Determine role-aware pipeline execution plan."
    )
    parser.add_argument(
        "--samples",
        required=True,
        type=Path,
        help="Path to samples.tsv",
    )
    parser.add_argument(
        "--project-config",
        required=True,
        type=Path,
        help="Path to the project YAML configuration.",
    )
    parser.add_argument(
        "--format",
        choices={"text", "json", "shell"},
        default="text",
        help="Output format",
    )
    return parser.parse_args()


def is_empty(value: str) -> bool:
    return value.strip() in EMPTY_VALUES


def is_included(value: str) -> bool:
    return value.strip().lower() in TRUE_VALUES


def read_included_samples(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        return [
            row
            for row in reader
            if is_included(row.get("include", "false"))
        ]


def build_execution_plan(
    samples: list[dict[str, str]],
    project_config: dict[str, Any],
) -> dict[str, Any]:
    role_counts = {
        "method_development": 0,
        "validation_control": 0,
        "workflow_application": 0,
    }

    sample_ids_by_role: dict[str, list[str]] = {
        role: []
        for role in role_counts
    }

    workflow_application_with_vcf: list[str] = []
    workflow_application_without_vcf: list[str] = []
    method_development_run_ids: set[str] = set()

    for sample in samples:
        sample_id = sample["sample_id"].strip()
        role = sample["analysis_role"].strip()

        if role not in role_counts:
            continue

        role_counts[role] += 1
        sample_ids_by_role[role].append(sample_id)

        if role == "method_development":
            run_id = sample.get("run_id", "").strip()
            if run_id:
                method_development_run_ids.add(run_id)

        if role == "workflow_application":
            if is_empty(sample["small_variant_vcf"]):
                workflow_application_without_vcf.append(sample_id)
            else:
                workflow_application_with_vcf.append(sample_id)

    modules = {
        "input_validation": True,
        "alignment_qc": bool(samples),
        "target_metrics": bool(samples),
        "coverage_analysis": bool(samples),
        "coverage_summaries": bool(samples),
        "method_development_summary": (
            role_counts["method_development"] > 0
        ),
        "method_development_reproducibility": (
            len(method_development_run_ids) == 2
        ),
        "giab_validation": (
            role_counts["validation_control"] > 0
        ),
        "variant_application": bool(
            workflow_application_with_vcf
        ),
        "off_target_analysis": bool(
            project_config.get("modules", {}).get(
                "off_target_analysis",
                False,
            )
        ),
        "igv_batch_preparation": bool(
            project_config.get("modules", {}).get(
                "igv_batch_preparation",
                False,
            )
        ),
        "sv_phasing_summary": bool(
            project_config.get("modules", {}).get(
                "sv_phasing_summary",
                False,
            )
            and role_counts["workflow_application"] > 0
        ),
        "igv_review_batch_preparation": bool(
            project_config.get("modules", {}).get(
                "igv_review_batch_preparation",
                False,
            )
        ) and role_counts["workflow_application"] > 0,
    }

    return {
        "included_sample_count": len(samples),
        "role_counts": role_counts,
        "sample_ids_by_role": sample_ids_by_role,
        "workflow_application_with_vcf": (
            workflow_application_with_vcf
        ),
        "workflow_application_without_vcf": (
            workflow_application_without_vcf
        ),
        "modules": modules,
    }


def print_text_plan(plan: dict[str, Any]) -> None:
    print("Role-aware execution plan")
    print("=========================")
    print()
    print(f"Included samples: {plan['included_sample_count']}")
    print()

    print("Sample roles:")
    for role, count in plan["role_counts"].items():
        print(f"  {role}: {count}")

    print()
    print("Selected modules:")

    for module_name, enabled in plan["modules"].items():
        status = "RUN" if enabled else "SKIP"
        print(f"  [{status}] {module_name}")

    samples_without_vcf = (
        plan["workflow_application_without_vcf"]
    )

    if samples_without_vcf:
        print()
        print(
            "Workflow-application samples without a small-variant VCF:"
        )
        for sample_id in samples_without_vcf:
            print(f"  - {sample_id}")
        print(
            "  Common BAM-based QC will run, but variant "
            "summarisation will be skipped for these samples."
        )


def print_shell_plan(plan: dict[str, Any]) -> None:
    """
    Print module decisions as Bash-compatible variable assignments.

    Only fixed internal module names are converted to shell variables.
    Sample identifiers and file paths are not emitted in shell format.
    """

    for module_name, enabled in plan["modules"].items():
        variable_name = (
            "RUN_"
            + module_name.upper()
        )

        value = "true" if enabled else "false"
        print(f"{variable_name}={value}")


def main() -> int:
    args = parse_arguments()

    if not args.project_config.is_file():
        raise SystemExit(
            f"ERROR: Project configuration not found: {args.project_config}"
        )

    with args.project_config.open("r", encoding="utf-8") as handle:
        project_config = yaml.safe_load(handle) or {}

    samples = read_included_samples(args.samples)
    plan = build_execution_plan(samples, project_config)

    if args.format == "json":
        print(json.dumps(plan, indent=2))
    elif args.format == "shell":
        print_shell_plan(plan)
    else:
        print_text_plan(plan)

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
