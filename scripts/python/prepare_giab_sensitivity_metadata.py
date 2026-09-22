#!/usr/bin/env python3

"""Build temporary validation metadata for optional GIAB sensitivity runs."""

from __future__ import annotations

import argparse
import csv
import sys
from pathlib import Path


SAMPLE_FIELDS = [
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

SENSITIVITY_FIELDS = {
    "benchmark_id",
    "source_sample_id",
    "truth_vcf",
    "truth_vcf_index",
    "truth_bed",
    "truth_version",
    "include",
    "notes",
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Prepare metadata for optional GIAB truth-set sensitivity analyses."
    )
    parser.add_argument("--samples", required=True, type=Path)
    parser.add_argument("--sensitivity-config", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    return parser.parse_args()


def read_table(path: Path) -> tuple[list[str], list[dict[str, str]]]:
    if not path.is_file():
        raise ValueError(f"Input table not found: {path}")
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames is None:
            raise ValueError(f"Input table has no header: {path}")
        return reader.fieldnames, list(reader)


def require_columns(
    observed: list[str], required: set[str], label: str
) -> None:
    missing = required - set(observed)
    if missing:
        raise ValueError(
            f"{label} is missing columns: " + ", ".join(sorted(missing))
        )


def is_true(value: str) -> bool:
    return value.strip().lower() == "true"


def require_file(value: str, field: str, benchmark_id: str) -> str:
    if not value or value.upper() == "NA":
        raise ValueError(f"{benchmark_id}: {field} is missing")
    path = Path(value)
    if not path.is_file():
        raise ValueError(f"{benchmark_id}: {field} not found: {path}")
    return str(path)


def build_rows(
    samples: list[dict[str, str]],
    sensitivity_rows: list[dict[str, str]],
) -> list[dict[str, str]]:
    source_by_id = {row["sample_id"]: row for row in samples}
    output: list[dict[str, str]] = []
    seen_benchmark_ids: set[str] = set()

    for row in sensitivity_rows:
        if not is_true(row.get("include", "")):
            continue

        benchmark_id = row["benchmark_id"].strip()
        source_id = row["source_sample_id"].strip()
        if not benchmark_id:
            raise ValueError("Included sensitivity row has an empty benchmark_id")
        if benchmark_id in seen_benchmark_ids:
            raise ValueError(f"Duplicate benchmark_id: {benchmark_id}")
        seen_benchmark_ids.add(benchmark_id)

        if source_id not in source_by_id:
            raise ValueError(f"Source sample not found: {source_id}")
        source = source_by_id[source_id]
        if source.get("analysis_role") != "validation_control":
            raise ValueError(f"{source_id} is not a validation_control sample")
        if not is_true(source.get("include", "")):
            raise ValueError(f"{source_id} is not included in samples metadata")

        truth_version = row["truth_version"].strip()
        if not truth_version or truth_version.upper() == "NA":
            raise ValueError(f"{benchmark_id}: truth_version is missing")

        output.append(
            {
                "sample_id": benchmark_id,
                "analysis_role": "validation_control",
                "run_id": source.get("run_id", "NA"),
                "bam_path": source.get("bam_path", "NA"),
                "bai_path": source.get("bai_path", "NA"),
                "small_variant_vcf": require_file(
                    source.get("small_variant_vcf", ""),
                    "small_variant_vcf",
                    benchmark_id,
                ),
                "small_variant_vcf_index": require_file(
                    source.get("small_variant_vcf_index", ""),
                    "small_variant_vcf_index",
                    benchmark_id,
                ),
                "clinvar_vcf": "NA",
                "clinvar_vcf_index": "NA",
                "truth_vcf": require_file(
                    row.get("truth_vcf", ""), "truth_vcf", benchmark_id
                ),
                "truth_vcf_index": require_file(
                    row.get("truth_vcf_index", ""),
                    "truth_vcf_index",
                    benchmark_id,
                ),
                "truth_bed": require_file(
                    row.get("truth_bed", ""), "truth_bed", benchmark_id
                ),
                "truth_version": truth_version,
                "sex": source.get("sex", "unknown"),
                "include": "true",
                "notes": (
                    f"Sensitivity benchmark for {source_id}; "
                    + row.get("notes", "").strip()
                ).rstrip("; "),
            }
        )

    if not output:
        raise ValueError(
            "No sensitivity analyses are enabled; set include=true after "
            "verifying the alternative truth-set paths"
        )
    return output


def write_rows(path: Path, rows: list[dict[str, str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=SAMPLE_FIELDS,
            delimiter="\t",
            lineterminator="\n",
        )
        writer.writeheader()
        writer.writerows(rows)


def main() -> int:
    args = parse_args()
    try:
        sample_fields, samples = read_table(args.samples)
        sensitivity_fields, sensitivity = read_table(args.sensitivity_config)
        require_columns(sample_fields, set(SAMPLE_FIELDS), "Samples table")
        require_columns(
            sensitivity_fields, SENSITIVITY_FIELDS, "Sensitivity config"
        )
        rows = build_rows(samples, sensitivity)
        write_rows(args.output, rows)
    except ValueError as exc:
        print(f"[ERROR] {exc}", file=sys.stderr)
        return 1

    print(f"[INFO] Enabled sensitivity benchmarks: {len(rows)}")
    print(f"[INFO] Temporary metadata: {args.output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
