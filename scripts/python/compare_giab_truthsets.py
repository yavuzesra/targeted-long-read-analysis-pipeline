#!/usr/bin/env python3

"""Compare primary and sensitivity GIAB results by variant type."""

from __future__ import annotations

import argparse
import csv
import sys
from pathlib import Path


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--sensitivity-config", required=True, type=Path)
    parser.add_argument("--primary-summary", required=True, type=Path)
    parser.add_argument("--primary-manifest", required=True, type=Path)
    parser.add_argument("--sensitivity-summary", required=True, type=Path)
    parser.add_argument("--sensitivity-manifest", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    return parser.parse_args()


def read_tsv(path: Path) -> list[dict[str, str]]:
    if not path.is_file():
        raise ValueError(f"Required table not found: {path}")
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames is None:
            raise ValueError(f"Table has no header: {path}")
        return list(reader)


def included_configs(path: Path) -> list[dict[str, str]]:
    rows = read_tsv(path)
    selected = [
        row for row in rows if row.get("include", "").strip().lower() == "true"
    ]
    if not selected:
        raise ValueError("Sensitivity config has no enabled rows")
    return selected


def metrics_by_sample(
    rows: list[dict[str, str]], sample_id: str
) -> dict[str, dict[str, str]]:
    selected = {
        row["variant_type"].upper(): row
        for row in rows
        if row.get("sample_id") == sample_id
    }
    missing = {"SNP", "INDEL"} - set(selected)
    if missing:
        raise ValueError(
            f"Missing benchmark rows for {sample_id}: "
            + ", ".join(sorted(missing))
        )
    return selected


def manifest_row(
    rows: list[dict[str, str]], sample_id: str
) -> dict[str, str]:
    selected = [row for row in rows if row.get("sample_id") == sample_id]
    if len(selected) != 1:
        raise ValueError(
            f"Expected one manifest row for {sample_id}, found {len(selected)}"
        )
    return selected[0]


def output_rows(
    analysis_type: str,
    benchmark_id: str,
    source_sample_id: str,
    metrics: dict[str, dict[str, str]],
    manifest: dict[str, str],
) -> list[dict[str, str]]:
    rows: list[dict[str, str]] = []
    for variant_type in ("SNP", "INDEL"):
        item = metrics[variant_type]
        rows.append(
            {
                "analysis_type": analysis_type,
                "benchmark_id": benchmark_id,
                "source_sample_id": source_sample_id,
                "truth_version": item["truth_version"],
                "benchmarkable_bases": manifest["benchmarkable_bases"],
                "variant_type": variant_type,
                "truth_total": item["truth_total"],
                "truth_tp": item["truth_tp"],
                "truth_fn": item["truth_fn"],
                "query_total": item["query_total"],
                "query_tp": item["query_tp"],
                "query_fp": item["query_fp"],
                "query_unk": item["query_unk"],
                "precision": item["precision"],
                "recall": item["recall"],
                "f1_score": item["f1_score"],
            }
        )
    return rows


def build_comparison(args: argparse.Namespace) -> list[dict[str, str]]:
    configs = included_configs(args.sensitivity_config)
    primary_summary = read_tsv(args.primary_summary)
    primary_manifest = read_tsv(args.primary_manifest)
    sensitivity_summary = read_tsv(args.sensitivity_summary)
    sensitivity_manifest = read_tsv(args.sensitivity_manifest)
    output: list[dict[str, str]] = []

    for config in configs:
        benchmark_id = config["benchmark_id"]
        source_id = config["source_sample_id"]
        output.extend(
            output_rows(
                "primary",
                source_id,
                source_id,
                metrics_by_sample(primary_summary, source_id),
                manifest_row(primary_manifest, source_id),
            )
        )
        output.extend(
            output_rows(
                "sensitivity",
                benchmark_id,
                source_id,
                metrics_by_sample(sensitivity_summary, benchmark_id),
                manifest_row(sensitivity_manifest, benchmark_id),
            )
        )
    return output


def write_tsv(path: Path, rows: list[dict[str, str]]) -> None:
    fieldnames = [
        "analysis_type",
        "benchmark_id",
        "source_sample_id",
        "truth_version",
        "benchmarkable_bases",
        "variant_type",
        "truth_total",
        "truth_tp",
        "truth_fn",
        "query_total",
        "query_tp",
        "query_fp",
        "query_unk",
        "precision",
        "recall",
        "f1_score",
    ]
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=fieldnames, delimiter="\t", lineterminator="\n"
        )
        writer.writeheader()
        writer.writerows(rows)


def main() -> int:
    args = parse_args()
    try:
        rows = build_comparison(args)
        write_tsv(args.output, rows)
    except (KeyError, ValueError) as exc:
        print(f"[ERROR] {exc}", file=sys.stderr)
        return 1
    print(f"[INFO] Truth-set comparison rows: {len(rows)}")
    print(f"[INFO] Output: {args.output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
