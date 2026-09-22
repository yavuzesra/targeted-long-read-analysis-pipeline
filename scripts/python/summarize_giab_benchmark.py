#!/usr/bin/env python3

"""
Create a standardized GIAB small-variant benchmark summary.

The script reads giab_benchmark_runs.tsv and extracts the PASS rows for
SNPs and INDELs from each hap.py summary CSV.

No combined average is calculated because controls may use different
GIAB truth-set versions.
"""

from __future__ import annotations

import argparse
import csv
import sys
from pathlib import Path


VARIANT_TYPES = ("SNP", "INDEL")


def parse_arguments() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Summarize hap.py GIAB benchmark results."
    )
    parser.add_argument(
        "--run-table",
        required=True,
        type=Path,
        help="giab_benchmark_runs.tsv produced by the benchmark runner",
    )
    parser.add_argument(
        "--output",
        required=True,
        type=Path,
        help="Combined benchmark summary TSV",
    )
    return parser.parse_args()


def read_run_table(path: Path) -> list[dict[str, str]]:
    if not path.is_file():
        raise ValueError(f"Benchmark run table not found: {path}")

    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")

        required = {
            "sample_id",
            "truth_version",
            "sex",
            "summary_csv",
            "status",
        }

        if reader.fieldnames is None:
            raise ValueError("Benchmark run table has no header.")

        missing = required - set(reader.fieldnames)

        if missing:
            raise ValueError(
                "Benchmark run table is missing columns: "
                + ", ".join(sorted(missing))
            )

        return list(reader)


def read_happy_pass_rows(
    summary_path: Path,
) -> dict[str, dict[str, str]]:
    if not summary_path.is_file():
        raise ValueError(f"hap.py summary CSV not found: {summary_path}")

    selected: dict[str, dict[str, str]] = {}

    with summary_path.open(
        "r",
        encoding="utf-8",
        newline="",
    ) as handle:
        reader = csv.DictReader(handle)

        if reader.fieldnames is None:
            raise ValueError(
                f"hap.py summary has no header: {summary_path}"
            )

        for row in reader:
            variant_type = row.get("Type", "").strip().upper()
            filter_status = row.get("Filter", "").strip().upper()

            if (
                variant_type in VARIANT_TYPES
                and filter_status == "PASS"
            ):
                selected[variant_type] = row

    missing_types = set(VARIANT_TYPES) - set(selected)

    if missing_types:
        raise ValueError(
            f"Missing PASS rows in {summary_path}: "
            + ", ".join(sorted(missing_types))
        )

    return selected


def numeric_value(
    row: dict[str, str],
    column: str,
) -> str:
    value = row.get(column, "").strip()

    if value == "":
        raise ValueError(f"Missing hap.py value: {column}")

    return value


def query_true_positives(row: dict[str, str]) -> str:
    """Return QUERY.TP, deriving it for hap.py summaries that omit the column."""
    value = row.get("QUERY.TP", "").strip()
    if value:
        return value

    query_total = int(float(numeric_value(row, "QUERY.TOTAL")))
    query_fp = int(float(numeric_value(row, "QUERY.FP")))
    query_unk = int(float(numeric_value(row, "QUERY.UNK")))
    query_tp = query_total - query_fp - query_unk
    if query_tp < 0:
        raise ValueError("Derived QUERY.TP is negative")
    return str(query_tp)


def build_summary(
    run_rows: list[dict[str, str]],
) -> list[dict[str, str]]:
    output_rows: list[dict[str, str]] = []

    for run_row in run_rows:
        if run_row["status"].strip().lower() != "completed":
            continue

        summary_path = Path(run_row["summary_csv"])
        happy_rows = read_happy_pass_rows(summary_path)

        for variant_type in VARIANT_TYPES:
            row = happy_rows[variant_type]

            output_rows.append(
                {
                    "sample_id": run_row["sample_id"],
                    "truth_version": run_row["truth_version"],
                    "sex": run_row["sex"],
                    "variant_type": variant_type,
                    "filter": "PASS",
                    "truth_total": numeric_value(
                        row,
                        "TRUTH.TOTAL",
                    ),
                    "truth_tp": numeric_value(
                        row,
                        "TRUTH.TP",
                    ),
                    "truth_fn": numeric_value(
                        row,
                        "TRUTH.FN",
                    ),
                    "query_total": numeric_value(
                        row,
                        "QUERY.TOTAL",
                    ),
                    "query_tp": query_true_positives(row),
                    "query_fp": numeric_value(
                        row,
                        "QUERY.FP",
                    ),
                    "query_unk": numeric_value(
                        row,
                        "QUERY.UNK",
                    ),
                    "recall": numeric_value(
                        row,
                        "METRIC.Recall",
                    ),
                    "precision": numeric_value(
                        row,
                        "METRIC.Precision",
                    ),
                    "f1_score": numeric_value(
                        row,
                        "METRIC.F1_Score",
                    ),
                }
            )

    if not output_rows:
        raise ValueError(
            "No completed GIAB benchmark rows were available."
        )

    return output_rows


def write_summary(
    rows: list[dict[str, str]],
    output_path: Path,
) -> None:
    fieldnames = [
        "sample_id",
        "truth_version",
        "sex",
        "variant_type",
        "filter",
        "truth_total",
        "truth_tp",
        "truth_fn",
        "query_total",
        "query_tp",
        "query_fp",
        "query_unk",
        "recall",
        "precision",
        "f1_score",
    ]

    output_path.parent.mkdir(parents=True, exist_ok=True)

    with output_path.open(
        "w",
        encoding="utf-8",
        newline="",
    ) as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=fieldnames,
            delimiter="\t",
        )
        writer.writeheader()
        writer.writerows(rows)


def main() -> int:
    args = parse_arguments()

    try:
        run_rows = read_run_table(args.run_table)
        summary_rows = build_summary(run_rows)
        write_summary(summary_rows, args.output)
    except ValueError as exc:
        print(f"[ERROR] {exc}", file=sys.stderr)
        return 1

    print(
        f"[INFO] GIAB benchmark summary rows: "
        f"{len(summary_rows)}"
    )
    print(f"[INFO] Output: {args.output}")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
