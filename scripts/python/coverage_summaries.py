#!/usr/bin/env python3

"""
Create gene-level coverage summaries from region-level coverage results.

The script joins region-level coverage metrics to gene annotations using
BED 0-based, half-open coordinates. BAM files are not read again.

A target interval may be assigned to more than one gene. In that case, the
same interval-level coverage metrics are reported independently for each gene.
"""

from __future__ import annotations

import argparse
import csv
import math
import statistics
import sys
from collections import defaultdict
from pathlib import Path
from typing import Any


def parse_arguments() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Create gene-level coverage summaries."
    )
    parser.add_argument(
        "--region-coverage",
        required=True,
        type=Path,
        help="Region-level coverage TSV",
    )
    parser.add_argument(
        "--annotation",
        required=True,
        type=Path,
        help="Gene-to-target annotation TSV",
    )
    parser.add_argument(
        "--per-sample-output",
        required=True,
        type=Path,
        help="Gene-level coverage output per sample",
    )
    parser.add_argument(
        "--cohort-output",
        required=True,
        type=Path,
        help="Gene-level summary across samples",
    )
    return parser.parse_args()


def read_tsv(path: Path) -> list[dict[str, str]]:
    if not path.is_file():
        raise ValueError(f"Input file not found: {path}")

    with path.open(
        "r",
        encoding="utf-8",
        newline="",
    ) as handle:
        reader = csv.DictReader(handle, delimiter="\t")

        if reader.fieldnames is None:
            raise ValueError(f"TSV file has no header: {path}")

        rows = list(reader)

    if not rows:
        raise ValueError(f"TSV file contains no data rows: {path}")

    return rows


def make_target_key(
    chromosome: str,
    start: str | int,
    end: str | int,
) -> str:
    return f"{chromosome}:{int(start)}-{int(end)}"


def parse_float(value: str, label: str) -> float:
    try:
        number = float(value)
    except ValueError as exc:
        raise ValueError(
            f"Invalid numeric value for {label}: '{value}'"
        ) from exc

    if not math.isfinite(number):
        raise ValueError(
            f"Non-finite numeric value for {label}: '{value}'"
        )

    return number


def parse_int(value: str, label: str) -> int:
    try:
        return int(float(value))
    except ValueError as exc:
        raise ValueError(
            f"Invalid integer value for {label}: '{value}'"
        ) from exc


def build_annotation_index(
    annotation_rows: list[dict[str, str]],
) -> dict[str, list[dict[str, str]]]:
    required = {
        "gene_symbol",
        "chromosome",
        "bed_start_0based",
        "bed_end_0based_exclusive",
        "target_key",
    }

    missing = required - set(annotation_rows[0])

    if missing:
        raise ValueError(
            "Annotation table is missing columns: "
            + ", ".join(sorted(missing))
        )

    index: dict[str, list[dict[str, str]]] = defaultdict(list)

    for row in annotation_rows:
        calculated_key = make_target_key(
            row["chromosome"],
            row["bed_start_0based"],
            row["bed_end_0based_exclusive"],
        )

        if calculated_key != row["target_key"]:
            raise ValueError(
                "Annotation target_key does not match coordinates: "
                f"{row['target_key']} versus {calculated_key}"
            )

        index[row["target_key"]].append(row)

    return dict(index)


def validate_region_columns(
    region_rows: list[dict[str, str]],
) -> list[int]:
    required = {
        "sample_id",
        "chromosome",
        "bed_start_0based",
        "bed_end_0based_exclusive",
        "region_length_bp",
        "target_bases",
        "mean_depth",
        "median_depth",
        "minimum_depth",
        "maximum_depth",
        "depth_standard_deviation",
        "depth_coefficient_of_variation",
        "zero_depth_bases",
        "zero_depth_percentage",
    }

    missing = required - set(region_rows[0])

    if missing:
        raise ValueError(
            "Region coverage table is missing columns: "
            + ", ".join(sorted(missing))
        )

    threshold_columns = [
        column
        for column in region_rows[0]
        if column.startswith("percentage_ge_")
        and column.endswith("x")
    ]

    if not threshold_columns:
        raise ValueError(
            "No percentage_ge_<threshold>x columns were found."
        )

    thresholds = []

    for column in threshold_columns:
        threshold_text = (
            column
            .removeprefix("percentage_ge_")
            .removesuffix("x")
        )

        thresholds.append(int(threshold_text))

    return sorted(thresholds)


def create_gene_rows(
    region_rows: list[dict[str, str]],
    annotation_index: dict[str, list[dict[str, str]]],
    thresholds: list[int],
) -> list[dict[str, Any]]:
    output_rows: list[dict[str, Any]] = []
    matched_targets: set[str] = set()

    for region in region_rows:
        target_key = make_target_key(
            region["chromosome"],
            region["bed_start_0based"],
            region["bed_end_0based_exclusive"],
        )

        annotations = annotation_index.get(target_key)

        if not annotations:
            raise ValueError(
                "No gene annotation was found for coverage target: "
                f"{target_key}"
            )

        matched_targets.add(target_key)

        region_length = parse_int(
            region["region_length_bp"],
            "region_length_bp",
        )
        target_bases = parse_int(
            region["target_bases"],
            "target_bases",
        )

        if region_length != target_bases:
            raise ValueError(
                f"Target length mismatch for {target_key}: "
                f"region_length_bp={region_length}, "
                f"target_bases={target_bases}"
            )

        for annotation in annotations:
            result: dict[str, Any] = {
                "sample_id": region["sample_id"],
                "gene_symbol": annotation["gene_symbol"],
                "target_key": target_key,
                "chromosome": region["chromosome"],
                "bed_start_0based": int(
                    region["bed_start_0based"]
                ),
                "bed_end_0based_exclusive": int(
                    region["bed_end_0based_exclusive"]
                ),
                "genomic_start_1based": int(
                    region["bed_start_0based"]
                ) + 1,
                "genomic_end_1based_inclusive": int(
                    region["bed_end_0based_exclusive"]
                ),
                "target_length_bp": target_bases,
                "mean_depth": parse_float(
                    region["mean_depth"],
                    "mean_depth",
                ),
                "median_depth": parse_float(
                    region["median_depth"],
                    "median_depth",
                ),
                "minimum_depth": parse_int(
                    region["minimum_depth"],
                    "minimum_depth",
                ),
                "maximum_depth": parse_int(
                    region["maximum_depth"],
                    "maximum_depth",
                ),
                "depth_standard_deviation": parse_float(
                    region["depth_standard_deviation"],
                    "depth_standard_deviation",
                ),
                "depth_coefficient_of_variation": (
                    region[
                        "depth_coefficient_of_variation"
                    ]
                ),
                "zero_depth_bases": parse_int(
                    region["zero_depth_bases"],
                    "zero_depth_bases",
                ),
                "zero_depth_percentage": parse_float(
                    region["zero_depth_percentage"],
                    "zero_depth_percentage",
                ),
                "annotation_match_status": annotation.get(
                    "match_status",
                    "",
                ),
                "annotation_gene_coverage_percentage": (
                    annotation.get(
                        "coverage_percentage",
                        "",
                    )
                ),
            }

            for threshold in thresholds:
                result[f"bases_ge_{threshold}x"] = parse_int(
                    region[f"bases_ge_{threshold}x"],
                    f"bases_ge_{threshold}x",
                )
                result[f"percentage_ge_{threshold}x"] = (
                    parse_float(
                        region[f"percentage_ge_{threshold}x"],
                        f"percentage_ge_{threshold}x",
                    )
                )

            output_rows.append(result)

    annotation_targets = set(annotation_index)
    missing_coverage_targets = annotation_targets - matched_targets

    if missing_coverage_targets:
        raise ValueError(
            "Annotation targets absent from region coverage: "
            + ", ".join(sorted(missing_coverage_targets))
        )

    return output_rows


def create_cohort_summary(
    gene_rows: list[dict[str, Any]],
    thresholds: list[int],
) -> list[dict[str, Any]]:
    rows_by_gene: dict[str, list[dict[str, Any]]] = defaultdict(list)

    for row in gene_rows:
        rows_by_gene[str(row["gene_symbol"])].append(row)

    summaries: list[dict[str, Any]] = []

    for gene_symbol in sorted(rows_by_gene):
        rows = rows_by_gene[gene_symbol]

        target_keys = {
            str(row["target_key"])
            for row in rows
        }

        sample_ids = {
            str(row["sample_id"])
            for row in rows
        }

        result: dict[str, Any] = {
            "gene_symbol": gene_symbol,
            "number_of_samples": len(sample_ids),
            "number_of_target_intervals": len(target_keys),
            "target_keys": ";".join(sorted(target_keys)),
            "mean_of_sample_mean_depth": statistics.fmean(
                float(row["mean_depth"])
                for row in rows
            ),
            "minimum_sample_mean_depth": min(
                float(row["mean_depth"])
                for row in rows
            ),
            "maximum_sample_mean_depth": max(
                float(row["mean_depth"])
                for row in rows
            ),
            "mean_zero_depth_percentage": statistics.fmean(
                float(row["zero_depth_percentage"])
                for row in rows
            ),
        }

        for threshold in thresholds:
            values = [
                float(row[f"percentage_ge_{threshold}x"])
                for row in rows
            ]

            result[
                f"mean_percentage_ge_{threshold}x"
            ] = statistics.fmean(values)

            result[
                f"minimum_percentage_ge_{threshold}x"
            ] = min(values)

            result[
                f"maximum_percentage_ge_{threshold}x"
            ] = max(values)

        summaries.append(result)

    return summaries


def format_value(value: Any) -> str:
    if isinstance(value, float):
        return f"{value:.6f}"

    return str(value)


def write_tsv(
    rows: list[dict[str, Any]],
    output_path: Path,
) -> None:
    if not rows:
        raise ValueError(
            f"No output rows were generated for {output_path}"
        )

    output_path.parent.mkdir(parents=True, exist_ok=True)

    with output_path.open(
        "w",
        encoding="utf-8",
        newline="",
    ) as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=list(rows[0]),
            delimiter="\t",
        )
        writer.writeheader()

        for row in rows:
            writer.writerow(
                {
                    key: format_value(value)
                    for key, value in row.items()
                }
            )


def main() -> int:
    args = parse_arguments()

    try:
        region_rows = read_tsv(args.region_coverage)
        annotation_rows = read_tsv(args.annotation)

        thresholds = validate_region_columns(region_rows)

        annotation_index = build_annotation_index(
            annotation_rows
        )

        gene_rows = create_gene_rows(
            region_rows=region_rows,
            annotation_index=annotation_index,
            thresholds=thresholds,
        )

        cohort_rows = create_cohort_summary(
            gene_rows=gene_rows,
            thresholds=thresholds,
        )

        write_tsv(
            gene_rows,
            args.per_sample_output,
        )

        write_tsv(
            cohort_rows,
            args.cohort_output,
        )

    except ValueError as exc:
        print(f"[ERROR] {exc}", file=sys.stderr)
        return 1

    genes = {
        str(row["gene_symbol"])
        for row in gene_rows
    }

    samples = {
        str(row["sample_id"])
        for row in gene_rows
    }

    print(f"[INFO] Samples: {len(samples)}")
    print(f"[INFO] Genes: {len(genes)}")
    print(
        f"[INFO] Gene-sample records: {len(gene_rows)}"
    )
    print(
        f"[INFO] Per-sample output: "
        f"{args.per_sample_output}"
    )
    print(
        f"[INFO] Cohort output: "
        f"{args.cohort_output}"
    )

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
