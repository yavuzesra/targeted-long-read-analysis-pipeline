#!/usr/bin/env python3

"""Create method-development summaries from completed shared-QC tables.

The module joins sample metadata with alignment, target-enrichment, and
panel-coverage outputs. It does not read BAM files or repeat upstream analyses.
"""

from __future__ import annotations

import argparse
import csv
import re
import statistics
from pathlib import Path


TRUE_VALUES = {"true", "yes", "1"}
GROUP_PATTERN = re.compile(
    r"^(?P<technical_group>.+)_(?P<technical_replicate>[0-9]+)$"
)

METRICS = (
    "primary_mapped_reads",
    "mean_primary_mapped_read_length_bp",
    "n50_primary_mapped_read_length_bp",
    "mean_alignment_accuracy_pct",
    "on_target_record_percentage",
    "on_target_aligned_base_fraction_percent",
    "target_base_enrichment_factor",
    "mean_depth",
    "percentage_ge_10x",
    "percentage_ge_20x",
    "zero_depth_percentage",
)

LIBRARY_COLUMNS = (
    "sample_id",
    "run_id",
    "technical_group",
    "technical_replicate",
    *METRICS,
    "panel_qc_pass",
)

GROUP_COLUMNS = (
    "summary_level",
    "run_id",
    "technical_group",
    "library_count",
    "panel_qc_pass_count",
    "panel_qc_pass_percentage",
    *(
        column
        for metric in METRICS
        for column in (f"mean_{metric}", f"sd_{metric}")
    ),
)


def parse_arguments() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Summarise completed shared-QC results for "
            "method-development technical libraries."
        )
    )
    parser.add_argument("--samples", required=True, type=Path)
    parser.add_argument("--alignment-qc", required=True, type=Path)
    parser.add_argument("--target-metrics", required=True, type=Path)
    parser.add_argument("--panel-coverage", required=True, type=Path)
    parser.add_argument("--library-output", required=True, type=Path)
    parser.add_argument("--group-output", required=True, type=Path)
    return parser.parse_args()


def read_tsv(path: Path) -> list[dict[str, str]]:
    if not path.is_file():
        raise ValueError(f"Input file not found: {path}")

    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames is None:
            raise ValueError(f"TSV file has no header: {path}")
        return list(reader)


def index_rows(
    rows: list[dict[str, str]],
    source_label: str,
) -> dict[str, dict[str, str]]:
    indexed: dict[str, dict[str, str]] = {}
    for row in rows:
        sample_id = row.get("sample_id", "").strip()
        if not sample_id:
            raise ValueError(
                f"{source_label} contains a row without sample_id."
            )
        if sample_id in indexed:
            raise ValueError(
                f"{source_label} contains duplicate sample_id: {sample_id}"
            )
        indexed[sample_id] = row
    return indexed


def method_samples(
    rows: list[dict[str, str]],
) -> list[dict[str, str]]:
    selected = [
        row
        for row in rows
        if row.get("analysis_role", "").strip()
        == "method_development"
        and row.get("include", "").strip().lower() in TRUE_VALUES
    ]
    if not selected:
        raise ValueError(
            "No included method_development samples were found."
        )
    return selected


def parse_group(sample_id: str) -> tuple[str, int]:
    match = GROUP_PATTERN.fullmatch(sample_id)
    if match is None:
        raise ValueError(
            "Method-development sample_id must end with an underscore "
            f"and a numeric technical replicate: {sample_id}"
        )
    return (
        match.group("technical_group"),
        int(match.group("technical_replicate")),
    )


def require_value(
    row: dict[str, str],
    column: str,
    sample_id: str,
    source_label: str,
) -> str:
    value = row.get(column, "").strip()
    if value in {"", "NA", "N/A", "."}:
        raise ValueError(
            f"{source_label}: missing {column} for {sample_id}."
        )
    return value


def create_library_rows(
    samples: list[dict[str, str]],
    alignment_rows: list[dict[str, str]],
    target_rows: list[dict[str, str]],
    coverage_rows: list[dict[str, str]],
) -> list[dict[str, object]]:
    alignment = index_rows(alignment_rows, "alignment QC")
    targets = index_rows(target_rows, "target metrics")
    coverage = index_rows(coverage_rows, "panel coverage")

    library_rows: list[dict[str, object]] = []
    for sample in method_samples(samples):
        sample_id = sample["sample_id"].strip()
        run_id = sample.get("run_id", "").strip()
        if not run_id:
            raise ValueError(f"Missing run_id for {sample_id}.")

        missing_sources = [
            label
            for label, indexed in (
                ("alignment QC", alignment),
                ("target metrics", targets),
                ("panel coverage", coverage),
            )
            if sample_id not in indexed
        ]
        if missing_sources:
            raise ValueError(
                f"Missing shared-QC row for {sample_id}: "
                + ", ".join(missing_sources)
            )

        technical_group, technical_replicate = parse_group(sample_id)
        joined = {
            **alignment[sample_id],
            **targets[sample_id],
            **coverage[sample_id],
        }
        row: dict[str, object] = {
            "sample_id": sample_id,
            "run_id": run_id,
            "technical_group": technical_group,
            "technical_replicate": technical_replicate,
        }
        for metric in METRICS:
            row[metric] = float(
                require_value(joined, metric, sample_id, "shared QC")
            )
        row["panel_qc_pass"] = require_value(
            joined,
            "panel_qc_pass",
            sample_id,
            "panel coverage",
        ).lower()
        library_rows.append(row)

    return sorted(
        library_rows,
        key=lambda row: (
            str(row["run_id"]),
            str(row["technical_group"]),
            int(row["technical_replicate"]),
        ),
    )


def format_number(value: float) -> str:
    return f"{value:.6f}"


def summarise_group(
    rows: list[dict[str, object]],
    summary_level: str,
    run_id: str,
    technical_group: str,
) -> dict[str, object]:
    result: dict[str, object] = {
        "summary_level": summary_level,
        "run_id": run_id,
        "technical_group": technical_group,
        "library_count": len(rows),
    }
    pass_count = sum(
        str(row["panel_qc_pass"]).lower() == "true"
        for row in rows
    )
    result["panel_qc_pass_count"] = pass_count
    result["panel_qc_pass_percentage"] = format_number(
        100.0 * pass_count / len(rows)
    )

    for metric in METRICS:
        values = [float(row[metric]) for row in rows]
        result[f"mean_{metric}"] = format_number(
            statistics.fmean(values)
        )
        result[f"sd_{metric}"] = format_number(
            statistics.stdev(values) if len(values) > 1 else 0.0
        )
    return result


def create_group_rows(
    library_rows: list[dict[str, object]],
) -> list[dict[str, object]]:
    by_run: dict[str, list[dict[str, object]]] = {}
    by_technical_group: dict[
        tuple[str, str], list[dict[str, object]]
    ] = {}

    for row in library_rows:
        run_id = str(row["run_id"])
        technical_group = str(row["technical_group"])
        by_run.setdefault(run_id, []).append(row)
        by_technical_group.setdefault(
            (run_id, technical_group), []
        ).append(row)

    summaries = [
        summarise_group(rows, "experiment", run_id, "ALL")
        for run_id, rows in sorted(by_run.items())
    ]
    summaries.extend(
        summarise_group(
            rows,
            "technical_group",
            run_id,
            technical_group,
        )
        for (run_id, technical_group), rows
        in sorted(by_technical_group.items())
    )
    return summaries


def write_tsv(
    path: Path,
    rows: list[dict[str, object]],
    fieldnames: tuple[str, ...],
) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=fieldnames,
            delimiter="\t",
            extrasaction="ignore",
        )
        writer.writeheader()
        writer.writerows(rows)


def main() -> int:
    args = parse_arguments()
    library_rows = create_library_rows(
        samples=read_tsv(args.samples),
        alignment_rows=read_tsv(args.alignment_qc),
        target_rows=read_tsv(args.target_metrics),
        coverage_rows=read_tsv(args.panel_coverage),
    )
    group_rows = create_group_rows(library_rows)

    write_tsv(args.library_output, library_rows, LIBRARY_COLUMNS)
    write_tsv(args.group_output, group_rows, GROUP_COLUMNS)

    print(
        "[INFO] Method-development libraries summarised: "
        f"{len(library_rows)}"
    )
    print(f"[INFO] Library output: {args.library_output}")
    print(f"[INFO] Group output: {args.group_output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
