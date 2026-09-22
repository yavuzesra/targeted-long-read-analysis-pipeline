#!/usr/bin/env python3

"""
Calculate alignment-record and target-base enrichment metrics.

Alignment-record metrics
------------------------
Eligible record:
    An alignment record passing the configured panel-analysis SAM flag filter.

On-target record:
    An eligible alignment record overlapping at least one merged target
    interval by one or more reference bases.

Base/depth metrics
------------------
Target depth sum:
    Sum of alignment depth across the merged target territory.

Whole-reference depth sum:
    Sum of alignment depth across the complete reference territory.

Target-base enrichment factor:
    Mean target depth divided by mean genome-wide depth.

Sequencing-run summaries
------------------------
Samples are grouped using the run_id field in samples.tsv. Run-level values
are unweighted arithmetic means of the corresponding sample-level metrics.

The original BAM files are read directly. No permanent filtered BAM files are
created.
"""

from __future__ import annotations

import argparse
import csv
import statistics
import subprocess
import sys
import tempfile
from pathlib import Path
from typing import Any, TextIO

import yaml


TRUE_VALUES = {"true", "yes", "1"}

PER_SAMPLE_COLUMNS = (
    "sample_id",
    "analysis_role",
    "run_id",
    "bam_path",
    "panel_analysis_exclude_flags",
    "eligible_mapped_primary_records",
    "on_target_records",
    "off_target_records",
    "on_target_record_percentage",
    "off_target_record_percentage",
    "target_length_bp",
    "reference_length_bp",
    "target_depth_sum",
    "total_depth_sum",
    "on_target_aligned_base_fraction_percent",
    "mean_target_depth",
    "mean_genome_wide_depth",
    "target_base_enrichment_factor",
)

RUN_COLUMNS = (
    "run_id",
    "n_samples",
    "mean_on_target_record_percentage",
    "mean_on_target_aligned_base_fraction_percent",
    "mean_target_base_enrichment_factor",
)


def parse_arguments() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Calculate per-sample alignment-record and target-base "
            "enrichment metrics and sequencing-run arithmetic means."
        )
    )
    parser.add_argument(
        "--project-config",
        required=True,
        type=Path,
        help="Path to project.yaml.",
    )
    parser.add_argument(
        "--samples",
        required=True,
        type=Path,
        help="Path to samples.tsv.",
    )
    parser.add_argument(
        "--filtering-config",
        required=True,
        type=Path,
        help="Path to filtering.yaml.",
    )
    parser.add_argument(
        "--output",
        required=True,
        type=Path,
        help="Per-sample output TSV path.",
    )
    parser.add_argument(
        "--run-output",
        required=True,
        type=Path,
        help="Sequencing-run summary output TSV path.",
    )
    return parser.parse_args()


def run_command(command: list[str]) -> str:
    result = subprocess.run(
        command,
        capture_output=True,
        text=True,
        check=False,
    )

    if result.returncode != 0:
        message = result.stderr.strip() or result.stdout.strip()
        raise RuntimeError(
            f"Command failed ({result.returncode}): "
            f"{' '.join(command)}\n{message}"
        )

    return result.stdout


def load_yaml(path: Path) -> dict[str, Any]:
    if not path.is_file():
        raise ValueError(f"YAML file not found: {path}")

    with path.open("r", encoding="utf-8") as handle:
        data = yaml.safe_load(handle)

    if not isinstance(data, dict):
        raise ValueError(f"YAML root must be a mapping: {path}")

    return data


def read_resources(project_config: Path) -> tuple[Path, Path]:
    config = load_yaml(project_config)

    try:
        resources = config["resources"]
        reference_fai = Path(
            str(resources["reference_fai"])
        ).expanduser()
        target_bed = Path(
            str(resources["target_bed"])
        ).expanduser()
    except (KeyError, TypeError) as exc:
        raise ValueError(
            "project.yaml must define resources.reference_fai "
            "and resources.target_bed."
        ) from exc

    if not reference_fai.is_file():
        raise ValueError(
            f"Reference FASTA index not found: {reference_fai}"
        )

    if not target_bed.is_file():
        raise ValueError(f"Target BED not found: {target_bed}")

    return reference_fai, target_bed


def read_panel_filter(filtering_config: Path) -> int:
    config = load_yaml(filtering_config)

    try:
        return int(
            config["filters"]["panel_analysis"]["exclude_flags"]
        )
    except (KeyError, TypeError, ValueError) as exc:
        raise ValueError(
            "filtering.yaml does not contain a valid "
            "filters.panel_analysis.exclude_flags value."
        ) from exc


def read_included_samples(
    samples_path: Path,
) -> list[dict[str, str]]:
    if not samples_path.is_file():
        raise ValueError(
            f"Sample metadata not found: {samples_path}"
        )

    with samples_path.open(
        "r",
        encoding="utf-8",
        newline="",
    ) as handle:
        reader = csv.DictReader(handle, delimiter="\t")

        required = {
            "sample_id",
            "analysis_role",
            "run_id",
            "bam_path",
            "include",
        }

        if reader.fieldnames is None:
            raise ValueError("samples.tsv has no header.")

        missing = required - set(reader.fieldnames)

        if missing:
            raise ValueError(
                "samples.tsv is missing required columns: "
                + ", ".join(sorted(missing))
            )

        samples = [
            row
            for row in reader
            if row["include"].strip().lower() in TRUE_VALUES
        ]

    if not samples:
        raise ValueError("No samples with include=true were found.")

    seen_sample_ids: set[str] = set()

    for sample in samples:
        sample_id = sample["sample_id"].strip()

        if not sample_id:
            raise ValueError(
                "An included sample has an empty sample_id."
            )

        if sample_id in seen_sample_ids:
            raise ValueError(
                f"Duplicate included sample_id: {sample_id}"
            )

        seen_sample_ids.add(sample_id)

        for column in ("analysis_role", "run_id", "bam_path"):
            if not sample[column].strip():
                raise ValueError(
                    f"Missing {column} for {sample_id}."
                )

    return samples


def read_reference_lengths(
    reference_fai: Path,
) -> dict[str, int]:
    lengths: dict[str, int] = {}

    with reference_fai.open("r", encoding="utf-8") as handle:
        for line_number, line in enumerate(handle, start=1):
            fields = line.rstrip("\n").split("\t")

            if len(fields) < 2:
                raise ValueError(
                    f"Invalid FAI line {line_number}: "
                    f"{line.rstrip()}"
                )

            chromosome = fields[0]

            try:
                length = int(fields[1])
            except ValueError as exc:
                raise ValueError(
                    f"Invalid chromosome length in FAI line "
                    f"{line_number}: {fields[1]}"
                ) from exc

            if length <= 0:
                raise ValueError(
                    f"Non-positive chromosome length in FAI line "
                    f"{line_number}: {length}"
                )

            lengths[chromosome] = length

    if not lengths:
        raise ValueError("Reference FAI contains no sequences.")

    return lengths


def read_bed_intervals(
    target_bed: Path,
) -> list[tuple[str, int, int]]:
    intervals: list[tuple[str, int, int]] = []

    with target_bed.open("r", encoding="utf-8") as handle:
        for line_number, line in enumerate(handle, start=1):
            stripped = line.strip()

            if not stripped or stripped.startswith(
                ("#", "track", "browser")
            ):
                continue

            fields = stripped.split("\t")

            if len(fields) < 3:
                raise ValueError(
                    f"BED line {line_number} has fewer than "
                    "three columns."
                )

            chromosome = fields[0]

            try:
                start = int(fields[1])
                end = int(fields[2])
            except ValueError as exc:
                raise ValueError(
                    f"BED line {line_number} contains "
                    "non-integer coordinates."
                ) from exc

            if start < 0:
                raise ValueError(
                    f"BED line {line_number} has a negative start."
                )

            if end <= start:
                raise ValueError(
                    f"BED line {line_number} must satisfy end > start."
                )

            intervals.append((chromosome, start, end))

    if not intervals:
        raise ValueError("Target BED contains no valid intervals.")

    return intervals


def merge_intervals(
    intervals: list[tuple[str, int, int]],
) -> list[tuple[str, int, int]]:
    sorted_intervals = sorted(
        intervals,
        key=lambda item: (item[0], item[1], item[2]),
    )

    merged: list[tuple[str, int, int]] = []

    for chromosome, start, end in sorted_intervals:
        if not merged:
            merged.append((chromosome, start, end))
            continue

        last_chromosome, last_start, last_end = merged[-1]

        if chromosome == last_chromosome and start <= last_end:
            merged[-1] = (
                last_chromosome,
                last_start,
                max(last_end, end),
            )
        else:
            merged.append((chromosome, start, end))

    return merged


def validate_bed_against_reference(
    intervals: list[tuple[str, int, int]],
    reference_lengths: dict[str, int],
) -> None:
    errors: list[str] = []

    for chromosome, start, end in intervals:
        if chromosome not in reference_lengths:
            errors.append(
                f"BED chromosome '{chromosome}' is absent "
                "from the reference FAI."
            )
            continue

        chromosome_length = reference_lengths[chromosome]

        if end > chromosome_length:
            errors.append(
                f"BED interval {chromosome}:{start}-{end} "
                f"exceeds chromosome length {chromosome_length}."
            )

    if errors:
        raise ValueError(" ".join(errors))


def write_merged_bed(
    intervals: list[tuple[str, int, int]],
    handle: TextIO,
) -> None:
    for chromosome, start, end in intervals:
        handle.write(f"{chromosome}\t{start}\t{end}\n")

    handle.flush()


def get_bam_chromosomes(bam_path: Path) -> set[str]:
    header = run_command(
        ["samtools", "view", "-H", str(bam_path)]
    )

    chromosomes: set[str] = set()

    for line in header.splitlines():
        if not line.startswith("@SQ"):
            continue

        for field in line.split("\t"):
            if field.startswith("SN:"):
                chromosomes.add(field.removeprefix("SN:"))
                break

    if not chromosomes:
        raise ValueError(
            f"No @SQ reference sequences were found in: {bam_path}"
        )

    return chromosomes


def validate_bed_against_bam(
    bed_chromosomes: set[str],
    bam_path: Path,
) -> None:
    bam_chromosomes = get_bam_chromosomes(bam_path)
    missing = sorted(bed_chromosomes - bam_chromosomes)

    if missing:
        raise ValueError(
            f"BAM and BED chromosome naming are incompatible for "
            f"{bam_path}. Missing BAM chromosomes: "
            + ", ".join(missing)
        )


def samtools_count(
    bam_path: Path,
    exclude_flags: int,
    target_bed: Path | None = None,
) -> int:
    command = [
        "samtools",
        "view",
        "-c",
        "-F",
        str(exclude_flags),
    ]

    if target_bed is not None:
        command.extend(["-L", str(target_bed)])

    command.append(str(bam_path))

    output = run_command(command).strip()

    try:
        return int(output)
    except ValueError as exc:
        raise RuntimeError(
            f"Invalid samtools count for {bam_path}: {output}"
        ) from exc


def samtools_depth_sum(
    bam_path: Path,
    exclude_flags: int,
    target_bed: Path | None = None,
) -> int:
    samtools_command = [
        "samtools",
        "depth",
        "-G",
        str(exclude_flags),
        "-q",
        "0",
        "-Q",
        "0",
    ]

    if target_bed is not None:
        samtools_command.extend(["-b", str(target_bed)])

    samtools_command.append(str(bam_path))

    awk_command = [
        "awk",
        "{sum += $3} END {printf \"%.0f\", sum + 0}",
    ]

    samtools_process = subprocess.Popen(
        samtools_command,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True,
    )

    if samtools_process.stdout is None:
        raise RuntimeError(
            "Unable to open samtools depth output stream."
        )

    awk_process = subprocess.Popen(
        awk_command,
        stdin=samtools_process.stdout,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True,
    )

    samtools_process.stdout.close()

    awk_stdout, awk_stderr = awk_process.communicate()

    samtools_stderr = (
        samtools_process.stderr.read()
        if samtools_process.stderr is not None
        else ""
    )

    samtools_returncode = samtools_process.wait()

    if samtools_returncode != 0:
        raise RuntimeError(
            f"samtools depth failed for {bam_path}: "
            f"{samtools_stderr.strip()}"
        )

    if awk_process.returncode != 0:
        raise RuntimeError(
            f"awk depth summation failed for {bam_path}: "
            f"{awk_stderr.strip()}"
        )

    output = awk_stdout.strip()

    try:
        return int(output)
    except ValueError as exc:
        raise RuntimeError(
            f"Invalid depth sum for {bam_path}: {output}"
        ) from exc


def calculate_record_values(
    eligible_records: int,
    on_target_records: int,
) -> dict[str, int | float]:
    if eligible_records <= 0:
        raise ValueError(
            "Eligible alignment-record count must be greater than zero."
        )

    if on_target_records < 0:
        raise ValueError(
            "On-target alignment-record count must not be negative."
        )

    if on_target_records > eligible_records:
        raise ValueError(
            "On-target record count exceeds eligible record count."
        )

    off_target_records = eligible_records - on_target_records

    on_target_percentage = (
        100.0 * on_target_records / eligible_records
    )

    return {
        "eligible_mapped_primary_records": eligible_records,
        "on_target_records": on_target_records,
        "off_target_records": off_target_records,
        "on_target_record_percentage": on_target_percentage,
        "off_target_record_percentage":
            100.0 - on_target_percentage,
    }


def calculate_base_values(
    target_depth_sum: int,
    total_depth_sum: int,
    target_length_bp: int,
    reference_length_bp: int,
) -> dict[str, float]:
    if target_length_bp <= 0:
        raise ValueError(
            "Target territory length must be greater than zero."
        )

    if reference_length_bp <= 0:
        raise ValueError(
            "Reference territory length must be greater than zero."
        )

    if total_depth_sum <= 0:
        raise ValueError(
            "Whole-reference depth sum must be greater than zero."
        )

    if target_depth_sum < 0:
        raise ValueError(
            "Target depth sum must not be negative."
        )

    if target_depth_sum > total_depth_sum:
        raise ValueError(
            "Target depth sum exceeds whole-reference depth sum."
        )

    aligned_base_fraction_percent = (
        100.0 * target_depth_sum / total_depth_sum
    )

    mean_target_depth = (
        target_depth_sum / target_length_bp
    )

    mean_genome_wide_depth = (
        total_depth_sum / reference_length_bp
    )

    target_base_enrichment_factor = (
        mean_target_depth / mean_genome_wide_depth
    )

    return {
        "on_target_aligned_base_fraction_percent":
            aligned_base_fraction_percent,
        "mean_target_depth": mean_target_depth,
        "mean_genome_wide_depth": mean_genome_wide_depth,
        "target_base_enrichment_factor":
            target_base_enrichment_factor,
    }


def calculate_sample_metrics(
    sample: dict[str, str],
    exclude_flags: int,
    merged_target_bed: Path,
    target_length_bp: int,
    reference_length_bp: int,
    bed_chromosomes: set[str],
) -> dict[str, int | float | str]:
    sample_id = sample["sample_id"].strip()
    analysis_role = sample["analysis_role"].strip()
    run_id = sample["run_id"].strip()
    bam_path = Path(sample["bam_path"]).expanduser()

    if not bam_path.is_file():
        raise ValueError(
            f"Sample '{sample_id}' BAM not found: {bam_path}"
        )

    quickcheck = subprocess.run(
        ["samtools", "quickcheck", "-v", str(bam_path)],
        capture_output=True,
        text=True,
        check=False,
    )

    if quickcheck.returncode != 0:
        message = (
            quickcheck.stderr.strip()
            or quickcheck.stdout.strip()
            or "samtools quickcheck failed."
        )
        raise ValueError(
            f"Sample '{sample_id}' failed BAM integrity check: "
            f"{message}"
        )

    validate_bed_against_bam(
        bed_chromosomes=bed_chromosomes,
        bam_path=bam_path,
    )

    eligible_records = samtools_count(
        bam_path=bam_path,
        exclude_flags=exclude_flags,
    )

    on_target_records = samtools_count(
        bam_path=bam_path,
        exclude_flags=exclude_flags,
        target_bed=merged_target_bed,
    )

    total_depth_sum = samtools_depth_sum(
        bam_path=bam_path,
        exclude_flags=exclude_flags,
    )

    target_depth_sum = samtools_depth_sum(
        bam_path=bam_path,
        exclude_flags=exclude_flags,
        target_bed=merged_target_bed,
    )

    record_values = calculate_record_values(
        eligible_records=eligible_records,
        on_target_records=on_target_records,
    )

    base_values = calculate_base_values(
        target_depth_sum=target_depth_sum,
        total_depth_sum=total_depth_sum,
        target_length_bp=target_length_bp,
        reference_length_bp=reference_length_bp,
    )

    return {
        "sample_id": sample_id,
        "analysis_role": analysis_role,
        "run_id": run_id,
        "bam_path": str(bam_path),
        "panel_analysis_exclude_flags": exclude_flags,
        **record_values,
        "target_length_bp": target_length_bp,
        "reference_length_bp": reference_length_bp,
        "target_depth_sum": target_depth_sum,
        "total_depth_sum": total_depth_sum,
        **base_values,
    }


def summarise_runs(
    sample_results: list[dict[str, int | float | str]],
) -> list[dict[str, int | float | str]]:
    grouped: dict[
        str,
        list[dict[str, int | float | str]],
    ] = {}

    for row in sample_results:
        run_id = str(row["run_id"])
        grouped.setdefault(run_id, []).append(row)

    summaries: list[dict[str, int | float | str]] = []

    for run_id, rows in sorted(grouped.items()):
        summaries.append(
            {
                "run_id": run_id,
                "n_samples": len(rows),
                "mean_on_target_record_percentage":
                    statistics.fmean(
                        float(
                            row[
                                "on_target_record_percentage"
                            ]
                        )
                        for row in rows
                    ),
                "mean_on_target_aligned_base_fraction_percent":
                    statistics.fmean(
                        float(
                            row[
                                "on_target_aligned_base_fraction_percent"
                            ]
                        )
                        for row in rows
                    ),
                "mean_target_base_enrichment_factor":
                    statistics.fmean(
                        float(
                            row[
                                "target_base_enrichment_factor"
                            ]
                        )
                        for row in rows
                    ),
            }
        )

    return summaries


def format_value(value: int | float | str) -> str:
    if isinstance(value, float):
        return f"{value:.6f}"

    return str(value)


def write_tsv(
    output_path: Path,
    rows: list[dict[str, int | float | str]],
    fieldnames: tuple[str, ...],
) -> None:
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
            extrasaction="ignore",
        )

        writer.writeheader()

        for row in rows:
            writer.writerow(
                {
                    field: format_value(row[field])
                    for field in fieldnames
                }
            )


def main() -> int:
    args = parse_arguments()

    try:
        reference_fai, target_bed = read_resources(
            args.project_config
        )

        exclude_flags = read_panel_filter(
            args.filtering_config
        )

        samples = read_included_samples(args.samples)

        reference_lengths = read_reference_lengths(
            reference_fai
        )

        intervals = read_bed_intervals(target_bed)

        validate_bed_against_reference(
            intervals=intervals,
            reference_lengths=reference_lengths,
        )

        merged_intervals = merge_intervals(intervals)

        target_length_bp = sum(
            end - start
            for _, start, end in merged_intervals
        )

        reference_length_bp = sum(
            reference_lengths.values()
        )

        bed_chromosomes = {
            chromosome
            for chromosome, _, _ in merged_intervals
        }

        print(
            f"[INFO] Included samples: {len(samples)}"
        )
        print(
            f"[INFO] Original target intervals: {len(intervals)}"
        )
        print(
            f"[INFO] Merged target intervals: "
            f"{len(merged_intervals)}"
        )
        print(
            f"[INFO] Merged target territory: "
            f"{target_length_bp} bp"
        )
        print(
            f"[INFO] Reference territory: "
            f"{reference_length_bp} bp"
        )
        print(
            f"[INFO] Excluded SAM flags: {exclude_flags}"
        )

        with tempfile.NamedTemporaryFile(
            mode="w",
            encoding="utf-8",
            suffix=".merged_targets.bed",
            delete=False,
        ) as temporary_bed:
            write_merged_bed(
                intervals=merged_intervals,
                handle=temporary_bed,
            )
            merged_target_bed = Path(temporary_bed.name)

        try:
            sample_results = []

            for sample in samples:
                sample_id = sample["sample_id"].strip()
                run_id = sample["run_id"].strip()

                print(
                    "[INFO] Calculating target metrics: "
                    f"{sample_id} ({run_id})",
                    flush=True,
                )

                sample_results.append(
                    calculate_sample_metrics(
                        sample=sample,
                        exclude_flags=exclude_flags,
                        merged_target_bed=merged_target_bed,
                        target_length_bp=target_length_bp,
                        reference_length_bp=reference_length_bp,
                        bed_chromosomes=bed_chromosomes,
                    )
                )
        finally:
            merged_target_bed.unlink(missing_ok=True)

        run_results = summarise_runs(sample_results)

        write_tsv(
            output_path=args.output,
            rows=sample_results,
            fieldnames=PER_SAMPLE_COLUMNS,
        )

        write_tsv(
            output_path=args.run_output,
            rows=run_results,
            fieldnames=RUN_COLUMNS,
        )

    except (ValueError, RuntimeError, OSError) as exc:
        print(f"[ERROR] {exc}", file=sys.stderr)
        return 1

    print(
        f"[INFO] Per-sample target metrics completed for "
        f"{len(sample_results)} sample(s)."
    )
    print(f"[INFO] Per-sample output: {args.output}")
    print(f"[INFO] Run-level output: {args.run_output}")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
