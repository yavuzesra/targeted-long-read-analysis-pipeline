#!/usr/bin/env python3

"""
Calculate alignment quality-control metrics from original wf-alignment BAM files.

The script:
    - reads included samples from samples.tsv;
    - uses SAM flag masks defined in filtering.yaml;
    - does not create permanent filtered BAM files;
    - writes one TSV row per sample.

Required external tool:
    samtools
"""

from __future__ import annotations

import argparse
import csv
import math
import statistics
import subprocess
import sys
from pathlib import Path
from typing import Iterable

import yaml


TRUE_VALUES = {"true", "yes", "1"}
EMPTY_VALUES = {"", "NA", "N/A", "none", "None", "."}


def parse_arguments() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Calculate BAM-based alignment QC metrics."
    )
    parser.add_argument(
        "--samples",
        required=True,
        type=Path,
        help="Path to samples.tsv",
    )
    parser.add_argument(
        "--filtering-config",
        required=True,
        type=Path,
        help="Path to filtering.yaml",
    )
    parser.add_argument(
        "--output",
        required=True,
        type=Path,
        help="Output TSV path",
    )
    return parser.parse_args()


def is_included(value: str) -> bool:
    return value.strip().lower() in TRUE_VALUES


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


def samtools_count(
    bam_path: Path,
    *,
    include_flags: int | None = None,
    exclude_flags: int | None = None,
) -> int:
    command = ["samtools", "view", "-c"]

    if include_flags is not None:
        command.extend(["-f", str(include_flags)])

    if exclude_flags is not None:
        command.extend(["-F", str(exclude_flags)])

    command.append(str(bam_path))

    output = run_command(command).strip()
    return int(output)


def read_filter_masks(path: Path) -> tuple[int, int, int]:
    if not path.is_file():
        raise ValueError(f"Filtering configuration not found: {path}")

    with path.open("r", encoding="utf-8") as handle:
        config = yaml.safe_load(handle)

    try:
        filters = config["filters"]

        primary_total = int(
            filters["primary_total"]["exclude_flags"]
        )
        primary_mapped = int(
            filters["primary_mapped"]["exclude_flags"]
        )
        panel_analysis = int(
            filters["panel_analysis"]["exclude_flags"]
        )
    except (KeyError, TypeError, ValueError) as exc:
        raise ValueError(
            "filtering.yaml does not contain valid filter definitions."
        ) from exc

    return primary_total, primary_mapped, panel_analysis


def read_included_samples(path: Path) -> list[dict[str, str]]:
    if not path.is_file():
        raise ValueError(f"Sample metadata not found: {path}")

    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")

        required_columns = {
            "sample_id",
            "bam_path",
            "include",
        }

        if reader.fieldnames is None:
            raise ValueError("samples.tsv has no header.")

        missing = required_columns - set(reader.fieldnames)

        if missing:
            raise ValueError(
                "samples.tsv is missing required columns: "
                + ", ".join(sorted(missing))
            )

        samples = [
            row
            for row in reader
            if is_included(row["include"])
        ]

    if not samples:
        raise ValueError("No samples with include=true were found.")

    return samples


def validate_bam(bam_path: Path) -> None:
    if not bam_path.is_file():
        raise ValueError(f"BAM file not found: {bam_path}")

    result = subprocess.run(
        ["samtools", "quickcheck", "-v", str(bam_path)],
        capture_output=True,
        text=True,
        check=False,
    )

    if result.returncode != 0:
        message = result.stderr.strip() or result.stdout.strip()
        raise ValueError(
            f"BAM integrity check failed for {bam_path}: {message}"
        )


def get_sort_order(bam_path: Path) -> str:
    header = run_command(
        ["samtools", "view", "-H", str(bam_path)]
    )

    for line in header.splitlines():
        if not line.startswith("@HD"):
            continue

        for field in line.split("\t"):
            if field.startswith("SO:"):
                return field.removeprefix("SO:")

    return "unknown"


def iter_primary_read_lengths(
    bam_path: Path,
    exclude_flags: int,
) -> Iterable[int]:
    process = subprocess.Popen(
        [
            "samtools",
            "view",
            "-F",
            str(exclude_flags),
            str(bam_path),
        ],
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True,
    )

    if process.stdout is None:
        raise RuntimeError(
            f"Unable to read samtools output for {bam_path}"
        )

    for line in process.stdout:
        fields = line.rstrip("\n").split("\t")

        if len(fields) < 10:
            continue

        sequence = fields[9]

        if sequence != "*":
            yield len(sequence)

    stderr = ""

    if process.stderr is not None:
        stderr = process.stderr.read().strip()

    return_code = process.wait()

    if return_code != 0:
        raise RuntimeError(
            f"samtools view failed for {bam_path}: {stderr}"
        )


def calculate_n50(lengths: list[int]) -> int | None:
    if not lengths:
        return None

    total_bases = sum(lengths)
    cumulative = 0

    for length in sorted(lengths, reverse=True):
        cumulative += length

        if cumulative >= total_bases / 2:
            return length

    return None


def calculate_read_length_metrics(
    bam_path: Path,
    exclude_flags: int,
) -> tuple[float | None, float | None, int | None]:
    lengths = list(
        iter_primary_read_lengths(
            bam_path=bam_path,
            exclude_flags=exclude_flags,
        )
    )

    if not lengths:
        return None, None, None

    mean_length = statistics.fmean(lengths)
    median_length = statistics.median(lengths)
    n50_length = calculate_n50(lengths)

    return mean_length, float(median_length), n50_length


def parse_samtools_stats(
    bam_path: Path,
    exclude_flags: int,
) -> dict[str, float | None]:
    """
    Calculate samtools statistics from streamed filtered BAM records.

    The filtered BAM stream is not written to disk. The exclusion mask is
    supplied by filtering.yaml so that alignment accuracy uses the same
    record set as the other panel analyses.
    """
    view_process = subprocess.Popen(
        [
            "samtools",
            "view",
            "-u",
            "-F",
            str(exclude_flags),
            str(bam_path),
        ],
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
    )

    if view_process.stdout is None:
        raise RuntimeError(
            f"Unable to stream filtered BAM records for {bam_path}"
        )

    stats_process = subprocess.run(
        ["samtools", "stats", "-"],
        stdin=view_process.stdout,
        capture_output=True,
        text=True,
        check=False,
    )

    # Allow samtools view to receive a closed pipe if samtools stats fails.
    view_process.stdout.close()

    view_stderr_bytes = b""
    if view_process.stderr is not None:
        view_stderr_bytes = view_process.stderr.read()

    view_return_code = view_process.wait()
    view_stderr = view_stderr_bytes.decode(
        "utf-8",
        errors="replace",
    ).strip()

    if view_return_code != 0:
        raise RuntimeError(
            f"samtools view failed for {bam_path}: {view_stderr}"
        )

    if stats_process.returncode != 0:
        message = (
            stats_process.stderr.strip()
            or stats_process.stdout.strip()
        )
        raise RuntimeError(
            f"samtools stats failed for {bam_path}: {message}"
        )

    error_rate: float | None = None

    for line in stats_process.stdout.splitlines():
        if not line.startswith("SN"):
            continue

        fields = line.split("\t")

        if len(fields) < 3:
            continue

        metric_name = fields[1].rstrip(":")
        metric_value = fields[2]

        if metric_name == "error rate":
            try:
                error_rate = float(metric_value)
            except ValueError:
                error_rate = None

    alignment_accuracy: float | None = None

    if error_rate is not None and math.isfinite(error_rate):
        alignment_accuracy = 100.0 * (1.0 - error_rate)

    return {
        "error_rate": error_rate,
        "alignment_accuracy_pct": alignment_accuracy,
    }


def format_value(
    value: int | float | str | None,
    decimals: int = 4,
) -> str:
    if value is None:
        return "NA"

    if isinstance(value, float):
        return f"{value:.{decimals}f}"

    return str(value)


def calculate_sample_metrics(
    sample: dict[str, str],
    primary_total_mask: int,
    primary_mapped_mask: int,
    panel_analysis_mask: int,
) -> dict[str, int | float | str | None]:
    sample_id = sample["sample_id"].strip()
    bam_path = Path(sample["bam_path"]).expanduser()

    validate_bam(bam_path)

    sort_order = get_sort_order(bam_path)

    if sort_order != "coordinate":
        raise ValueError(
            f"Sample '{sample_id}' BAM is not coordinate-sorted. "
            f"Detected sort order: {sort_order}"
        )

    total_records = samtools_count(bam_path)

    primary_total = samtools_count(
        bam_path,
        exclude_flags=primary_total_mask,
    )

    primary_mapped = samtools_count(
        bam_path,
        exclude_flags=primary_mapped_mask,
    )

    primary_unmapped = samtools_count(
        bam_path,
        include_flags=4,
        exclude_flags=primary_total_mask,
    )

    secondary = samtools_count(
        bam_path,
        include_flags=256,
    )

    supplementary = samtools_count(
        bam_path,
        include_flags=2048,
    )

    duplicate = samtools_count(
        bam_path,
        include_flags=1024,
    )

    qc_failed = samtools_count(
        bam_path,
        include_flags=512,
    )

    panel_analysis_records = samtools_count(
        bam_path,
        exclude_flags=panel_analysis_mask,
    )

    count_difference = (
        primary_total
        - primary_mapped
        - primary_unmapped
    )

    if count_difference != 0:
        raise ValueError(
            f"Sample '{sample_id}' failed primary-read count "
            f"consistency check: difference={count_difference}"
        )

    mapping_rate: float | None = None

    if primary_total > 0:
        mapping_rate = 100.0 * primary_mapped / primary_total

    mean_length, median_length, n50_length = (
        calculate_read_length_metrics(
            bam_path=bam_path,
            exclude_flags=primary_mapped_mask,
        )
    )

    stats = parse_samtools_stats(
        bam_path=bam_path,
        exclude_flags=panel_analysis_mask,
    )

    return {
        "sample_id": sample_id,
        "bam_path": str(bam_path),
        "bam_sort_order": sort_order,
        "total_alignment_records": total_records,
        "primary_total_reads": primary_total,
        "primary_mapped_reads": primary_mapped,
        "primary_unmapped_reads": primary_unmapped,
        "primary_mapping_rate_pct": mapping_rate,
        "secondary_alignment_records": secondary,
        "supplementary_alignment_records": supplementary,
        "duplicate_flagged_records": duplicate,
        "qc_failed_records": qc_failed,
        "panel_analysis_retained_records": (
            panel_analysis_records
        ),
        "mean_primary_mapped_read_length_bp": mean_length,
        "median_primary_mapped_read_length_bp": median_length,
        "n50_primary_mapped_read_length_bp": n50_length,
        "samtools_error_rate": stats["error_rate"],
        "mean_alignment_accuracy_pct": (
            stats["alignment_accuracy_pct"]
        ),
    }


def write_results(
    rows: list[dict[str, int | float | str | None]],
    output_path: Path,
) -> None:
    output_path.parent.mkdir(parents=True, exist_ok=True)

    fieldnames = [
        "sample_id",
        "bam_path",
        "bam_sort_order",
        "total_alignment_records",
        "primary_total_reads",
        "primary_mapped_reads",
        "primary_unmapped_reads",
        "primary_mapping_rate_pct",
        "secondary_alignment_records",
        "supplementary_alignment_records",
        "duplicate_flagged_records",
        "qc_failed_records",
        "panel_analysis_retained_records",
        "mean_primary_mapped_read_length_bp",
        "median_primary_mapped_read_length_bp",
        "n50_primary_mapped_read_length_bp",
        "samtools_error_rate",
        "mean_alignment_accuracy_pct",
    ]

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

        for row in rows:
            formatted_row = {
                key: format_value(value)
                for key, value in row.items()
            }
            writer.writerow(formatted_row)


def main() -> int:
    args = parse_arguments()

    try:
        (
            primary_total_mask,
            primary_mapped_mask,
            panel_analysis_mask,
        ) = read_filter_masks(args.filtering_config)

        samples = read_included_samples(args.samples)

        results = []

        for sample in samples:
            sample_id = sample["sample_id"].strip()
            print(
                f"[INFO] Calculating alignment QC: {sample_id}",
                flush=True,
            )

            metrics = calculate_sample_metrics(
                sample=sample,
                primary_total_mask=primary_total_mask,
                primary_mapped_mask=primary_mapped_mask,
                panel_analysis_mask=panel_analysis_mask,
            )

            results.append(metrics)

        write_results(results, args.output)

    except (ValueError, RuntimeError) as exc:
        print(f"[ERROR] {exc}", file=sys.stderr)
        return 1

    print(
        f"[INFO] Alignment QC completed for "
        f"{len(results)} sample(s)."
    )
    print(f"[INFO] Output: {args.output}")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
