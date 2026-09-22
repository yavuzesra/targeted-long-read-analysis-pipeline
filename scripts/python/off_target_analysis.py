#!/usr/bin/env python3
"""
Perform an exploratory base-level off-target analysis.

The script measures reference-aligned bases retained by the panel-analysis
filter outside the target BED intervals. It does not provide a formal QC
acceptance criterion.

Retained alignments exclude:
- unmapped records
- secondary alignments
- QC-failed records
- duplicate records
- supplementary alignments

Off-target bases are calculated from bedGraph coverage after subtracting
the merged target BED intervals. The bedtools genomecov --split option
excludes reference deletions and counts aligned reference blocks.
"""

from __future__ import annotations

import argparse
import csv
import subprocess
import sys
from collections import defaultdict
from pathlib import Path
from typing import TextIO


TRUE_VALUES = {"true", "yes", "1"}
PANEL_EXCLUDE_FLAGS = 3844


def parse_arguments() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Calculate exploratory base-level off-target distributions "
            "for selected samples."
        )
    )
    parser.add_argument(
        "--samples",
        required=True,
        type=Path,
        help="Selected samples.tsv file.",
    )
    parser.add_argument(
        "--target-bed",
        required=True,
        type=Path,
        help="Custom target BED file.",
    )
    parser.add_argument(
        "--reference-fai",
        required=True,
        type=Path,
        help="Reference FASTA index file.",
    )
    parser.add_argument(
        "--work-root",
        required=True,
        type=Path,
        help="Directory for intermediate files.",
    )
    parser.add_argument(
        "--output-root",
        required=True,
        type=Path,
        help="Directory for final result tables.",
    )
    parser.add_argument(
        "--window-size",
        type=int,
        default=1_000_000,
        help="Genomic window size in base pairs. Default: 1,000,000.",
    )
    parser.add_argument(
        "--top-windows",
        type=int,
        default=50,
        help="Number of highest-density windows to report. Default: 50.",
    )
    return parser.parse_args()


def run_command(
    command: list[str],
    *,
    stdout_path: Path | None = None,
) -> None:
    if stdout_path is None:
        subprocess.run(command, check=True)
        return

    stdout_path.parent.mkdir(parents=True, exist_ok=True)

    with stdout_path.open("w", encoding="utf-8") as handle:
        subprocess.run(
            command,
            check=True,
            stdout=handle,
            text=True,
        )


def run_pipeline_to_file(
    producer_command: list[str],
    consumer_command: list[str],
    output_path: Path,
) -> None:
    """Pipe binary output between two commands and write the final output."""
    output_path.parent.mkdir(parents=True, exist_ok=True)

    with output_path.open("wb") as handle:
        producer = subprocess.Popen(
            producer_command,
            stdout=subprocess.PIPE,
        )

        if producer.stdout is None:
            producer.kill()
            raise RuntimeError(
                f"Could not capture output from: {' '.join(producer_command)}"
            )

        consumer = subprocess.Popen(
            consumer_command,
            stdin=producer.stdout,
            stdout=handle,
        )

        producer.stdout.close()

        consumer_returncode = consumer.wait()
        producer_returncode = producer.wait()

    if producer_returncode != 0:
        raise subprocess.CalledProcessError(
            producer_returncode,
            producer_command,
        )

    if consumer_returncode != 0:
        raise subprocess.CalledProcessError(
            consumer_returncode,
            consumer_command,
        )


def require_file(path: Path, label: str) -> None:
    if not path.is_file():
        raise FileNotFoundError(f"{label} not found: {path}")

    if path.stat().st_size == 0:
        raise ValueError(f"{label} is empty: {path}")


def is_included(value: str) -> bool:
    return value.strip().lower() in TRUE_VALUES


def read_samples(path: Path) -> list[dict[str, str]]:
    require_file(path, "Samples metadata")

    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")

        required = {
            "sample_id",
            "bam_path",
            "bai_path",
            "include",
        }

        if reader.fieldnames is None:
            raise ValueError("Samples metadata has no header.")

        missing = required - set(reader.fieldnames)

        if missing:
            raise ValueError(
                "Samples metadata is missing required columns: "
                + ", ".join(sorted(missing))
            )

        samples = [
            row
            for row in reader
            if is_included(row.get("include", "false"))
        ]

    if not samples:
        raise ValueError("No included samples were found.")

    return samples


def write_genome_sizes(
    reference_fai: Path,
    output_path: Path,
) -> dict[str, int]:
    chromosome_lengths: dict[str, int] = {}

    with reference_fai.open("r", encoding="utf-8") as source:
        with output_path.open("w", encoding="utf-8") as target:
            for line in source:
                fields = line.rstrip("\n").split("\t")

                if len(fields) < 2:
                    continue

                chromosome = fields[0]
                length = int(fields[1])

                chromosome_lengths[chromosome] = length
                target.write(f"{chromosome}\t{length}\n")

    if not chromosome_lengths:
        raise ValueError(
            f"No chromosome lengths were read from {reference_fai}"
        )

    return chromosome_lengths


def bedgraph_total(path: Path) -> int:
    total = 0

    with path.open("r", encoding="utf-8") as handle:
        for line in handle:
            chromosome, start, end, depth = line.rstrip("\n").split("\t")[:4]

            del chromosome

            interval_length = int(end) - int(start)
            total += interval_length * int(float(depth))

    return total


def chromosome_bases(path: Path) -> dict[str, int]:
    totals: dict[str, int] = defaultdict(int)

    with path.open("r", encoding="utf-8") as handle:
        for line in handle:
            chromosome, start, end, depth = line.rstrip("\n").split("\t")[:4]

            interval_length = int(end) - int(start)
            totals[chromosome] += interval_length * int(float(depth))

    return dict(totals)


def calculate_window_bases(
    intersect_path: Path,
) -> dict[tuple[str, int, int], int]:
    window_totals: dict[tuple[str, int, int], int] = defaultdict(int)

    with intersect_path.open("r", encoding="utf-8") as handle:
        for line in handle:
            fields = line.rstrip("\n").split("\t")

            chromosome = fields[0]
            start = int(fields[1])
            end = int(fields[2])
            depth_text = fields[6]
            overlap = int(fields[7])

            if depth_text == "." or overlap <= 0:
                window_totals[(chromosome, start, end)] += 0
                continue

            depth = int(float(depth_text))
            window_totals[(chromosome, start, end)] += depth * overlap

    return dict(window_totals)


def summarise_base_counts(
    total_retained_aligned_bases: int,
    off_target_aligned_bases: int,
) -> dict[str, int | float]:
    """Calculate mutually consistent on-target and off-target totals."""
    if total_retained_aligned_bases < 0:
        raise ValueError(
            "Total retained aligned bases cannot be negative."
        )

    if off_target_aligned_bases < 0:
        raise ValueError(
            "Off-target aligned bases cannot be negative."
        )

    if off_target_aligned_bases > total_retained_aligned_bases:
        raise ValueError(
            "Off-target aligned bases exceed total retained aligned bases."
        )

    on_target_aligned_bases = (
        total_retained_aligned_bases - off_target_aligned_bases
    )

    if total_retained_aligned_bases > 0:
        on_target_fraction_percent = (
            100.0
            * on_target_aligned_bases
            / total_retained_aligned_bases
        )
        off_target_fraction_percent = (
            100.0 - on_target_fraction_percent
        )
    else:
        on_target_fraction_percent = 0.0
        off_target_fraction_percent = 0.0

    return {
        "total_retained_aligned_bases": total_retained_aligned_bases,
        "on_target_aligned_bases": on_target_aligned_bases,
        "off_target_aligned_bases": off_target_aligned_bases,
        "on_target_fraction_percent": on_target_fraction_percent,
        "off_target_fraction_percent": off_target_fraction_percent,
    }


def write_tsv(
    path: Path,
    fieldnames: list[str],
    rows: list[dict[str, object]],
) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)

    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=fieldnames,
            delimiter="\t",
            lineterminator="\n",
        )
        writer.writeheader()
        writer.writerows(rows)


def main() -> int:
    args = parse_arguments()

    require_file(args.target_bed, "Target BED")
    require_file(args.reference_fai, "Reference FASTA index")

    if args.window_size <= 0:
        raise ValueError("--window-size must be greater than zero.")

    if args.top_windows <= 0:
        raise ValueError("--top-windows must be greater than zero.")

    samples = read_samples(args.samples)

    args.work_root.mkdir(parents=True, exist_ok=True)
    args.output_root.mkdir(parents=True, exist_ok=True)

    shared_work = args.work_root / "shared"
    shared_work.mkdir(parents=True, exist_ok=True)

    genome_sizes = shared_work / "genome.sizes"
    merged_target = shared_work / "target_merged.bed"
    sorted_target = shared_work / "target_sorted.bed"
    genome_windows = shared_work / (
        f"genome_{args.window_size}_bp_windows.bed"
    )

    chromosome_lengths = write_genome_sizes(
        args.reference_fai,
        genome_sizes,
    )

    run_command(
        [
            "bedtools",
            "sort",
            "-i",
            str(args.target_bed),
        ],
        stdout_path=sorted_target,
    )

    run_command(
        [
            "bedtools",
            "merge",
            "-i",
            str(sorted_target),
        ],
        stdout_path=merged_target,
    )

    run_command(
        [
            "bedtools",
            "makewindows",
            "-g",
            str(genome_sizes),
            "-w",
            str(args.window_size),
        ],
        stdout_path=genome_windows,
    )

    sample_summary_rows: list[dict[str, object]] = []
    chromosome_rows: list[dict[str, object]] = []
    window_rows: list[dict[str, object]] = []

    for sample in samples:
        sample_id = sample["sample_id"].strip()
        bam_path = Path(sample["bam_path"].strip())
        bai_path = Path(sample["bai_path"].strip())

        require_file(bam_path, f"{sample_id} BAM")
        require_file(bai_path, f"{sample_id} BAM index")

        print(f"[INFO] Processing off-target bases: {sample_id}")

        sample_work = args.work_root / sample_id
        sample_work.mkdir(parents=True, exist_ok=True)

        retained_bedgraph = sample_work / (
            f"{sample_id}.retained.coverage.bedgraph"
        )
        off_target_bedgraph = sample_work / (
            f"{sample_id}.off_target.coverage.bedgraph"
        )
        window_intersections = sample_work / (
            f"{sample_id}.off_target.window_intersections.tsv"
        )

        run_pipeline_to_file(
            [
                "samtools",
                "view",
                "-u",
                "-F",
                str(PANEL_EXCLUDE_FLAGS),
                str(bam_path),
            ],
            [
                "bedtools",
                "genomecov",
                "-ibam",
                "stdin",
                "-bg",
                "-split",
            ],
            retained_bedgraph,
        )

        run_command(
            [
                "bedtools",
                "subtract",
                "-a",
                str(retained_bedgraph),
                "-b",
                str(merged_target),
            ],
            stdout_path=off_target_bedgraph,
        )

        run_command(
            [
                "bedtools",
                "intersect",
                "-a",
                str(genome_windows),
                "-b",
                str(off_target_bedgraph),
                "-wao",
            ],
            stdout_path=window_intersections,
        )

        total_retained_aligned_bases = bedgraph_total(
            retained_bedgraph
        )
        off_target_aligned_bases = bedgraph_total(
            off_target_bedgraph
        )

        base_summary = summarise_base_counts(
            total_retained_aligned_bases,
            off_target_aligned_bases,
        )

        on_target_aligned_bases = int(
            base_summary["on_target_aligned_bases"]
        )

        sample_summary_rows.append(
            {
                "sample_id": sample_id,
                "total_retained_aligned_bases": (
                    total_retained_aligned_bases
                ),
                "on_target_aligned_bases": on_target_aligned_bases,
                "off_target_aligned_bases": off_target_aligned_bases,
                "on_target_fraction_percent": (
                    f"{base_summary['on_target_fraction_percent']:.6f}"
                ),
                "off_target_fraction_percent": (
                    f"{base_summary['off_target_fraction_percent']:.6f}"
                ),
            }
        )

        per_chromosome = chromosome_bases(off_target_bedgraph)

        for chromosome, chromosome_length in chromosome_lengths.items():
            off_target_bases = per_chromosome.get(chromosome, 0)

            mean_depth = (
                off_target_bases / chromosome_length
                if chromosome_length > 0
                else 0.0
            )

            aligned_bases_per_mb = mean_depth * 1_000_000

            fraction_of_sample = (
                off_target_bases / off_target_aligned_bases
                if off_target_aligned_bases > 0
                else 0.0
            )

            normalised_per_million_primary = (
                off_target_bases
                / total_retained_aligned_bases
                * 1_000_000
                if total_retained_aligned_bases > 0
                else 0.0
            )

            chromosome_rows.append(
                {
                    "sample_id": sample_id,
                    "chromosome": chromosome,
                    "chromosome_length_bp": chromosome_length,
                    "off_target_aligned_bases": off_target_bases,
                    "off_target_aligned_bases_per_mb": (
                        f"{aligned_bases_per_mb:.6f}"
                    ),
                    "mean_off_target_depth": f"{mean_depth:.9f}",
                    "fraction_of_sample_off_target_bases": (
                        f"{fraction_of_sample:.9f}"
                    ),
                    "off_target_bases_per_million_total_retained_aligned_bases": (
                        f"{normalised_per_million_primary:.6f}"
                    ),
                }
            )

        sample_windows = calculate_window_bases(window_intersections)

        for (chromosome, start, end), off_target_bases in sample_windows.items():
            window_length = end - start

            mean_depth = (
                off_target_bases / window_length
                if window_length > 0
                else 0.0
            )

            fraction_of_sample = (
                off_target_bases / off_target_aligned_bases
                if off_target_aligned_bases > 0
                else 0.0
            )

            normalised_per_million_primary = (
                off_target_bases
                / total_retained_aligned_bases
                * 1_000_000
                if total_retained_aligned_bases > 0
                else 0.0
            )

            window_rows.append(
                {
                    "sample_id": sample_id,
                    "chromosome": chromosome,
                    "window_start_0based": start,
                    "window_end_0based": end,
                    "window_length_bp": window_length,
                    "off_target_aligned_bases": off_target_bases,
                    "mean_off_target_depth": f"{mean_depth:.9f}",
                    "fraction_of_sample_off_target_bases": (
                        f"{fraction_of_sample:.9f}"
                    ),
                    "off_target_bases_per_million_total_retained_aligned_bases": (
                        f"{normalised_per_million_primary:.6f}"
                    ),
                }
            )

    sample_summary_path = (
        args.output_root / "off_target_base_summary_per_sample.tsv"
    )
    chromosome_path = (
        args.output_root
        / "off_target_aligned_bases_per_chromosome_per_sample.tsv"
    )
    windows_path = (
        args.output_root
        / "off_target_aligned_bases_1Mb_windows_per_sample.tsv"
    )

    write_tsv(
        sample_summary_path,
        [
            "sample_id",
            "total_retained_aligned_bases",
            "on_target_aligned_bases",
            "off_target_aligned_bases",
            "on_target_fraction_percent",
            "off_target_fraction_percent",
        ],
        sample_summary_rows,
    )

    write_tsv(
        chromosome_path,
        [
            "sample_id",
            "chromosome",
            "chromosome_length_bp",
            "off_target_aligned_bases",
            "off_target_aligned_bases_per_mb",
            "mean_off_target_depth",
            "fraction_of_sample_off_target_bases",
            "off_target_bases_per_million_total_retained_aligned_bases",
        ],
        chromosome_rows,
    )

    write_tsv(
        windows_path,
        [
            "sample_id",
            "chromosome",
            "window_start_0based",
            "window_end_0based",
            "window_length_bp",
            "off_target_aligned_bases",
            "mean_off_target_depth",
            "fraction_of_sample_off_target_bases",
            "off_target_bases_per_million_total_retained_aligned_bases",
        ],
        window_rows,
    )

    hotspot_groups: dict[
        tuple[str, int, int, int],
        list[dict[str, object]],
    ] = defaultdict(list)

    for row in window_rows:
        key = (
            str(row["chromosome"]),
            int(row["window_start_0based"]),
            int(row["window_end_0based"]),
            int(row["window_length_bp"]),
        )
        hotspot_groups[key].append(row)

    hotspot_rows: list[dict[str, object]] = []

    for key, rows in hotspot_groups.items():
        chromosome, start, end, window_length = key

        total_bases = sum(
            int(row["off_target_aligned_bases"])
            for row in rows
        )

        mean_bases = total_bases / len(rows)

        mean_depth = sum(
            float(row["mean_off_target_depth"])
            for row in rows
        ) / len(rows)

        mean_normalised_density = sum(
            float(
                row[
                    "off_target_bases_per_million_total_retained_aligned_bases"
                ]
            )
            for row in rows
        ) / len(rows)

        samples_with_signal = sum(
            int(row["off_target_aligned_bases"]) > 0
            for row in rows
        )

        hotspot_rows.append(
            {
                "chromosome": chromosome,
                "window_start_0based": start,
                "window_end_0based": end,
                "window_length_bp": window_length,
                "sample_count": len(rows),
                "samples_with_off_target_signal": samples_with_signal,
                "total_off_target_aligned_bases": total_bases,
                "mean_off_target_aligned_bases": (
                    f"{mean_bases:.6f}"
                ),
                "mean_off_target_depth": f"{mean_depth:.9f}",
                "mean_off_target_bases_per_million_total_retained_aligned_bases": (
                    f"{mean_normalised_density:.6f}"
                ),
            }
        )

    hotspot_rows.sort(
        key=lambda row: float(
            row[
                "mean_off_target_bases_per_million_total_retained_aligned_bases"
            ]
        ),
        reverse=True,
    )

    hotspot_path = (
        args.output_root / "off_target_1Mb_hotspots_overall.tsv"
    )
    top_hotspot_path = (
        args.output_root
        / f"top{args.top_windows}_off_target_1Mb_hotspots.tsv"
    )

    hotspot_fields = [
        "chromosome",
        "window_start_0based",
        "window_end_0based",
        "window_length_bp",
        "sample_count",
        "samples_with_off_target_signal",
        "total_off_target_aligned_bases",
        "mean_off_target_aligned_bases",
        "mean_off_target_depth",
        "mean_off_target_bases_per_million_total_retained_aligned_bases",
    ]

    write_tsv(
        hotspot_path,
        hotspot_fields,
        hotspot_rows,
    )

    write_tsv(
        top_hotspot_path,
        hotspot_fields,
        hotspot_rows[: args.top_windows],
    )

    print()
    print(f"[INFO] Off-target samples analysed: {len(samples)}")
    print(f"[INFO] Sample summary: {sample_summary_path}")
    print(f"[INFO] Chromosome summary: {chromosome_path}")
    print(f"[INFO] Window table: {windows_path}")
    print(f"[INFO] Hotspot table: {hotspot_path}")
    print(f"[INFO] Top hotspots: {top_hotspot_path}")

    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except (
        FileNotFoundError,
        ValueError,
        subprocess.CalledProcessError,
    ) as error:
        print(f"ERROR: {error}", file=sys.stderr)
        raise SystemExit(1)
