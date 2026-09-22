#!/usr/bin/env python3

"""
Calculate target-region depth and panel-wide coverage metrics.

The script reads original wf-alignment BAM files directly. Alignment records
are filtered according to filtering.yaml and are streamed to samtools depth.
No permanent filtered BAM file is created.

Coordinate conventions:
    - BED coordinates are 0-based and half-open: [start, end)
    - samtools depth positions are 1-based
    - reported genomic coordinates are also provided as 1-based inclusive
"""

from __future__ import annotations

import argparse
import bisect
import csv
import math
import statistics
import subprocess
import sys
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

import yaml


TRUE_VALUES = {"true", "yes", "1"}


@dataclass
class TargetRegion:
    region_id: str
    chromosome: str
    start: int
    end: int
    depths: list[int] = field(default_factory=list)

    @property
    def length(self) -> int:
        return self.end - self.start


def parse_arguments() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Calculate target depth and coverage uniformity."
    )
    parser.add_argument(
        "--project-config",
        required=True,
        type=Path,
        help="Path to project.yaml",
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
        "--region-output",
        required=True,
        type=Path,
        help="Region-level output TSV",
    )
    parser.add_argument(
        "--panel-output",
        required=True,
        type=Path,
        help="Panel-level output TSV",
    )
    return parser.parse_args()


def load_yaml(path: Path) -> dict[str, Any]:
    if not path.is_file():
        raise ValueError(f"YAML file not found: {path}")

    with path.open("r", encoding="utf-8") as handle:
        data = yaml.safe_load(handle)

    if not isinstance(data, dict):
        raise ValueError(f"YAML root must be a mapping: {path}")

    return data


def read_configuration(
    project_config: Path,
    filtering_config: Path,
) -> tuple[Path, int, list[int], int, float]:
    project = load_yaml(project_config)
    filtering = load_yaml(filtering_config)

    try:
        target_bed = Path(
            str(project["resources"]["target_bed"])
        ).expanduser()

        coverage = project["coverage"]
        thresholds = [int(value) for value in coverage["thresholds"]]
        qc_depth = int(coverage["panel_qc_threshold_depth"])
        qc_required = float(
            coverage["panel_qc_required_percentage"]
        )

        exclude_flags = int(
            filtering["filters"]["panel_analysis"]["exclude_flags"]
        )
    except (KeyError, TypeError, ValueError) as exc:
        raise ValueError(
            "Invalid project or filtering configuration."
        ) from exc

    if not target_bed.is_file():
        raise ValueError(f"Target BED not found: {target_bed}")

    if not thresholds:
        raise ValueError("At least one coverage threshold is required.")

    if qc_depth not in thresholds:
        thresholds.append(qc_depth)

    thresholds = sorted(set(thresholds))

    return (
        target_bed,
        exclude_flags,
        thresholds,
        qc_depth,
        qc_required,
    )


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

        required = {"sample_id", "bam_path", "include"}

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

    return samples


def read_target_regions(target_bed: Path) -> list[TargetRegion]:
    regions: list[TargetRegion] = []

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

            if start < 0 or end <= start:
                raise ValueError(
                    f"Invalid BED interval on line {line_number}: "
                    f"{chromosome}:{start}-{end}"
                )

            if len(fields) >= 4 and fields[3].strip():
                region_id = fields[3].strip()
            else:
                region_id = (
                    f"region_{line_number:03d}_"
                    f"{chromosome}_{start}_{end}"
                )

            regions.append(
                TargetRegion(
                    region_id=region_id,
                    chromosome=chromosome,
                    start=start,
                    end=end,
                )
            )

    if not regions:
        raise ValueError("Target BED contains no valid intervals.")

    regions.sort(
        key=lambda region: (
            region.chromosome,
            region.start,
            region.end,
        )
    )

    previous: TargetRegion | None = None

    for region in regions:
        if (
            previous is not None
            and region.chromosome == previous.chromosome
            and region.start < previous.end
        ):
            raise ValueError(
                "Overlapping BED intervals are not supported by this "
                "region-level implementation. Merge the BED or resolve "
                f"the overlap between '{previous.region_id}' and "
                f"'{region.region_id}'."
            )

        previous = region

    return regions


def build_region_index(
    regions: list[TargetRegion],
) -> tuple[
    dict[str, list[TargetRegion]],
    dict[str, list[int]],
]:
    regions_by_chromosome: dict[str, list[TargetRegion]] = {}

    for region in regions:
        regions_by_chromosome.setdefault(
            region.chromosome,
            [],
        ).append(region)

    starts_by_chromosome = {
        chromosome: [region.start for region in chromosome_regions]
        for chromosome, chromosome_regions
        in regions_by_chromosome.items()
    }

    return regions_by_chromosome, starts_by_chromosome


def initialise_depth_arrays(regions: list[TargetRegion]) -> None:
    for region in regions:
        region.depths = [0] * region.length


def assign_depth(
    chromosome: str,
    position_1based: int,
    depth: int,
    regions_by_chromosome: dict[str, list[TargetRegion]],
    starts_by_chromosome: dict[str, list[int]],
) -> None:
    chromosome_regions = regions_by_chromosome.get(chromosome)

    if chromosome_regions is None:
        return

    position_0based = position_1based - 1
    starts = starts_by_chromosome[chromosome]

    index = bisect.bisect_right(starts, position_0based) - 1

    if index < 0:
        return

    region = chromosome_regions[index]

    if not (region.start <= position_0based < region.end):
        return

    offset = position_0based - region.start
    region.depths[offset] = depth


def stream_depth(
    bam_path: Path,
    target_bed: Path,
    exclude_flags: int,
    regions_by_chromosome: dict[str, list[TargetRegion]],
    starts_by_chromosome: dict[str, list[int]],
) -> None:
    command = [
        "samtools",
        "depth",
        "-aa",
        "-b",
        str(target_bed),
        "-G",
        str(exclude_flags),
        str(bam_path),
    ]

    process = subprocess.Popen(
        command,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True,
    )

    if process.stdout is None:
        raise RuntimeError(
            f"Unable to read samtools depth output for {bam_path}"
        )

    for line in process.stdout:
        fields = line.rstrip("\n").split("\t")

        if len(fields) < 3:
            continue

        chromosome = fields[0]

        try:
            position_1based = int(fields[1])
            depth = int(fields[2])
        except ValueError as exc:
            raise RuntimeError(
                f"Invalid samtools depth output: {line.rstrip()}"
            ) from exc

        assign_depth(
            chromosome=chromosome,
            position_1based=position_1based,
            depth=depth,
            regions_by_chromosome=regions_by_chromosome,
            starts_by_chromosome=starts_by_chromosome,
        )

    stderr = ""

    if process.stderr is not None:
        stderr = process.stderr.read().strip()

    return_code = process.wait()

    if return_code != 0:
        raise RuntimeError(
            f"samtools depth failed for {bam_path}: {stderr}"
        )


def calculate_depth_statistics(
    depths: list[int],
    thresholds: list[int],
) -> dict[str, int | float]:
    if not depths:
        raise ValueError("Depth array is empty.")

    number_of_bases = len(depths)
    mean_depth = statistics.fmean(depths)
    median_depth = float(statistics.median(depths))
    minimum_depth = min(depths)
    maximum_depth = max(depths)

    if number_of_bases > 1:
        standard_deviation = statistics.pstdev(depths)
    else:
        standard_deviation = 0.0

    coefficient_of_variation = (
        standard_deviation / mean_depth
        if mean_depth > 0
        else math.nan
    )

    metrics: dict[str, int | float] = {
        "target_bases": number_of_bases,
        "mean_depth": mean_depth,
        "median_depth": median_depth,
        "minimum_depth": minimum_depth,
        "maximum_depth": maximum_depth,
        "depth_standard_deviation": standard_deviation,
        "depth_coefficient_of_variation": (
            coefficient_of_variation
        ),
        "zero_depth_bases": sum(
            depth == 0
            for depth in depths
        ),
    }

    metrics["zero_depth_percentage"] = (
        100.0
        * int(metrics["zero_depth_bases"])
        / number_of_bases
    )

    for threshold in thresholds:
        covered_bases = sum(
            depth >= threshold
            for depth in depths
        )

        metrics[f"bases_ge_{threshold}x"] = covered_bases
        metrics[f"percentage_ge_{threshold}x"] = (
            100.0 * covered_bases / number_of_bases
        )

    return metrics


def format_value(value: int | float | str | bool) -> str:
    if isinstance(value, bool):
        return str(value).lower()

    if isinstance(value, float):
        if math.isnan(value):
            return "NA"
        return f"{value:.6f}"

    return str(value)


def analyse_sample(
    sample: dict[str, str],
    target_bed: Path,
    exclude_flags: int,
    thresholds: list[int],
    qc_depth: int,
    qc_required: float,
    region_template: list[TargetRegion],
) -> tuple[
    list[dict[str, int | float | str]],
    dict[str, int | float | str | bool],
]:
    sample_id = sample["sample_id"].strip()
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
        )
        raise ValueError(
            f"Sample '{sample_id}' failed BAM integrity check: "
            f"{message}"
        )

    regions = [
        TargetRegion(
            region_id=region.region_id,
            chromosome=region.chromosome,
            start=region.start,
            end=region.end,
        )
        for region in region_template
    ]

    initialise_depth_arrays(regions)

    (
        regions_by_chromosome,
        starts_by_chromosome,
    ) = build_region_index(regions)

    stream_depth(
        bam_path=bam_path,
        target_bed=target_bed,
        exclude_flags=exclude_flags,
        regions_by_chromosome=regions_by_chromosome,
        starts_by_chromosome=starts_by_chromosome,
    )

    region_results: list[
        dict[str, int | float | str]
    ] = []

    panel_depths: list[int] = []

    for region in regions:
        metrics = calculate_depth_statistics(
            depths=region.depths,
            thresholds=thresholds,
        )

        panel_depths.extend(region.depths)

        result: dict[str, int | float | str] = {
            "sample_id": sample_id,
            "region_id": region.region_id,
            "chromosome": region.chromosome,
            "bed_start_0based": region.start,
            "bed_end_0based_exclusive": region.end,
            "genomic_start_1based": region.start + 1,
            "genomic_end_1based_inclusive": region.end,
            "region_length_bp": region.length,
        }

        result.update(metrics)
        region_results.append(result)

    panel_metrics = calculate_depth_statistics(
        depths=panel_depths,
        thresholds=thresholds,
    )

    panel_qc_percentage = float(
        panel_metrics[f"percentage_ge_{qc_depth}x"]
    )

    panel_result: dict[
        str,
        int | float | str | bool,
    ] = {
        "sample_id": sample_id,
        "bam_path": str(bam_path),
        "panel_analysis_exclude_flags": exclude_flags,
        "panel_qc_depth_threshold": qc_depth,
        "panel_qc_required_percentage": qc_required,
        "panel_qc_pass": (
            panel_qc_percentage >= qc_required
        ),
    }

    panel_result.update(panel_metrics)

    return region_results, panel_result


def write_tsv(
    rows: list[dict[str, Any]],
    output_path: Path,
) -> None:
    if not rows:
        raise ValueError(
            f"No rows available for output: {output_path}"
        )

    output_path.parent.mkdir(parents=True, exist_ok=True)

    fieldnames = list(rows[0].keys())

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
            writer.writerow(
                {
                    key: format_value(value)
                    for key, value in row.items()
                }
            )


def main() -> int:
    args = parse_arguments()

    try:
        (
            target_bed,
            exclude_flags,
            thresholds,
            qc_depth,
            qc_required,
        ) = read_configuration(
            project_config=args.project_config,
            filtering_config=args.filtering_config,
        )

        samples = read_included_samples(args.samples)
        regions = read_target_regions(target_bed)

        total_target_bases = sum(
            region.length
            for region in regions
        )

        print(f"[INFO] Target regions: {len(regions)}")
        print(
            f"[INFO] Total target size: "
            f"{total_target_bases} bp"
        )
        print(
            f"[INFO] Coverage thresholds: "
            + ", ".join(f"{value}x" for value in thresholds)
        )
        print(
            f"[INFO] Panel QC criterion: "
            f">={qc_required:.2f}% of target bases "
            f">={qc_depth}x"
        )

        all_region_results: list[
            dict[str, int | float | str]
        ] = []

        all_panel_results: list[
            dict[str, int | float | str | bool]
        ] = []

        for sample in samples:
            sample_id = sample["sample_id"].strip()

            print(
                f"[INFO] Calculating target coverage: {sample_id}",
                flush=True,
            )

            region_results, panel_result = analyse_sample(
                sample=sample,
                target_bed=target_bed,
                exclude_flags=exclude_flags,
                thresholds=thresholds,
                qc_depth=qc_depth,
                qc_required=qc_required,
                region_template=regions,
            )

            all_region_results.extend(region_results)
            all_panel_results.append(panel_result)

        write_tsv(
            rows=all_region_results,
            output_path=args.region_output,
        )

        write_tsv(
            rows=all_panel_results,
            output_path=args.panel_output,
        )

    except (ValueError, RuntimeError) as exc:
        print(f"[ERROR] {exc}", file=sys.stderr)
        return 1

    print(
        f"[INFO] Coverage analysis completed for "
        f"{len(all_panel_results)} sample(s)."
    )
    print(f"[INFO] Region output: {args.region_output}")
    print(f"[INFO] Panel output: {args.panel_output}")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
