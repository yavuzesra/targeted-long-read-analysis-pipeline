#!/usr/bin/env python3

"""Compare relative gene-level coverage across two method-development runs."""

from __future__ import annotations

import argparse
import csv
import math
import statistics
from collections import defaultdict
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import TwoSlopeNorm


TRUE_VALUES = {"true", "yes", "1"}


def parse_arguments() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Normalise gene mean depth to each library's panel mean depth, "
            "rank genes, compare two runs, and create a combined heatmap."
        )
    )
    parser.add_argument("--samples", required=True, type=Path)
    parser.add_argument("--gene-coverage", required=True, type=Path)
    parser.add_argument("--panel-coverage", required=True, type=Path)
    parser.add_argument("--output-dir", required=True, type=Path)
    return parser.parse_args()


def read_tsv(path: Path) -> list[dict[str, str]]:
    if not path.is_file():
        raise ValueError(f"Input file not found: {path}")
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames is None:
            raise ValueError(f"TSV file has no header: {path}")
        rows = list(reader)
    if not rows:
        raise ValueError(f"TSV file contains no data rows: {path}")
    return rows


def write_tsv(
    path: Path,
    rows: list[dict[str, object]],
    fieldnames: list[str],
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


def method_sample_runs(rows: list[dict[str, str]]) -> dict[str, str]:
    selected: dict[str, str] = {}
    for row in rows:
        if (
            row.get("analysis_role", "").strip() != "method_development"
            or row.get("include", "").strip().lower() not in TRUE_VALUES
        ):
            continue
        sample_id = row.get("sample_id", "").strip()
        run_id = row.get("run_id", "").strip()
        if not sample_id or not run_id:
            raise ValueError(
                "Included method-development rows require sample_id and run_id."
            )
        if sample_id in selected:
            raise ValueError(f"Duplicate sample_id in metadata: {sample_id}")
        selected[sample_id] = run_id

    run_ids = sorted(set(selected.values()))
    if len(run_ids) != 2:
        raise ValueError(
            "Coverage reproducibility requires exactly two included "
            f"method-development runs; found {len(run_ids)}."
        )
    return selected


def panel_depths(
    rows: list[dict[str, str]],
    sample_runs: dict[str, str],
) -> dict[str, float]:
    result: dict[str, float] = {}
    for row in rows:
        sample_id = row.get("sample_id", "").strip()
        if sample_id not in sample_runs:
            continue
        try:
            value = float(row.get("mean_depth", "nan"))
        except ValueError as exc:
            raise ValueError(
                f"Invalid panel mean depth for {sample_id}."
            ) from exc
        if not math.isfinite(value) or value <= 0:
            raise ValueError(
                f"Panel mean depth must be positive for {sample_id}."
            )
        if sample_id in result:
            raise ValueError(f"Duplicate panel row for {sample_id}.")
        result[sample_id] = value

    missing = sorted(set(sample_runs) - set(result))
    if missing:
        raise ValueError("Missing panel rows: " + ", ".join(missing))
    return result


def average_ranks(values: list[float], reverse: bool = False) -> list[float]:
    order = sorted(
        range(len(values)), key=lambda index: values[index], reverse=reverse
    )
    ranks = [0.0] * len(values)
    position = 0
    while position < len(order):
        end = position + 1
        while (
            end < len(order)
            and values[order[end]] == values[order[position]]
        ):
            end += 1
        rank = (position + 1 + end) / 2.0
        for index in order[position:end]:
            ranks[index] = rank
        position = end
    return ranks


def pearson_correlation(x: list[float], y: list[float]) -> float:
    if len(x) != len(y) or len(x) < 2:
        raise ValueError("Correlation requires paired vectors of equal length.")
    mean_x = statistics.fmean(x)
    mean_y = statistics.fmean(y)
    numerator = sum(
        (value_x - mean_x) * (value_y - mean_y)
        for value_x, value_y in zip(x, y)
    )
    denominator = math.sqrt(
        sum((value - mean_x) ** 2 for value in x)
        * sum((value - mean_y) ** 2 for value in y)
    )
    return numerator / denominator if denominator else math.nan


def create_outputs(
    sample_rows: list[dict[str, str]],
    gene_rows: list[dict[str, str]],
    panel_rows: list[dict[str, str]],
) -> tuple[
    list[dict[str, object]],
    list[dict[str, object]],
    list[dict[str, object]],
    list[dict[str, object]],
]:
    sample_runs = method_sample_runs(sample_rows)
    panel = panel_depths(panel_rows, sample_runs)
    run_ids = sorted(set(sample_runs.values()))
    library_rows: list[dict[str, object]] = []
    grouped: dict[tuple[str, str, str], list[float]] = defaultdict(list)
    target_keys: dict[str, set[str]] = defaultdict(set)
    seen: set[tuple[str, str]] = set()

    for row in gene_rows:
        sample_id = row.get("sample_id", "").strip()
        if sample_id not in sample_runs:
            continue
        gene = row.get("gene_symbol", "").strip()
        target_key = row.get("target_key", "").strip()
        if not gene or not target_key:
            raise ValueError(
                f"Gene row lacks gene_symbol or target_key: {sample_id}"
            )
        key = (sample_id, gene)
        if key in seen:
            raise ValueError(
                "Expected one target row per sample and gene; duplicate: "
                f"{sample_id}, {gene}"
            )
        seen.add(key)
        try:
            gene_depth = float(row.get("mean_depth", "nan"))
        except ValueError as exc:
            raise ValueError(
                f"Invalid gene mean depth for {sample_id}, {gene}."
            ) from exc
        if not math.isfinite(gene_depth) or gene_depth < 0:
            raise ValueError(
                f"Invalid gene mean depth for {sample_id}, {gene}."
            )
        relative = gene_depth / panel[sample_id]
        run_id = sample_runs[sample_id]
        grouped[(run_id, gene, target_key)].append(relative)
        target_keys[gene].add(target_key)
        library_rows.append(
            {
                "sample_id": sample_id,
                "run_id": run_id,
                "gene_symbol": gene,
                "target_key": target_key,
                "gene_mean_depth": f"{gene_depth:.6f}",
                "panel_mean_depth": f"{panel[sample_id]:.6f}",
                "relative_gene_coverage": f"{relative:.6f}",
            }
        )

    expected_samples = set(sample_runs)
    genes_by_sample: dict[str, set[str]] = defaultdict(set)
    for row in library_rows:
        genes_by_sample[str(row["sample_id"])].add(str(row["gene_symbol"]))
    missing_samples = sorted(expected_samples - set(genes_by_sample))
    if missing_samples:
        raise ValueError(
            "Missing gene-coverage rows for: " + ", ".join(missing_samples)
        )
    reference_genes = genes_by_sample[next(iter(sorted(expected_samples)))]
    inconsistent = sorted(
        sample_id
        for sample_id, genes in genes_by_sample.items()
        if genes != reference_genes
    )
    if inconsistent:
        raise ValueError(
            "Method-development libraries have different gene sets: "
            + ", ".join(inconsistent)
        )

    summary_rows: list[dict[str, object]] = []
    for (run_id, gene, target_key), values in sorted(grouped.items()):
        summary_rows.append(
            {
                "run_id": run_id,
                "gene_symbol": gene,
                "target_key": target_key,
                "library_count": len(values),
                "mean_relative_gene_coverage": (
                    f"{statistics.fmean(values):.6f}"
                ),
                "sd_relative_gene_coverage": (
                    f"{statistics.stdev(values):.6f}"
                    if len(values) > 1
                    else "0.000000"
                ),
            }
        )

    for run_id in run_ids:
        selected = [row for row in summary_rows if row["run_id"] == run_id]
        ranks = average_ranks(
            [float(row["mean_relative_gene_coverage"]) for row in selected],
            reverse=True,
        )
        for row, rank in zip(selected, ranks):
            row["coverage_rank_high_to_low"] = f"{rank:.1f}"

    by_run_gene = {
        (str(row["run_id"]), str(row["gene_symbol"])): row
        for row in summary_rows
    }
    comparison_rows: list[dict[str, object]] = []
    for gene in sorted(reference_genes):
        first = by_run_gene[(run_ids[0], gene)]
        second = by_run_gene[(run_ids[1], gene)]
        first_value = float(first["mean_relative_gene_coverage"])
        second_value = float(second["mean_relative_gene_coverage"])
        first_rank = float(first["coverage_rank_high_to_low"])
        second_rank = float(second["coverage_rank_high_to_low"])
        comparison_rows.append(
            {
                "gene_symbol": gene,
                "target_key": ";".join(sorted(target_keys[gene])),
                "run_1_id": run_ids[0],
                "run_1_mean_relative_coverage": f"{first_value:.6f}",
                "run_1_rank_high_to_low": first["coverage_rank_high_to_low"],
                "run_2_id": run_ids[1],
                "run_2_mean_relative_coverage": f"{second_value:.6f}",
                "run_2_rank_high_to_low": second["coverage_rank_high_to_low"],
                "mean_relative_coverage_across_runs": (
                    f"{statistics.fmean((first_value, second_value)):.6f}"
                ),
                "absolute_rank_difference": (
                    f"{abs(first_rank - second_rank):.1f}"
                ),
            }
        )

    comparison_rows.sort(
        key=lambda row: (
            -float(row["mean_relative_coverage_across_runs"]),
            str(row["gene_symbol"]),
        )
    )
    first_vector = [
        float(row["run_1_mean_relative_coverage"])
        for row in comparison_rows
    ]
    second_vector = [
        float(row["run_2_mean_relative_coverage"])
        for row in comparison_rows
    ]
    spearman = pearson_correlation(
        average_ranks(first_vector),
        average_ranks(second_vector),
    )
    correlation_rows = [
        {
            "run_1_id": run_ids[0],
            "run_2_id": run_ids[1],
            "number_of_genes": len(comparison_rows),
            "pearson_r": f"{pearson_correlation(first_vector, second_vector):.6f}",
            "spearman_rho": f"{spearman:.6f}",
            "normalization": "gene_mean_depth/library_panel_mean_depth",
        }
    ]

    library_rows.sort(
        key=lambda row: (
            str(row["run_id"]),
            str(row["sample_id"]),
            str(row["gene_symbol"]),
        )
    )
    summary_rows.sort(
        key=lambda row: (
            str(row["run_id"]),
            float(row["coverage_rank_high_to_low"]),
            str(row["gene_symbol"]),
        )
    )
    return library_rows, summary_rows, comparison_rows, correlation_rows


def create_heatmap(
    library_rows: list[dict[str, object]],
    comparison_rows: list[dict[str, object]],
    output_dir: Path,
) -> None:
    sample_to_run = {
        str(row["sample_id"]): str(row["run_id"])
        for row in library_rows
    }
    run_ids = sorted(set(sample_to_run.values()))
    samples = sorted(sample_to_run, key=lambda item: (sample_to_run[item], item))
    genes = [str(row["gene_symbol"]) for row in comparison_rows]
    lookup = {
        (str(row["gene_symbol"]), str(row["sample_id"])): float(
            row["relative_gene_coverage"]
        )
        for row in library_rows
    }
    matrix = np.asarray(
        [[lookup[(gene, sample)] for sample in samples] for gene in genes],
        dtype=float,
    )
    positive = matrix[matrix > 0]
    if positive.size == 0:
        raise ValueError("Heatmap requires at least one positive coverage value.")
    zero_floor = float(np.min(positive)) / 2.0
    log_matrix = np.log2(np.where(matrix > 0, matrix, zero_floor))

    value_rows: list[dict[str, object]] = []
    for gene_index, gene in enumerate(genes):
        for sample_index, sample in enumerate(samples):
            value_rows.append(
                {
                    "gene_symbol": gene,
                    "sample_id": sample,
                    "run_id": sample_to_run[sample],
                    "relative_gene_coverage": f"{matrix[gene_index, sample_index]:.6f}",
                    "log2_relative_gene_coverage": (
                        f"{log_matrix[gene_index, sample_index]:.6f}"
                    ),
                    "zero_value_floor_applied": str(
                        matrix[gene_index, sample_index] == 0
                    ).lower(),
                }
            )
    write_tsv(
        output_dir / "gene_relative_coverage_heatmap_values.tsv",
        value_rows,
        [
            "gene_symbol",
            "sample_id",
            "run_id",
            "relative_gene_coverage",
            "log2_relative_gene_coverage",
            "zero_value_floor_applied",
        ],
    )

    lower_limit = min(float(np.min(log_matrix)), -0.25)
    upper_limit = max(float(np.max(log_matrix)), 0.25)
    figure_width = max(11.0, 0.62 * len(samples))
    figure_height = max(9.0, 0.30 * len(genes))
    figure, axis = plt.subplots(figsize=(figure_width, figure_height))
    image = axis.imshow(
        log_matrix,
        aspect="auto",
        cmap="RdBu_r",
        norm=TwoSlopeNorm(
            vmin=lower_limit,
            vcenter=0.0,
            vmax=upper_limit,
        ),
        interpolation="nearest",
    )
    axis.set_xticks(np.arange(len(samples)))
    axis.set_xticklabels(samples, rotation=55, ha="right", fontsize=8)
    axis.set_yticks(np.arange(len(genes)))
    axis.set_yticklabels(genes, fontsize=8)
    axis.set_xlabel("Method-development technical library")
    axis.set_ylabel("Panel gene")
    axis.set_title(
        "Relative gene-level coverage across method-development experiments",
        fontsize=13,
        fontweight="bold",
        pad=48,
    )

    technical_groups = [sample.rsplit("_", 1)[0] for sample in samples]
    group_segments: list[tuple[int, int, str]] = []
    segment_start = 0
    for index in range(1, len(technical_groups) + 1):
        if (
            index == len(technical_groups)
            or technical_groups[index] != technical_groups[segment_start]
        ):
            group_segments.append(
                (segment_start, index, technical_groups[segment_start])
            )
            segment_start = index

    for start, end, group in group_segments:
        if end < len(samples):
            axis.axvline(
                end - 0.5,
                color="#6F6F6F",
                linewidth=0.7,
                linestyle="--",
            )
        display_group = group.split("_", 1)[-1]
        axis.text(
            (start + end - 1) / 2,
            1.005,
            display_group,
            transform=axis.get_xaxis_transform(),
            ha="center",
            va="bottom",
            fontsize=8,
            clip_on=False,
        )

    first_run_count = sum(
        sample_to_run[sample] == run_ids[0] for sample in samples
    )
    if 0 < first_run_count < len(samples):
        axis.axvline(first_run_count - 0.5, color="black", linewidth=1.2)
        first_midpoint = (first_run_count - 1) / 2
        second_midpoint = first_run_count + (len(samples) - first_run_count - 1) / 2
        axis.text(
            first_midpoint,
            1.045,
            run_ids[0],
            transform=axis.get_xaxis_transform(),
            ha="center",
            va="bottom",
            fontsize=10,
            fontweight="bold",
            clip_on=False,
        )
        axis.text(
            second_midpoint,
            1.045,
            run_ids[1],
            transform=axis.get_xaxis_transform(),
            ha="center",
            va="bottom",
            fontsize=10,
            fontweight="bold",
            clip_on=False,
        )

    colour_bar = figure.colorbar(image, ax=axis, pad=0.02)
    colour_bar.set_label(
        "log2(gene mean depth / library panel mean depth)"
    )
    figure.tight_layout()
    figure.savefig(
        output_dir / "gene_relative_coverage_heatmap.png",
        dpi=300,
        bbox_inches="tight",
    )
    plt.close(figure)


def main() -> int:
    args = parse_arguments()
    output_dir = args.output_dir
    output_dir.mkdir(parents=True, exist_ok=True)
    outputs = create_outputs(
        read_tsv(args.samples),
        read_tsv(args.gene_coverage),
        read_tsv(args.panel_coverage),
    )
    specifications = [
        (
            "gene_relative_coverage_per_library.tsv",
            [
                "sample_id",
                "run_id",
                "gene_symbol",
                "target_key",
                "gene_mean_depth",
                "panel_mean_depth",
                "relative_gene_coverage",
            ],
        ),
        (
            "gene_relative_coverage_experiment_summary.tsv",
            [
                "run_id",
                "gene_symbol",
                "target_key",
                "library_count",
                "mean_relative_gene_coverage",
                "sd_relative_gene_coverage",
                "coverage_rank_high_to_low",
            ],
        ),
        (
            "gene_relative_coverage_experiment_comparison.tsv",
            [
                "gene_symbol",
                "target_key",
                "run_1_id",
                "run_1_mean_relative_coverage",
                "run_1_rank_high_to_low",
                "run_2_id",
                "run_2_mean_relative_coverage",
                "run_2_rank_high_to_low",
                "mean_relative_coverage_across_runs",
                "absolute_rank_difference",
            ],
        ),
        (
            "gene_relative_coverage_correlation.tsv",
            [
                "run_1_id",
                "run_2_id",
                "number_of_genes",
                "pearson_r",
                "spearman_rho",
                "normalization",
            ],
        ),
    ]
    for rows, (filename, fieldnames) in zip(outputs, specifications):
        path = output_dir / filename
        write_tsv(path, rows, fieldnames)
        print(f"[INFO] Output: {path}")
    create_heatmap(outputs[0], outputs[2], output_dir)
    print(f"[INFO] Output: {output_dir / 'gene_relative_coverage_heatmap.png'}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
