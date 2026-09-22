#!/usr/bin/env python3

"""Compare mean base-wise depth across two method-development runs."""

from __future__ import annotations

import argparse
import csv
import math
import shutil
import subprocess
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import yaml


TRUE_VALUES = {"true", "yes", "1"}
EMPTY_VALUES = {"", "NA", "N/A", "none", "None", "."}


def parse_arguments() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Calculate base-wise depth reproducibility for every panel gene "
            "and plot objectively selected high- and low-coverage genes."
        )
    )
    parser.add_argument("--samples", required=True, type=Path)
    parser.add_argument("--gene-coverage", required=True, type=Path)
    parser.add_argument("--gene-comparison", required=True, type=Path)
    parser.add_argument("--filtering-config", required=True, type=Path)
    parser.add_argument("--output-dir", required=True, type=Path)
    parser.add_argument("--samtools", default="samtools")
    parser.add_argument("--top-n", type=int, default=3)
    parser.add_argument(
        "--maximum-plotted-points",
        type=int,
        default=50000,
        help=(
            "Maximum deterministic points drawn per gene. All target bases "
            "remain included in the statistics. Applies only when "
            "--plot-style scatter is selected."
        ),
    )
    parser.add_argument(
        "--plot-style",
        choices=("density", "scatter"),
        default="density",
        help=(
            "Figure representation. Density is the default and includes all "
            "target bases in logarithmically coloured hexagonal bins. Scatter "
            "retains the previous deterministic-subset rendering for optional "
            "debugging."
        ),
    )
    parser.add_argument(
        "--density-gridsize",
        type=int,
        default=60,
        help="Number of hexagonal bins along the x axis for density figures.",
    )
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
    if len(set(selected.values())) != 2:
        raise ValueError(
            "Base-wise reproducibility requires exactly two included "
            "method-development runs."
        )
    return selected


def bam_paths(
    sample_rows: list[dict[str, str]],
    sample_runs: dict[str, str],
) -> dict[str, Path]:
    result: dict[str, Path] = {}
    for row in sample_rows:
        sample_id = row.get("sample_id", "").strip()
        if sample_id not in sample_runs:
            continue
        raw_path = row.get("bam_path", "").strip()
        if raw_path in EMPTY_VALUES:
            raise ValueError(f"Missing bam_path for {sample_id}.")
        path = Path(raw_path).expanduser()
        if not path.is_file():
            raise ValueError(f"BAM file not found for {sample_id}: {path}")

        configured_index = row.get("bai_path", "").strip()
        index_candidates = []
        if configured_index not in EMPTY_VALUES:
            index_candidates.append(Path(configured_index).expanduser())
        index_candidates.extend(
            [Path(f"{path}.bai"), path.with_suffix(".bai"), Path(f"{path}.csi")]
        )
        if not any(candidate.is_file() for candidate in index_candidates):
            raise ValueError(f"BAM index not found for {sample_id}: {path}")
        result[sample_id] = path

    missing = sorted(set(sample_runs) - set(result))
    if missing:
        raise ValueError("Missing BAM rows for: " + ", ".join(missing))
    return result


def panel_exclude_flags(path: Path) -> int:
    if not path.is_file():
        raise ValueError(f"Filtering configuration not found: {path}")
    with path.open("r", encoding="utf-8") as handle:
        config = yaml.safe_load(handle) or {}
    try:
        value = int(config["filters"]["panel_analysis"]["exclude_flags"])
    except (KeyError, TypeError, ValueError) as exc:
        raise ValueError(
            "filtering.yaml must define filters.panel_analysis.exclude_flags."
        ) from exc
    if value < 0:
        raise ValueError("Panel-analysis exclude_flags cannot be negative.")
    return value


def target_gene_mapping(
    rows: list[dict[str, str]],
    sample_runs: dict[str, str],
) -> dict[str, dict[str, object]]:
    targets: dict[str, dict[str, object]] = {}
    for row in rows:
        if row.get("sample_id", "").strip() not in sample_runs:
            continue
        target_key = row.get("target_key", "").strip()
        gene = row.get("gene_symbol", "").strip()
        try:
            coordinates = (
                row.get("chromosome", "").strip(),
                int(row.get("genomic_start_1based", "0")),
                int(row.get("genomic_end_1based_inclusive", "0")),
            )
        except ValueError as exc:
            raise ValueError(
                f"Invalid target coordinates for {target_key}."
            ) from exc
        if (
            not target_key
            or not gene
            or not coordinates[0]
            or coordinates[1] < 1
            or coordinates[2] < coordinates[1]
        ):
            raise ValueError(f"Invalid target-gene row: {target_key}, {gene}")
        if target_key not in targets:
            targets[target_key] = {
                "chromosome": coordinates[0],
                "start": coordinates[1],
                "end": coordinates[2],
                "genes": set(),
            }
        target = targets[target_key]
        if (target["chromosome"], target["start"], target["end"]) != coordinates:
            raise ValueError(f"Inconsistent coordinates for {target_key}.")
        genes = target["genes"]
        assert isinstance(genes, set)
        genes.add(gene)
    if not targets:
        raise ValueError("No method-development target-gene mappings found.")
    return targets


def experiment_sample_order(
    sample_runs: dict[str, str],
) -> tuple[list[str], list[str], list[str]]:
    run_ids = sorted(set(sample_runs.values()))
    first = sorted(
        sample for sample, run in sample_runs.items() if run == run_ids[0]
    )
    second = sorted(
        sample for sample, run in sample_runs.items() if run == run_ids[1]
    )
    if not first or not second:
        raise ValueError("Each method-development run requires libraries.")
    return run_ids, first, second


def rankdata(values: np.ndarray) -> np.ndarray:
    order = np.argsort(values, kind="mergesort")
    ranks = np.empty(len(values), dtype=float)
    position = 0
    while position < len(order):
        end = position + 1
        while (
            end < len(order)
            and values[order[end]] == values[order[position]]
        ):
            end += 1
        ranks[order[position:end]] = (position + 1 + end) / 2.0
        position = end
    return ranks


def correlation(x: np.ndarray, y: np.ndarray) -> float:
    if len(x) < 2 or np.ptp(x) == 0 or np.ptp(y) == 0:
        return math.nan
    return float(np.corrcoef(x, y)[0, 1])


def regression_statistics(x: np.ndarray, y: np.ndarray) -> dict[str, float]:
    pearson = correlation(x, y)
    spearman = correlation(rankdata(x), rankdata(y))
    if np.ptp(x) == 0:
        slope = math.nan
        intercept = math.nan
    else:
        slope, intercept = np.polyfit(x, y, 1)
    return {
        "pearson_r": pearson,
        "spearman_rho": spearman,
        "linear_regression_r_squared": (
            pearson**2 if math.isfinite(pearson) else math.nan
        ),
        "linear_regression_slope": float(slope),
        "linear_regression_intercept": float(intercept),
    }


def run_samtools_depth(
    samtools: str,
    exclude_flags: int,
    chromosome: str,
    start: int,
    end: int,
    ordered_paths: list[Path],
) -> np.ndarray:
    region = f"{chromosome}:{start}-{end}"
    command = [
        samtools,
        "depth",
        "-aa",
        "-d",
        "0",
        "-G",
        str(exclude_flags),
        "-r",
        region,
        *[str(path) for path in ordered_paths],
    ]
    completed = subprocess.run(
        command,
        check=False,
        capture_output=True,
        text=True,
    )
    if completed.returncode != 0:
        raise RuntimeError(
            f"samtools depth failed for {region}: {completed.stderr.strip()}"
        )

    expected_positions = end - start + 1
    positions: list[int] = []
    depths: list[list[float]] = []
    for line in completed.stdout.splitlines():
        fields = line.split("\t")
        if len(fields) != 2 + len(ordered_paths):
            raise ValueError(
                f"Unexpected samtools depth column count for {region}."
            )
        positions.append(int(fields[1]))
        depths.append([float(value) for value in fields[2:]])

    if len(positions) != expected_positions:
        raise ValueError(
            f"Expected {expected_positions} positions for {region}, "
            f"received {len(positions)}."
        )
    expected = np.arange(start, end + 1)
    if not np.array_equal(np.asarray(positions), expected):
        raise ValueError(f"Non-contiguous depth positions for {region}.")
    return np.asarray(depths, dtype=float)


def deterministic_plot_subset(
    x: np.ndarray,
    y: np.ndarray,
    maximum_points: int,
) -> tuple[np.ndarray, np.ndarray]:
    if maximum_points < 1:
        raise ValueError("maximum-plotted-points must be positive.")
    if len(x) <= maximum_points:
        return x, y
    indices = np.linspace(0, len(x) - 1, maximum_points, dtype=int)
    return x[indices], y[indices]


def selected_genes(
    rows: list[dict[str, str]],
    top_n: int,
) -> tuple[list[str], list[str]]:
    if top_n < 1:
        raise ValueError("top-n must be positive.")
    ranked = sorted(
        rows,
        key=lambda row: (
            -float(row["mean_relative_coverage_across_runs"]),
            row["gene_symbol"],
        ),
    )
    if len(ranked) < top_n * 2:
        raise ValueError("Not enough genes for disjoint high and low panels.")
    high = [row["gene_symbol"] for row in ranked[:top_n]]
    low = [row["gene_symbol"] for row in ranked[-top_n:]][::-1]
    return high, low


def decorate_gene_panel(
    axis: plt.Axes,
    gene: str,
    x: np.ndarray,
    y: np.ndarray,
    stats: dict[str, float],
    run_ids: list[str],
) -> None:
    slope = stats["linear_regression_slope"]
    intercept = stats["linear_regression_intercept"]
    if math.isfinite(slope) and math.isfinite(intercept):
        line_x = np.asarray([float(np.min(x)), float(np.max(x))])
        axis.plot(
            line_x,
            slope * line_x + intercept,
            color="#D95F02",
            linewidth=1.2,
        )
    axis.set_title(gene, fontsize=11, fontweight="bold")
    axis.set_xlabel(f"{run_ids[0]} mean depth per base", fontsize=9)
    axis.set_ylabel(f"{run_ids[1]} mean depth per base", fontsize=9)
    axis.grid(alpha=0.18, linewidth=0.5)
    axis.text(
        0.03,
        0.97,
        (
            f"n = {len(x):,}\n"
            f"Pearson r = {stats['pearson_r']:.3f}\n"
            f"Spearman ρ = {stats['spearman_rho']:.3f}\n"
            f"R² = {stats['linear_regression_r_squared']:.3f}"
        ),
        transform=axis.transAxes,
        va="top",
        ha="left",
        fontsize=8,
        bbox={
            "facecolor": "white",
            "edgecolor": "#B8B8B8",
            "alpha": 0.9,
        },
    )


def draw_scatter_panel(
    axis: plt.Axes,
    gene: str,
    x: np.ndarray,
    y: np.ndarray,
    stats: dict[str, float],
    run_ids: list[str],
    maximum_points: int,
) -> int:
    """Draw the previous deterministic-subset scatter representation."""
    plot_x, plot_y = deterministic_plot_subset(x, y, maximum_points)
    axis.scatter(
        plot_x,
        plot_y,
        s=2.0,
        alpha=0.16,
        color="#1674B8",
        edgecolors="none",
        rasterized=True,
    )
    decorate_gene_panel(axis, gene, x, y, stats, run_ids)
    return len(plot_x)


def draw_density_panel(
    axis: plt.Axes,
    gene: str,
    x: np.ndarray,
    y: np.ndarray,
    stats: dict[str, float],
    run_ids: list[str],
    gridsize: int,
) -> int:
    """Draw every paired target base in logarithmically coloured hexagons."""
    if gridsize < 1:
        raise ValueError("density-gridsize must be positive.")
    density = axis.hexbin(
        x,
        y,
        gridsize=gridsize,
        bins="log",
        mincnt=1,
        cmap="viridis",
        linewidths=0.0,
        rasterized=True,
    )
    colour_bar = axis.get_figure().colorbar(density, ax=axis, pad=0.02)
    colour_bar.set_label("target bases per bin (log scale)", fontsize=8)
    colour_bar.ax.tick_params(labelsize=7)
    decorate_gene_panel(axis, gene, x, y, stats, run_ids)
    return len(x)


def save_selected_figure(
    path: Path,
    title: str,
    genes: list[str],
    profiles: dict[str, tuple[np.ndarray, np.ndarray, dict[str, float]]],
    run_ids: list[str],
    maximum_points: int,
    plot_style: str,
    density_gridsize: int,
) -> None:
    if plot_style not in {"density", "scatter"}:
        raise ValueError(f"Unsupported plot style: {plot_style}")
    figure, axes = plt.subplots(1, len(genes), figsize=(5.1 * len(genes), 4.8))
    for axis, gene in zip(np.atleast_1d(axes), genes):
        x, y, stats = profiles[gene]
        if plot_style == "density":
            draw_density_panel(
                axis,
                gene,
                x,
                y,
                stats,
                run_ids,
                density_gridsize,
            )
        else:
            draw_scatter_panel(
                axis,
                gene,
                x,
                y,
                stats,
                run_ids,
                maximum_points,
            )
    figure.suptitle(title, fontsize=14, fontweight="bold", color="#0B6FB8")
    figure.tight_layout(rect=(0, 0, 1, 0.93))
    figure.savefig(path, dpi=300, bbox_inches="tight")
    plt.close(figure)


def format_statistic(value: float) -> str:
    return f"{value:.6f}" if math.isfinite(value) else "NA"


def main() -> int:
    args = parse_arguments()
    resolved_samtools = shutil.which(args.samtools)
    if resolved_samtools is None:
        raise ValueError(f"samtools executable not found: {args.samtools}")

    sample_rows = read_tsv(args.samples)
    gene_rows = read_tsv(args.gene_coverage)
    comparison_rows = read_tsv(args.gene_comparison)
    sample_runs = method_sample_runs(sample_rows)
    bams = bam_paths(sample_rows, sample_runs)
    targets = target_gene_mapping(gene_rows, sample_runs)
    exclude_flags = panel_exclude_flags(args.filtering_config)
    run_ids, first_samples, second_samples = experiment_sample_order(sample_runs)
    ordered_samples = first_samples + second_samples
    ordered_paths = [bams[sample] for sample in ordered_samples]
    first_count = len(first_samples)

    profiles: dict[
        str, tuple[np.ndarray, np.ndarray, dict[str, float]]
    ] = {}
    statistic_rows: list[dict[str, object]] = []
    for target_key, target in sorted(targets.items()):
        depth_matrix = run_samtools_depth(
            resolved_samtools,
            exclude_flags,
            str(target["chromosome"]),
            int(target["start"]),
            int(target["end"]),
            ordered_paths,
        )
        first_mean = np.mean(depth_matrix[:, :first_count], axis=1)
        second_mean = np.mean(depth_matrix[:, first_count:], axis=1)
        stats = regression_statistics(first_mean, second_mean)
        genes = target["genes"]
        assert isinstance(genes, set)
        shared_target = len(genes) > 1
        for gene in sorted(genes):
            profiles[gene] = (first_mean, second_mean, stats)
            statistic_rows.append(
                {
                    "gene_symbol": gene,
                    "target_key": target_key,
                    "shared_target_interval": str(shared_target).lower(),
                    "genes_sharing_target_interval": ";".join(sorted(genes)),
                    "run_1_id": run_ids[0],
                    "run_1_library_count": len(first_samples),
                    "run_2_id": run_ids[1],
                    "run_2_library_count": len(second_samples),
                    "target_base_count": len(first_mean),
                    **{
                        key: format_statistic(value)
                        for key, value in stats.items()
                    },
                }
            )
        print(f"[INFO] Base-wise reproducibility calculated: {target_key}")

    comparison_genes = {row["gene_symbol"] for row in comparison_rows}
    if comparison_genes != set(profiles):
        raise ValueError(
            "Gene sets differ between relative-coverage and base-wise inputs."
        )

    args.output_dir.mkdir(parents=True, exist_ok=True)
    statistics_path = args.output_dir / "basewise_reproducibility_by_gene.tsv"
    write_tsv(
        statistics_path,
        sorted(statistic_rows, key=lambda row: str(row["gene_symbol"])),
        [
            "gene_symbol",
            "target_key",
            "shared_target_interval",
            "genes_sharing_target_interval",
            "run_1_id",
            "run_1_library_count",
            "run_2_id",
            "run_2_library_count",
            "target_base_count",
            "pearson_r",
            "spearman_rho",
            "linear_regression_r_squared",
            "linear_regression_slope",
            "linear_regression_intercept",
        ],
    )

    high_genes, low_genes = selected_genes(comparison_rows, args.top_n)
    selection_rows = [
        {"panel": panel, "selection_rank": rank, "gene_symbol": gene}
        for panel, genes in (
            ("high_relative_coverage", high_genes),
            ("low_relative_coverage", low_genes),
        )
        for rank, gene in enumerate(genes, start=1)
    ]
    write_tsv(
        args.output_dir / "basewise_figure_gene_selection.tsv",
        selection_rows,
        ["panel", "selection_rank", "gene_symbol"],
    )
    save_selected_figure(
        args.output_dir / "basewise_reproducibility_high_coverage_genes.png",
        "Base-wise coverage reproducibility: highest relative coverage",
        high_genes,
        profiles,
        run_ids,
        args.maximum_plotted_points,
        args.plot_style,
        args.density_gridsize,
    )
    save_selected_figure(
        args.output_dir / "basewise_reproducibility_low_coverage_genes.png",
        "Base-wise coverage reproducibility: lowest relative coverage",
        low_genes,
        profiles,
        run_ids,
        args.maximum_plotted_points,
        args.plot_style,
        args.density_gridsize,
    )

    print(f"[INFO] Statistics output: {statistics_path}")
    print(f"[INFO] High-coverage figure genes: {', '.join(high_genes)}")
    print(f"[INFO] Low-coverage figure genes: {', '.join(low_genes)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
