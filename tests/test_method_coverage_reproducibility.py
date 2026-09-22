from __future__ import annotations

import importlib.util
from pathlib import Path

import numpy as np
import pytest


PROJECT_ROOT = Path(__file__).resolve().parents[1]


def load_module(name: str, relative_path: str):
    path = PROJECT_ROOT / relative_path
    spec = importlib.util.spec_from_file_location(name, path)
    assert spec is not None
    assert spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


gene_summary = load_module(
    "method_gene_reproducibility",
    "scripts/python/summarize_method_coverage_reproducibility.py",
)
basewise = load_module(
    "method_basewise_reproducibility",
    "scripts/python/plot_basewise_coverage_reproducibility.py",
)


def sample(sample_id: str, run_id: str) -> dict[str, str]:
    return {
        "sample_id": sample_id,
        "run_id": run_id,
        "analysis_role": "method_development",
        "include": "true",
    }


def gene_row(
    sample_id: str,
    gene: str,
    target: str,
    mean_depth: float,
) -> dict[str, str]:
    return {
        "sample_id": sample_id,
        "gene_symbol": gene,
        "target_key": target,
        "mean_depth": str(mean_depth),
    }


def test_gene_coverage_is_normalised_within_each_library() -> None:
    samples = [
        sample("E1_A_1", "experiment1"),
        sample("E1_A_2", "experiment1"),
        sample("E2_A_1", "experiment2"),
        sample("E2_A_2", "experiment2"),
    ]
    panels = [
        {"sample_id": "E1_A_1", "mean_depth": "10"},
        {"sample_id": "E1_A_2", "mean_depth": "20"},
        {"sample_id": "E2_A_1", "mean_depth": "40"},
        {"sample_id": "E2_A_2", "mean_depth": "80"},
    ]
    genes = []
    for sample_id, panel_depth in (
        ("E1_A_1", 10),
        ("E1_A_2", 20),
        ("E2_A_1", 40),
        ("E2_A_2", 80),
    ):
        genes.extend(
            [
                gene_row(sample_id, "HIGH", "chr1:0-10", panel_depth * 2),
                gene_row(sample_id, "LOW", "chr2:0-10", panel_depth * 0.5),
            ]
        )

    library, summary, comparison, correlation = gene_summary.create_outputs(
        samples, genes, panels
    )

    assert len(library) == 8
    assert len(summary) == 4
    assert [row["gene_symbol"] for row in comparison] == ["HIGH", "LOW"]
    assert float(comparison[0]["mean_relative_coverage_across_runs"]) == 2.0
    assert float(comparison[1]["mean_relative_coverage_across_runs"]) == 0.5
    assert float(correlation[0]["pearson_r"]) == pytest.approx(1.0)
    assert float(correlation[0]["spearman_rho"]) == pytest.approx(1.0)


def test_heatmap_uses_one_common_matrix_and_writes_values(tmp_path: Path) -> None:
    libraries = [
        {
            "sample_id": "E1_A_1",
            "run_id": "experiment1",
            "gene_symbol": "HIGH",
            "relative_gene_coverage": "2.0",
        },
        {
            "sample_id": "E2_A_1",
            "run_id": "experiment2",
            "gene_symbol": "HIGH",
            "relative_gene_coverage": "1.5",
        },
        {
            "sample_id": "E1_A_1",
            "run_id": "experiment1",
            "gene_symbol": "LOW",
            "relative_gene_coverage": "0.5",
        },
        {
            "sample_id": "E2_A_1",
            "run_id": "experiment2",
            "gene_symbol": "LOW",
            "relative_gene_coverage": "0.25",
        },
    ]
    comparison = [
        {"gene_symbol": "HIGH"},
        {"gene_symbol": "LOW"},
    ]

    gene_summary.create_heatmap(libraries, comparison, tmp_path)

    assert (tmp_path / "gene_relative_coverage_heatmap.png").is_file()
    values = gene_summary.read_tsv(
        tmp_path / "gene_relative_coverage_heatmap_values.tsv"
    )
    assert len(values) == 4
    assert float(values[0]["log2_relative_gene_coverage"]) == pytest.approx(1.0)


def test_regression_statistics_use_all_paired_values() -> None:
    x = np.asarray([0.0, 1.0, 2.0, 3.0])
    y = np.asarray([1.0, 3.0, 5.0, 7.0])
    stats = basewise.regression_statistics(x, y)
    assert stats["pearson_r"] == pytest.approx(1.0)
    assert stats["spearman_rho"] == pytest.approx(1.0)
    assert stats["linear_regression_r_squared"] == pytest.approx(1.0)
    assert stats["linear_regression_slope"] == pytest.approx(2.0)
    assert stats["linear_regression_intercept"] == pytest.approx(1.0)


def test_selected_genes_are_determined_by_across_run_mean() -> None:
    rows = [
        {
            "gene_symbol": f"G{index}",
            "mean_relative_coverage_across_runs": str(value),
        }
        for index, value in enumerate([1.1, 0.2, 2.0, 0.5, 1.5, 0.1])
    ]
    high, low = basewise.selected_genes(rows, top_n=2)
    assert high == ["G2", "G4"]
    assert low == ["G5", "G1"]


def test_plot_subsampling_is_deterministic() -> None:
    x = np.arange(100, dtype=float)
    y = x * 2
    first_x, first_y = basewise.deterministic_plot_subset(x, y, 10)
    second_x, second_y = basewise.deterministic_plot_subset(x, y, 10)
    assert np.array_equal(first_x, second_x)
    assert np.array_equal(first_y, second_y)
    assert len(first_x) == 10


def test_density_panel_draws_all_target_bases_with_logarithmic_colour() -> None:
    x = np.arange(1000, dtype=float)
    y = x * 1.2
    stats = basewise.regression_statistics(x, y)
    figure, axis = basewise.plt.subplots()

    drawn = basewise.draw_density_panel(
        axis,
        "GENE",
        x,
        y,
        stats,
        ["experiment1", "experiment2"],
        gridsize=40,
    )

    density = axis.collections[0]
    assert drawn == len(x)
    assert np.sum(density.get_array()) == pytest.approx(len(x))
    assert density.norm.__class__.__name__ == "LogNorm"
    basewise.plt.close(figure)


def test_filtering_config_supplies_panel_flag_mask(tmp_path: Path) -> None:
    config = tmp_path / "filtering.yaml"
    config.write_text(
        "filters:\n  panel_analysis:\n    exclude_flags: 3844\n",
        encoding="utf-8",
    )
    assert basewise.panel_exclude_flags(config) == 3844
