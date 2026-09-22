from __future__ import annotations

import importlib.util
from pathlib import Path

import pytest


PROJECT_ROOT = Path(__file__).resolve().parents[1]
MODULE_PATH = (
    PROJECT_ROOT
    / "scripts"
    / "python"
    / "summarize_method_development.py"
)
SPEC = importlib.util.spec_from_file_location(
    "summarize_method_development",
    MODULE_PATH,
)
assert SPEC is not None
assert SPEC.loader is not None
summary = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(summary)


def make_sample(sample_id: str, run_id: str) -> dict[str, str]:
    return {
        "sample_id": sample_id,
        "analysis_role": "method_development",
        "run_id": run_id,
        "include": "true",
    }


def make_qc_row(
    sample_id: str,
    offset: float,
) -> dict[str, str]:
    values = {
        metric: str(index + offset)
        for index, metric in enumerate(summary.METRICS, start=1)
    }
    values["sample_id"] = sample_id
    values["panel_qc_pass"] = "false"
    return values


def test_parse_group_uses_final_numeric_suffix() -> None:
    assert summary.parse_group("E1_Sample20_4") == (
        "E1_Sample20",
        4,
    )


def test_create_library_and_group_rows() -> None:
    samples = [
        make_sample("E1_Sample20_1", "experiment1"),
        make_sample("E1_Sample20_2", "experiment1"),
    ]
    rows = [
        make_qc_row("E1_Sample20_1", 0.0),
        make_qc_row("E1_Sample20_2", 2.0),
    ]

    library_rows = summary.create_library_rows(
        samples=samples,
        alignment_rows=rows,
        target_rows=rows,
        coverage_rows=rows,
    )
    group_rows = summary.create_group_rows(library_rows)

    assert len(library_rows) == 2
    assert len(group_rows) == 2
    assert group_rows[0]["summary_level"] == "experiment"
    assert group_rows[0]["library_count"] == 2
    assert float(
        group_rows[0]["mean_primary_mapped_reads"]
    ) == pytest.approx(2.0)
    assert group_rows[1]["technical_group"] == "E1_Sample20"


def test_missing_shared_qc_sample_is_rejected() -> None:
    sample = make_sample("E1_Sample20_1", "experiment1")
    row = make_qc_row("different_sample", 0.0)

    with pytest.raises(ValueError, match="Missing shared-QC row"):
        summary.create_library_rows(
            samples=[sample],
            alignment_rows=[row],
            target_rows=[row],
            coverage_rows=[row],
        )
