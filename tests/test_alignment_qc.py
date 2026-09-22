from __future__ import annotations

import importlib.util
from pathlib import Path


PROJECT_ROOT = Path(__file__).resolve().parents[1]
MODULE_PATH = (
    PROJECT_ROOT / "scripts/python/alignment_qc.py"
)

SPEC = importlib.util.spec_from_file_location(
    "alignment_qc",
    MODULE_PATH,
)

if SPEC is None or SPEC.loader is None:
    raise RuntimeError(
        f"Unable to load alignment_qc module: {MODULE_PATH}"
    )

alignment_qc = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(alignment_qc)


def test_n50_standard_lengths() -> None:
    lengths = [100, 200, 300, 400]

    result = alignment_qc.calculate_n50(lengths)

    assert result == 300


def test_n50_single_read() -> None:
    result = alignment_qc.calculate_n50([1250])

    assert result == 1250


def test_n50_empty_input() -> None:
    result = alignment_qc.calculate_n50([])

    assert result is None


def test_format_float_uses_four_decimals() -> None:
    result = alignment_qc.format_value(97.72749)

    assert result == "97.7275"


def test_format_none_returns_na() -> None:
    result = alignment_qc.format_value(None)

    assert result == "NA"


def test_alignment_accuracy_from_error_rate() -> None:
    error_rate = 0.022725

    accuracy = 100.0 * (1.0 - error_rate)

    assert round(accuracy, 4) == 97.7275


def test_primary_count_consistency() -> None:
    primary_total = 36716
    primary_mapped = 36713
    primary_unmapped = 3

    difference = (
        primary_total
        - primary_mapped
        - primary_unmapped
    )

    assert difference == 0
