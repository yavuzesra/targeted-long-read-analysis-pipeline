from __future__ import annotations

import argparse
import csv
import importlib.util
from pathlib import Path


PROJECT_ROOT = Path(__file__).resolve().parents[1]
MODULE_PATH = (
    PROJECT_ROOT
    / "scripts"
    / "optional"
    / "summarize_sv_phasing.py"
)

SPEC = importlib.util.spec_from_file_location(
    "summarize_sv_phasing",
    MODULE_PATH,
)

assert SPEC is not None
assert SPEC.loader is not None

sv_phasing = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(sv_phasing)


def read_output_row(path: Path) -> dict[str, str]:
    with path.open(
        "r",
        encoding="utf-8",
        newline="",
    ) as handle:
        reader = csv.DictReader(
            handle,
            delimiter="\t",
        )
        return next(reader)


def test_default_row_contains_na_values() -> None:
    row = sv_phasing.default_row("CLINICAL_TEST")

    assert row["sample_id"] == "CLINICAL_TEST"

    for column in sv_phasing.OUTPUT_COLUMNS:
        if column != "sample_id":
            assert row[column] == "NA"


def test_sv_length_parsing() -> None:
    assert sv_phasing.parse_sv_length("-447") == [447]
    assert sv_phasing.parse_sv_length("100,-250") == [100, 250]
    assert sv_phasing.parse_sv_length(".") == []
    assert sv_phasing.parse_sv_length("invalid") == []


def test_genotype_classification() -> None:
    assert sv_phasing.is_heterozygous("0/1")
    assert sv_phasing.is_heterozygous("1|0")
    assert not sv_phasing.is_heterozygous("1/1")

    assert sv_phasing.is_homozygous_alternative("1/1")
    assert sv_phasing.is_homozygous_alternative("1|1")
    assert not sv_phasing.is_homozygous_alternative("0/1")


def test_all_missing_workflow_outputs_produce_na(
    tmp_path: Path,
    monkeypatch,
) -> None:
    sample_sheet = tmp_path / "samples.tsv"
    reference = tmp_path / "reference.fna"
    output = tmp_path / "summary.tsv"
    missing_directory = tmp_path / "missing_output"

    sample_sheet.write_text(
        "sample_id\toutput_directory\n"
        f"CLINICAL_TEST\t{missing_directory}\n",
        encoding="utf-8",
    )

    reference.write_text(
        ">chr1\nACGT\n",
        encoding="utf-8",
    )

    monkeypatch.setattr(
        sv_phasing,
        "parse_arguments",
        lambda: argparse.Namespace(
            sample_sheet=sample_sheet,
            reference_fasta=reference,
            output_tsv=output,
        ),
    )

    monkeypatch.setattr(
        sv_phasing,
        "check_tool",
        lambda _tool: None,
    )

    exit_code = sv_phasing.main()

    assert exit_code == 0
    assert output.is_file()

    row = read_output_row(output)

    assert row["sample_id"] == "CLINICAL_TEST"

    for column in sv_phasing.OUTPUT_COLUMNS:
        if column != "sample_id":
            assert row[column] == "NA"


def test_available_inputs_are_processed_independently(
    tmp_path: Path,
    monkeypatch,
) -> None:
    sample_sheet = tmp_path / "samples.tsv"
    reference = tmp_path / "reference.fna"
    output = tmp_path / "summary.tsv"
    workflow_directory = tmp_path / "workflow_output"

    workflow_directory.mkdir()

    sample_sheet.write_text(
        "sample_id\toutput_directory\n"
        f"CLINICAL_TEST\t{workflow_directory}\n",
        encoding="utf-8",
    )

    reference.write_text(
        ">chr1\nACGT\n",
        encoding="utf-8",
    )

    sv_vcf = workflow_directory / "CLINICAL_TEST.wf_sv.vcf.gz"
    snp_vcf = workflow_directory / "CLINICAL_TEST.wf_snp.vcf.gz"

    sv_vcf.touch()
    snp_vcf.touch()

    monkeypatch.setattr(
        sv_phasing,
        "parse_arguments",
        lambda: argparse.Namespace(
            sample_sheet=sample_sheet,
            reference_fasta=reference,
            output_tsv=output,
        ),
    )

    monkeypatch.setattr(
        sv_phasing,
        "check_tool",
        lambda _tool: None,
    )

    def fake_sv_summary(
        _path: Path,
        row: dict[str, str | int],
    ) -> None:
        row["pass_sv"] = 31
        row["ins"] = 16
        row["del"] = 15
        row["inv"] = 0

    def fake_snp_summary(
        _path: Path,
        row: dict[str, str | int],
    ) -> None:
        row["het_small_variants"] = 1596
        row["phased_het_small_variants"] = "NA"
        row["phased_het_pct"] = "NA"

    monkeypatch.setattr(
        sv_phasing,
        "summarize_sv",
        fake_sv_summary,
    )

    monkeypatch.setattr(
        sv_phasing,
        "summarize_snp_phasing",
        fake_snp_summary,
    )

    exit_code = sv_phasing.main()

    assert exit_code == 0

    row = read_output_row(output)

    assert row["pass_sv"] == "31"
    assert row["ins"] == "16"
    assert row["del"] == "15"
    assert row["inv"] == "0"

    assert row["het_small_variants"] == "1596"

    # No haplotagged CRAM was created.
    assert row["hp1_reads"] == "NA"
    assert row["hp2_reads"] == "NA"
    assert row["haplotagged_reads"] == "NA"
    assert row["total_reads"] == "NA"
    assert row["haplotagged_pct"] == "NA"
