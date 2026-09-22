#!/usr/bin/env python3

"""
Create a technical small-variant summary for one clinical sample.

The script reads a VCF or compressed VCF using bcftools and reports:
- total records
- PASS records
- SNVs
- insertions
- deletions
- complex or other variants
- heterozygous genotypes
- homozygous-alternate genotypes

The output is a tab-separated file with one row.
"""

from __future__ import annotations

import argparse
import csv
import subprocess
import sys
from pathlib import Path


def parse_arguments() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Summarize small variants for one clinical sample."
    )
    parser.add_argument(
        "--sample-id",
        required=True,
        help="Sample identifier used in the output.",
    )
    parser.add_argument(
        "--input-vcf",
        required=True,
        type=Path,
        help="Input VCF or VCF.GZ file.",
    )
    parser.add_argument(
        "--output-tsv",
        required=True,
        type=Path,
        help="Output TSV file.",
    )
    return parser.parse_args()


def run_bcftools_query(vcf_path: Path) -> list[str]:
    """
    Return VCF records in a compact tab-separated representation.

    Fields:
        FILTER, REF, ALT, GT
    """
    command = [
        "bcftools",
        "query",
        "-f",
        "%FILTER\\t%REF\\t%ALT[\\t%GT]\\n",
        str(vcf_path),
    ]

    try:
        result = subprocess.run(
            command,
            check=True,
            capture_output=True,
            text=True,
        )
    except FileNotFoundError:
        raise RuntimeError(
            "bcftools was not found in the active environment."
        )
    except subprocess.CalledProcessError as error:
        raise RuntimeError(
            "bcftools could not read the input VCF:\n"
            + error.stderr.strip()
        )

    return [
        line
        for line in result.stdout.splitlines()
        if line.strip()
    ]


def classify_variant(ref: str, alt: str) -> str:
    """
    Classify one biallelic REF/ALT pair.

    Multiallelic ALT values are classified as complex_or_other because
    one VCF record may contain more than one variant type.
    """
    if "," in alt or alt.startswith("<") or "[" in alt or "]" in alt:
        return "complex_or_other"

    if len(ref) == 1 and len(alt) == 1:
        return "snv"

    if len(alt) > len(ref):
        return "insertion"

    if len(ref) > len(alt):
        return "deletion"

    return "complex_or_other"


def normalize_genotype(genotype: str) -> str:
    return genotype.replace("|", "/")


def summarize_records(records: list[str]) -> dict[str, int]:
    summary = {
        "total_variants": 0,
        "pass_variants": 0,
        "snvs": 0,
        "insertions": 0,
        "deletions": 0,
        "complex_or_other": 0,
        "heterozygous": 0,
        "homozygous_alternate": 0,
        "missing_or_other_genotype": 0,
    }

    for line in records:
        fields = line.split("\t")

        if len(fields) < 4:
            raise RuntimeError(
                f"Unexpected bcftools query output: {line!r}"
            )

        variant_filter, ref, alt, genotype = fields[:4]

        summary["total_variants"] += 1

        if variant_filter in {"PASS", "."}:
            summary["pass_variants"] += 1

        variant_type = classify_variant(ref, alt)
        summary[
            {
                "snv": "snvs",
                "insertion": "insertions",
                "deletion": "deletions",
                "complex_or_other": "complex_or_other",
            }[variant_type]
        ] += 1

        gt = normalize_genotype(genotype)

        if gt in {"0/1", "1/0"}:
            summary["heterozygous"] += 1
        elif gt == "1/1":
            summary["homozygous_alternate"] += 1
        else:
            summary["missing_or_other_genotype"] += 1

    return summary


def write_summary(
    sample_id: str,
    summary: dict[str, int],
    output_path: Path,
) -> None:
    output_path.parent.mkdir(parents=True, exist_ok=True)

    fieldnames = [
        "sample_id",
        "total_variants",
        "pass_variants",
        "snvs",
        "insertions",
        "deletions",
        "complex_or_other",
        "heterozygous",
        "homozygous_alternate",
        "missing_or_other_genotype",
    ]

    row = {"sample_id": sample_id, **summary}

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
        writer.writerow(row)


def main() -> int:
    args = parse_arguments()

    input_vcf = args.input_vcf.expanduser().resolve()
    output_tsv = args.output_tsv.expanduser().resolve()

    if not input_vcf.is_file():
        print(
            f"ERROR: Input VCF does not exist: {input_vcf}",
            file=sys.stderr,
        )
        return 1

    try:
        records = run_bcftools_query(input_vcf)
        summary = summarize_records(records)
        write_summary(args.sample_id, summary, output_tsv)
    except RuntimeError as error:
        print(f"ERROR: {error}", file=sys.stderr)
        return 1

    print(f"Clinical small-variant summary created: {output_tsv}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
