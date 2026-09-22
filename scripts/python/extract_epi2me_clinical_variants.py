#!/usr/bin/env python3

"""
Extract an analysis-ready variant table from an EPI2ME Human Variation
ClinVar-annotated small-variant VCF.

The script creates:
1. A complete table containing all VCF records.
2. A review table containing records with clinically relevant or uncertain
   ClinVar classifications.

This script reports existing EPI2ME and ClinVar annotations. It does not
perform clinical classification or diagnostic interpretation.
"""

from __future__ import annotations

import argparse
import csv
import subprocess
import sys
from pathlib import Path
from typing import Dict, List, Tuple


ANN_FIELD_NAMES = [
    "allele",
    "consequence",
    "impact",
    "gene",
    "gene_id",
    "feature_type",
    "feature_id",
    "transcript_biotype",
    "rank",
    "hgvs_c",
    "hgvs_p",
    "cdna_position",
    "cds_position",
    "protein_position",
    "distance",
    "messages",
]

IMPACT_RANK = {
    "HIGH": 4,
    "MODERATE": 3,
    "LOW": 2,
    "MODIFIER": 1,
    "": 0,
    ".": 0,
}

REVIEW_TERMS = (
    "pathogenic",
    "likely_pathogenic",
    "uncertain_significance",
    "conflicting",
    "risk_factor",
    "association",
    "drug_response",
    "affects",
)


def parse_arguments() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Extract detailed annotated variants from an EPI2ME "
            "ClinVar VCF."
        )
    )
    parser.add_argument(
        "--sample-id",
        required=True,
        help="Sample identifier written to the output tables.",
    )
    parser.add_argument(
        "--input-vcf",
        required=True,
        type=Path,
        help="EPI2ME *.wf_snp_clinvar.vcf.gz input file.",
    )
    parser.add_argument(
        "--output-all",
        required=True,
        type=Path,
        help="Output TSV containing all annotated variants.",
    )
    parser.add_argument(
        "--output-review",
        required=True,
        type=Path,
        help="Output TSV containing variants selected for review.",
    )
    return parser.parse_args()


def run_bcftools_query(vcf_path: Path) -> List[str]:
    """
    Extract one tab-separated line per VCF record.

    Missing INFO values are returned as dots by bcftools.
    """
    query_format = (
        "%CHROM\\t"
        "%POS\\t"
        "%ID\\t"
        "%REF\\t"
        "%ALT\\t"
        "%FILTER\\t"
        "%QUAL\\t"
        "%INFO/ANN\\t"
        "%INFO/CLNSIG\\t"
        "%INFO/CLNDN\\t"
        "%INFO/CLNREVSTAT\\t"
        "%INFO/CLNVI\\t"
        "%INFO/RS\\t"
        "[%GT\\t%GQ\\t%DP\\t%AD\\t%AF]\\n"
    )

    command = [
        "bcftools",
        "query",
        "-f",
        query_format,
        str(vcf_path),
    ]

    try:
        result = subprocess.run(
            command,
            check=True,
            capture_output=True,
            text=True,
        )
    except FileNotFoundError as error:
        raise RuntimeError(
            "bcftools was not found in the active environment."
        ) from error
    except subprocess.CalledProcessError as error:
        raise RuntimeError(
            "bcftools could not query the input VCF:\n"
            + error.stderr.strip()
        ) from error

    return [
        line
        for line in result.stdout.splitlines()
        if line.strip()
    ]


def clean_value(value: str) -> str:
    """Convert VCF missing values to an empty string."""
    value = value.strip()

    if value in {"", ".", "NA", "N/A"}:
        return ""

    return value


def parse_ann_entry(entry: str) -> Dict[str, str]:
    """Parse one SnpEff-style ANN annotation entry."""
    values = entry.split("|")

    if len(values) < len(ANN_FIELD_NAMES):
        values.extend(
            [""] * (len(ANN_FIELD_NAMES) - len(values))
        )

    return {
        field_name: clean_value(value)
        for field_name, value in zip(
            ANN_FIELD_NAMES,
            values,
        )
    }


def annotation_score(annotation: Dict[str, str]) -> Tuple[int, int, int]:
    """
    Rank annotations for selecting one representative transcript.

    Priority:
    1. Higher predicted impact.
    2. Protein-coding transcript.
    3. Annotation with an HGVS protein description.
    """
    impact_score = IMPACT_RANK.get(
        annotation.get("impact", "").upper(),
        0,
    )

    protein_coding_score = int(
        annotation.get("transcript_biotype") == "protein_coding"
    )

    hgvs_p_score = int(bool(annotation.get("hgvs_p")))

    return (
        impact_score,
        protein_coding_score,
        hgvs_p_score,
    )


def select_representative_annotation(
    ann_value: str,
) -> Tuple[Dict[str, str], str, str]:
    """
    Select one representative annotation and retain all unique genes and
    consequences found in the ANN field.
    """
    ann_value = clean_value(ann_value)

    if not ann_value:
        empty = {
            field_name: ""
            for field_name in ANN_FIELD_NAMES
        }
        return empty, "", ""

    annotations = [
        parse_ann_entry(entry)
        for entry in ann_value.split(",")
        if entry.strip()
    ]

    if not annotations:
        empty = {
            field_name: ""
            for field_name in ANN_FIELD_NAMES
        }
        return empty, "", ""

    representative = max(
        annotations,
        key=annotation_score,
    )

    all_genes = sorted(
        {
            annotation["gene"]
            for annotation in annotations
            if annotation.get("gene")
        }
    )

    all_consequences = sorted(
        {
            annotation["consequence"]
            for annotation in annotations
            if annotation.get("consequence")
        }
    )

    return (
        representative,
        ",".join(all_genes),
        ",".join(all_consequences),
    )


def classify_variant(ref: str, alt: str) -> str:
    """Classify the VCF record by REF and ALT allele lengths."""
    if (
        "," in alt
        or alt.startswith("<")
        or "[" in alt
        or "]" in alt
    ):
        return "complex_or_multiallelic"

    if len(ref) == 1 and len(alt) == 1:
        return "SNV"

    if len(alt) > len(ref):
        return "insertion"

    if len(ref) > len(alt):
        return "deletion"

    return "complex_or_other"


def normalize_clinvar_text(value: str) -> str:
    """Normalize ClinVar text only for review-table matching."""
    return (
        clean_value(value)
        .lower()
        .replace(" ", "_")
        .replace("/", "_")
    )


def requires_review(clinvar_significance: str) -> bool:
    """
    Return True when the existing ClinVar classification merits review.

    Benign and likely benign records alone are not selected.
    Records without a ClinVar classification are not selected.
    """
    normalized = normalize_clinvar_text(
        clinvar_significance
    )

    if not normalized:
        return False

    return any(
        term in normalized
        for term in REVIEW_TERMS
    )


def parse_record(
    sample_id: str,
    line: str,
) -> Dict[str, str]:
    fields = line.split("\t")

    if len(fields) < 18:
        raise RuntimeError(
            "Unexpected bcftools output with "
            f"{len(fields)} columns:\n{line}"
        )

    (
        chromosome,
        position,
        record_id,
        ref,
        alt,
        variant_filter,
        quality,
        ann,
        clinvar_significance,
        clinvar_condition,
        clinvar_review_status,
        clinvar_variant_identifiers,
        dbsnp_rs,
        genotype,
        genotype_quality,
        depth,
        allele_depth,
        allele_fraction,
    ) = fields[:18]

    (
        representative,
        all_genes,
        all_consequences,
    ) = select_representative_annotation(ann)

    return {
        "sample_id": sample_id,
        "chromosome": clean_value(chromosome),
        "position_1based": clean_value(position),
        "record_id": clean_value(record_id),
        "ref": clean_value(ref),
        "alt": clean_value(alt),
        "variant_type": classify_variant(ref, alt),
        "filter": clean_value(variant_filter),
        "quality": clean_value(quality),
        "genotype": clean_value(genotype),
        "genotype_quality": clean_value(genotype_quality),
        "depth": clean_value(depth),
        "allele_depth": clean_value(allele_depth),
        "allele_fraction": clean_value(allele_fraction),
        "gene": representative.get("gene", ""),
        "all_annotated_genes": all_genes,
        "gene_id": representative.get("gene_id", ""),
        "transcript_id": representative.get("feature_id", ""),
        "transcript_biotype": representative.get(
            "transcript_biotype",
            "",
        ),
        "consequence": representative.get(
            "consequence",
            "",
        ),
        "all_consequences": all_consequences,
        "annotation_impact": representative.get(
            "impact",
            "",
        ),
        "hgvs_c": representative.get("hgvs_c", ""),
        "hgvs_p": representative.get("hgvs_p", ""),
        "clinvar_significance": clean_value(
            clinvar_significance
        ),
        "clinvar_condition": clean_value(
            clinvar_condition
        ),
        "clinvar_review_status": clean_value(
            clinvar_review_status
        ),
        "clinvar_variant_identifiers": clean_value(
            clinvar_variant_identifiers
        ),
        "dbsnp_rs": clean_value(dbsnp_rs),
    }


def write_tsv(
    rows: List[Dict[str, str]],
    output_path: Path,
    fieldnames: List[str],
) -> None:
    output_path.parent.mkdir(
        parents=True,
        exist_ok=True,
    )

    with output_path.open(
        "w",
        encoding="utf-8",
        newline="",
    ) as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=fieldnames,
            delimiter="\t",
            extrasaction="ignore",
        )
        writer.writeheader()
        writer.writerows(rows)


def main() -> int:
    args = parse_arguments()

    input_vcf = args.input_vcf.expanduser().resolve()
    output_all = args.output_all.expanduser().resolve()
    output_review = args.output_review.expanduser().resolve()

    if not input_vcf.is_file():
        print(
            f"ERROR: Input VCF does not exist: {input_vcf}",
            file=sys.stderr,
        )
        return 1

    fieldnames = [
        "sample_id",
        "chromosome",
        "position_1based",
        "record_id",
        "ref",
        "alt",
        "variant_type",
        "filter",
        "quality",
        "genotype",
        "genotype_quality",
        "depth",
        "allele_depth",
        "allele_fraction",
        "gene",
        "all_annotated_genes",
        "gene_id",
        "transcript_id",
        "transcript_biotype",
        "consequence",
        "all_consequences",
        "annotation_impact",
        "hgvs_c",
        "hgvs_p",
        "clinvar_significance",
        "clinvar_condition",
        "clinvar_review_status",
        "clinvar_variant_identifiers",
        "dbsnp_rs",
    ]

    try:
        query_lines = run_bcftools_query(input_vcf)

        all_rows = [
            parse_record(args.sample_id, line)
            for line in query_lines
        ]

        review_rows = [
            row
            for row in all_rows
            if requires_review(
                row["clinvar_significance"]
            )
        ]

        write_tsv(
            all_rows,
            output_all,
            fieldnames,
        )
        write_tsv(
            review_rows,
            output_review,
            fieldnames,
        )

    except RuntimeError as error:
        print(
            f"ERROR: {error}",
            file=sys.stderr,
        )
        return 1

    print(f"All annotated variants: {len(all_rows)}")
    print(f"Variants selected for review: {len(review_rows)}")
    print(f"Complete table: {output_all}")
    print(f"Review table:   {output_review}")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
