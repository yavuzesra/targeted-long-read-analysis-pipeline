#!/usr/bin/env python3
"""
Create an exploratory structural-variant and phasing summary from
EPI2ME wf-human-variation outputs.

Missing files are handled independently:
- missing wf_sv.vcf.gz -> SV metrics are NA
- missing wf_snp.vcf.gz -> small-variant phasing metrics are NA
- missing haplotagged.cram -> read-level haplotagging metrics are NA

Existing but unreadable files are treated as errors.
"""

from __future__ import annotations

import argparse
import csv
import subprocess
import sys
from collections import Counter
from pathlib import Path
from typing import Iterable


OUTPUT_COLUMNS = [
    "sample_id",
    "pass_sv",
    "ins",
    "del",
    "inv",
    "het_sv",
    "hom_alt_sv",
    "phased_sv",
    "min_abs_svlen",
    "max_abs_svlen",
    "het_small_variants",
    "phased_het_small_variants",
    "unphased_het_small_variants",
    "phased_het_pct",
    "hp1_reads",
    "hp2_reads",
    "haplotagged_reads",
    "total_reads",
    "haplotagged_pct",
    "multivariant_phase_blocks",
]


def parse_arguments() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Summarize exploratory structural-variant and phasing outputs."
        )
    )
    parser.add_argument("sample_sheet", type=Path)
    parser.add_argument("reference_fasta", type=Path)
    parser.add_argument("output_tsv", type=Path)
    return parser.parse_args()


def run_command(
    command: list[str],
    *,
    check: bool = True,
) -> subprocess.CompletedProcess[str]:
    return subprocess.run(
        command,
        text=True,
        encoding="utf-8",
        errors="replace",
        capture_output=True,
        check=check,
    )


def check_tool(tool: str) -> None:
    result = run_command([tool, "--version"], check=False)

    if result.returncode != 0:
        raise RuntimeError(f"Required tool is unavailable: {tool}")


def is_heterozygous(genotype: str) -> bool:
    return genotype in {"0/1", "1/0", "0|1", "1|0"}


def is_homozygous_alternative(genotype: str) -> bool:
    return genotype in {"1/1", "1|1"}


def parse_sv_length(value: str) -> list[int]:
    lengths: list[int] = []

    if not value or value == ".":
        return lengths

    for item in value.split(","):
        try:
            lengths.append(abs(int(float(item))))
        except ValueError:
            continue

    return lengths


def default_row(sample_id: str) -> dict[str, str | int]:
    return {column: "NA" for column in OUTPUT_COLUMNS} | {
        "sample_id": sample_id
    }


def summarize_sv(
    path: Path,
    row: dict[str, str | int],
) -> None:
    result = run_command(
        [
            "bcftools",
            "query",
            "-i",
            'FILTER="PASS"',
            "-f",
            r"%INFO/SVTYPE\t%INFO/SVLEN[\t%GT]\n",
            str(path),
        ],
        check=False,
    )

    if result.returncode != 0:
        raise RuntimeError(
            f"SV VCF is present but unreadable: {path}\n{result.stderr}"
        )

    type_counts: Counter[str] = Counter()
    pass_count = 0
    heterozygous = 0
    homozygous_alternative = 0
    phased = 0
    lengths: list[int] = []

    for line in result.stdout.splitlines():
        if not line.strip():
            continue

        fields = line.split("\t")
        sv_type = fields[0] if fields else ""
        sv_length = fields[1] if len(fields) > 1 else "."
        genotype = fields[2] if len(fields) > 2 else "."

        pass_count += 1
        type_counts[sv_type] += 1
        lengths.extend(parse_sv_length(sv_length))

        if is_heterozygous(genotype):
            heterozygous += 1

        if is_homozygous_alternative(genotype):
            homozygous_alternative += 1

        if "|" in genotype:
            phased += 1

    row.update(
        {
            "pass_sv": pass_count,
            "ins": type_counts["INS"],
            "del": type_counts["DEL"],
            "inv": type_counts["INV"],
            "het_sv": heterozygous,
            "hom_alt_sv": homozygous_alternative,
            "phased_sv": phased,
            "min_abs_svlen": min(lengths) if lengths else "NA",
            "max_abs_svlen": max(lengths) if lengths else "NA",
        }
    )


def snp_vcf_has_phase_set(path: Path) -> bool:
    result = run_command(
        ["bcftools", "view", "-h", str(path)],
        check=False,
    )

    if result.returncode != 0:
        raise RuntimeError(
            f"SNP VCF is present but unreadable: {path}\n{result.stderr}"
        )

    return any(
        line.startswith("##FORMAT=<ID=PS,")
        for line in result.stdout.splitlines()
    )


def summarize_snp_phasing(
    path: Path,
    row: dict[str, str | int],
) -> None:
    has_phase_set = snp_vcf_has_phase_set(path)

    if has_phase_set:
        format_string = r"%CHROM\t%POS[\t%GT\t%PS]\n"
    else:
        format_string = r"%CHROM\t%POS[\t%GT]\n"

    result = run_command(
        [
            "bcftools",
            "query",
            "-f",
            format_string,
            str(path),
        ],
        check=False,
    )

    if result.returncode != 0:
        raise RuntimeError(
            f"SNP VCF is present but unreadable: {path}\n{result.stderr}"
        )

    heterozygous = 0
    phased_heterozygous = 0
    unphased_heterozygous = 0
    phase_block_counts: Counter[tuple[str, str]] = Counter()

    for line in result.stdout.splitlines():
        if not line.strip():
            continue

        fields = line.split("\t")

        chromosome = fields[0]
        genotype = fields[2] if len(fields) > 2 else "."
        phase_set = fields[3] if len(fields) > 3 else "."

        if not is_heterozygous(genotype):
            continue

        heterozygous += 1

        if "|" in genotype:
            phased_heterozygous += 1

            if has_phase_set and phase_set not in {"", "."}:
                phase_block_counts[(chromosome, phase_set)] += 1
        else:
            unphased_heterozygous += 1

    row["het_small_variants"] = heterozygous

    if not has_phase_set:
        print(
            f"WARNING: No FORMAT/PS field for {path.name}; "
            "phase-assignment metrics set to NA.",
            file=sys.stderr,
        )
        return

    phased_rate = (
        100.0 * phased_heterozygous / heterozygous
        if heterozygous > 0
        else 0.0
    )

    multi_variant_blocks = sum(
        count >= 2 for count in phase_block_counts.values()
    )

    row.update(
        {
            "phased_het_small_variants": phased_heterozygous,
            "unphased_het_small_variants": unphased_heterozygous,
            "phased_het_pct": f"{phased_rate:.4f}",
            "multivariant_phase_blocks": multi_variant_blocks,
        }
    )


def summarize_haplotagged_cram(
    path: Path,
    reference_fasta: Path,
    row: dict[str, str | int],
) -> None:
    quickcheck = run_command(
        ["samtools", "quickcheck", "-v", str(path)],
        check=False,
    )

    if quickcheck.returncode != 0:
        raise RuntimeError(
            f"CRAM is present but failed samtools quickcheck: {path}\n"
            f"{quickcheck.stderr}"
        )

    result = run_command(
        [
            "samtools",
            "view",
            "-T",
            str(reference_fasta),
            str(path),
        ],
        check=False,
    )

    if result.returncode != 0:
        raise RuntimeError(
            f"CRAM could not be read: {path}\n{result.stderr}"
        )

    total_reads = 0
    hp1_reads = 0
    hp2_reads = 0

    for line in result.stdout.splitlines():
        if not line:
            continue

        total_reads += 1
        fields = line.split("\t")
        tags = fields[11:]

        if "HP:i:1" in tags:
            hp1_reads += 1
        elif "HP:i:2" in tags:
            hp2_reads += 1

    haplotagged_reads = hp1_reads + hp2_reads
    haplotagged_rate = (
        100.0 * haplotagged_reads / total_reads
        if total_reads > 0
        else 0.0
    )

    row.update(
        {
            "hp1_reads": hp1_reads,
            "hp2_reads": hp2_reads,
            "haplotagged_reads": haplotagged_reads,
            "total_reads": total_reads,
            "haplotagged_pct": f"{haplotagged_rate:.4f}",
        }
    )


def read_samples(path: Path) -> Iterable[tuple[str, Path]]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")

        required = {"sample_id", "output_directory"}

        if not reader.fieldnames or not required.issubset(reader.fieldnames):
            raise ValueError(
                "Sample sheet must contain sample_id and output_directory."
            )

        for record in reader:
            sample_id = record["sample_id"].strip()
            directory = record["output_directory"].strip()

            if sample_id:
                yield sample_id, Path(directory)


def main() -> int:
    args = parse_arguments()

    check_tool("bcftools")
    check_tool("samtools")

    if not args.sample_sheet.is_file():
        raise FileNotFoundError(
            f"Sample sheet not found: {args.sample_sheet}"
        )

    if not args.reference_fasta.is_file():
        raise FileNotFoundError(
            f"Reference FASTA not found: {args.reference_fasta}"
        )

    args.output_tsv.parent.mkdir(parents=True, exist_ok=True)

    rows: list[dict[str, str | int]] = []

    for sample_id, output_directory in read_samples(args.sample_sheet):
        print(f"Processing {sample_id}...", file=sys.stderr)

        row = default_row(sample_id)

        sv_vcf = output_directory / f"{sample_id}.wf_sv.vcf.gz"
        snp_vcf = output_directory / f"{sample_id}.wf_snp.vcf.gz"
        cram = output_directory / f"{sample_id}.haplotagged.cram"

        if sv_vcf.is_file():
            summarize_sv(sv_vcf, row)
        else:
            print(
                f"WARNING: Missing SV VCF for {sample_id}; "
                "SV metrics set to NA.",
                file=sys.stderr,
            )

        if snp_vcf.is_file():
            summarize_snp_phasing(snp_vcf, row)
        else:
            print(
                f"WARNING: Missing SNP VCF for {sample_id}; "
                "small-variant phasing metrics set to NA.",
                file=sys.stderr,
            )

        if cram.is_file():
            summarize_haplotagged_cram(
                cram,
                args.reference_fasta,
                row,
            )
        else:
            print(
                f"WARNING: Missing haplotagged CRAM for {sample_id}; "
                "read-level metrics set to NA.",
                file=sys.stderr,
            )

        rows.append(row)

    with args.output_tsv.open(
        "w",
        encoding="utf-8",
        newline="",
    ) as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=OUTPUT_COLUMNS,
            delimiter="\t",
            lineterminator="\n",
        )
        writer.writeheader()
        writer.writerows(rows)

    print(
        f"Exploratory SV and phasing summary created: "
        f"{args.output_tsv}",
        file=sys.stderr,
    )

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
