#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import re
import sys
from pathlib import Path


REQUIRED_REVIEW_COLUMNS = {
    "sample_id",
    "chromosome",
    "position_1based",
    "ref",
    "alt",
}

OPTIONAL_LABEL_COLUMNS = [
    "gene",
    "genotype",
    "clinvar_significance",
    "consequence",
]

BAM_COLUMN_CANDIDATES = [
    "bam",
    "bam_path",
    "input_bam",
    "alignment_bam",
    "sorted_aligned_bam",
    "analysis_bam",
    "primary_bam",
]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Prepare IGV batch files for review-candidate variants."
    )
    parser.add_argument("--samples", required=True, type=Path)
    parser.add_argument("--project-root", required=True, type=Path)
    parser.add_argument("--reference-fasta", required=True, type=Path)
    parser.add_argument("--analysis-role", default="workflow_application")
    parser.add_argument("--output-root", required=True, type=Path)
    parser.add_argument("--window-bp", type=int, default=100)
    return parser.parse_args()


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        return list(reader)


def truthy(value: str) -> bool:
    return str(value).strip().lower() not in {"", "0", "false", "no", "n"}


def sanitize_label(value: str) -> str:
    cleaned = re.sub(r"[^A-Za-z0-9._-]+", "_", str(value).strip())
    cleaned = cleaned.strip("._")
    return cleaned or "NA"


def normalize_path(raw_value: str, base_dir: Path) -> Path:
    path = Path(raw_value).expanduser()
    if path.is_absolute():
        return path
    return (base_dir / path).resolve()


def find_bam_path(row: dict[str, str], base_dir: Path) -> Path:
    for column in BAM_COLUMN_CANDIDATES:
        value = row.get(column, "").strip()
        if value and value.upper() != "NA":
            return normalize_path(value, base_dir)
    raise KeyError(
        "No BAM path column found. Tried: " + ", ".join(BAM_COLUMN_CANDIDATES)
    )


def build_locus(chromosome: str, position_1based: int, ref: str, alt: str, window_bp: int) -> str:
    allele_span = max(len(ref or ""), len(alt or ""), 1)
    start = max(1, position_1based - window_bp)
    end = position_1based + allele_span + window_bp - 1
    return f"{chromosome}:{start}-{end}"


def main() -> int:
    args = parse_args()

    samples_path = args.samples.resolve()
    project_root = args.project_root.resolve()
    reference_fasta = args.reference_fasta.resolve()
    output_root = args.output_root.resolve() / args.analysis_role

    if not samples_path.exists():
        print(f"ERROR: Samples file not found: {samples_path}", file=sys.stderr)
        return 1

    if not reference_fasta.exists():
        print(f"ERROR: Reference FASTA not found: {reference_fasta}", file=sys.stderr)
        return 1

    batch_dir = output_root / "batch_scripts"
    manifest_dir = output_root / "manifests"
    snapshots_root = output_root / "snapshots"

    batch_dir.mkdir(parents=True, exist_ok=True)
    manifest_dir.mkdir(parents=True, exist_ok=True)
    snapshots_root.mkdir(parents=True, exist_ok=True)

    batch_manifest_path = manifest_dir / "igv_review_batch_manifest.tsv"
    loci_manifest_path = manifest_dir / "review_candidate_loci.tsv"

    sample_rows = read_tsv(samples_path)

    manifest_headers = [
        "sample_id",
        "analysis_role",
        "bam_path",
        "review_candidates_tsv",
        "batch_script",
        "snapshot_directory",
        "variant_count",
    ]

    loci_headers = [
        "sample_id",
        "analysis_role",
        "variant_index",
        "chromosome",
        "position_1based",
        "ref",
        "alt",
        "locus",
        "snapshot_name",
        "gene",
        "genotype",
        "clinvar_significance",
        "consequence",
        "review_candidates_tsv",
        "bam_path",
    ]

    prepared_samples = 0
    total_loci = 0

    with batch_manifest_path.open("w", encoding="utf-8", newline="") as manifest_handle, \
         loci_manifest_path.open("w", encoding="utf-8", newline="") as loci_handle:

        manifest_writer = csv.DictWriter(
            manifest_handle,
            fieldnames=manifest_headers,
            delimiter="\t",
        )
        loci_writer = csv.DictWriter(
            loci_handle,
            fieldnames=loci_headers,
            delimiter="\t",
        )

        manifest_writer.writeheader()
        loci_writer.writeheader()

        for sample_row in sample_rows:
            sample_id = sample_row.get("sample_id", "").strip()
            analysis_role = sample_row.get("analysis_role", "").strip()
            include_value = sample_row.get("include", "true")

            if not sample_id:
                continue

            if analysis_role != args.analysis_role:
                continue

            if not truthy(include_value):
                continue

            try:
                bam_path = find_bam_path(sample_row, project_root)
            except KeyError as exc:
                print(f"ERROR: {sample_id}: {exc}", file=sys.stderr)
                return 1

            review_candidates_tsv = (
                project_root
                / "results"
                / "workflow_application"
                / sample_id
                / "small_variants"
                / "variants"
                / f"{sample_id}.review_candidates.tsv"
            )

            if not review_candidates_tsv.exists():
                print(
                    f"WARNING: Review-candidate TSV not found for {sample_id}: {review_candidates_tsv}",
                    file=sys.stderr,
                )
                continue

            with review_candidates_tsv.open("r", encoding="utf-8", newline="") as handle:
                reader = csv.DictReader(handle, delimiter="\t")
                fieldnames = set(reader.fieldnames or [])
                missing = REQUIRED_REVIEW_COLUMNS - fieldnames
                if missing:
                    print(
                        f"ERROR: {sample_id}: review_candidates TSV missing columns: {sorted(missing)}",
                        file=sys.stderr,
                    )
                    return 1

                variant_rows = list(reader)

            if not variant_rows:
                print(
                    f"WARNING: {sample_id}: review_candidates TSV is empty, skipping.",
                    file=sys.stderr,
                )
                continue

            sample_snapshot_dir = snapshots_root / sample_id
            sample_snapshot_dir.mkdir(parents=True, exist_ok=True)

            batch_script_path = batch_dir / f"{sample_id}.review_candidates.batch.txt"

            batch_lines = [
                "new",
                f"genome {reference_fasta}",
                f"load {bam_path}",
                f"snapshotDirectory {sample_snapshot_dir}",
                "maxPanelHeight 2000",
                "collapse",
            ]

            variant_count = 0

            for idx, variant in enumerate(variant_rows, start=1):
                chromosome = variant["chromosome"].strip()
                position_text = variant["position_1based"].strip()
                ref = variant.get("ref", "").strip()
                alt = variant.get("alt", "").strip()

                try:
                    position_1based = int(position_text)
                except ValueError:
                    print(
                        f"WARNING: {sample_id}: invalid position '{position_text}', skipping row {idx}.",
                        file=sys.stderr,
                    )
                    continue

                gene = variant.get("gene", "NA").strip() or "NA"
                genotype = variant.get("genotype", "NA").strip() or "NA"
                clinvar_significance = variant.get("clinvar_significance", "NA").strip() or "NA"
                consequence = variant.get("consequence", "NA").strip() or "NA"

                locus = build_locus(
                    chromosome=chromosome,
                    position_1based=position_1based,
                    ref=ref,
                    alt=alt,
                    window_bp=args.window_bp,
                )

                snapshot_name = (
                    f"{idx:03d}_"
                    f"{sanitize_label(gene)}_"
                    f"{sanitize_label(chromosome)}_"
                    f"{position_1based}_"
                    f"{sanitize_label(ref)}_"
                    f"{sanitize_label(alt)}.png"
                )

                batch_lines.extend(
                    [
                        f"goto {locus}",
                        "sort base",
                        f"snapshot {snapshot_name}",
                    ]
                )

                loci_writer.writerow(
                    {
                        "sample_id": sample_id,
                        "analysis_role": analysis_role,
                        "variant_index": idx,
                        "chromosome": chromosome,
                        "position_1based": position_1based,
                        "ref": ref,
                        "alt": alt,
                        "locus": locus,
                        "snapshot_name": snapshot_name,
                        "gene": gene,
                        "genotype": genotype,
                        "clinvar_significance": clinvar_significance,
                        "consequence": consequence,
                        "review_candidates_tsv": str(review_candidates_tsv),
                        "bam_path": str(bam_path),
                    }
                )

                variant_count += 1
                total_loci += 1

            if variant_count == 0:
                print(
                    f"WARNING: {sample_id}: no valid loci were produced, skipping batch file.",
                    file=sys.stderr,
                )
                continue

            batch_lines.append("exit")

            batch_script_path.write_text(
                "\n".join(batch_lines) + "\n",
                encoding="utf-8",
            )

            manifest_writer.writerow(
                {
                    "sample_id": sample_id,
                    "analysis_role": analysis_role,
                    "bam_path": str(bam_path),
                    "review_candidates_tsv": str(review_candidates_tsv),
                    "batch_script": str(batch_script_path),
                    "snapshot_directory": str(sample_snapshot_dir),
                    "variant_count": variant_count,
                }
            )

            prepared_samples += 1
            print(
                f"[INFO] Prepared review-candidate IGV batch: {sample_id} ({variant_count} loci)",
                file=sys.stderr,
            )

    print(
        f"[INFO] Review-candidate IGV preparation completed. "
        f"Samples prepared: {prepared_samples}; loci prepared: {total_loci}",
        file=sys.stderr,
    )
    print(f"[INFO] Batch manifest: {batch_manifest_path}", file=sys.stderr)
    print(f"[INFO] Loci manifest: {loci_manifest_path}", file=sys.stderr)

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
