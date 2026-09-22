#!/usr/bin/env python3
"""
Prepare IGV batch files for complete target-region inspection.

The script does not launch IGV. It creates:
- one IGV batch file per configured sample group;
- one locus manifest containing all target regions;
- one snapshot directory per sample group.

Coordinate conversion:
- BED coordinates are 0-based and half-open.
- IGV loci are 1-based and inclusive.
- IGV start = BED start + 1.
- IGV end = BED end.

This module is intended for technical visual inspection and does not perform
variant interpretation.
"""

from __future__ import annotations

import argparse
import csv
import re
import sys
from collections import defaultdict
from pathlib import Path


TRUE_VALUES = {"true", "yes", "1"}


def parse_arguments() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Prepare IGV batch files for complete inspection of target BED "
            "regions."
        )
    )
    parser.add_argument(
        "--samples",
        required=True,
        type=Path,
        help="Selected or complete samples.tsv file.",
    )
    parser.add_argument(
        "--groups",
        required=True,
        type=Path,
        help="IGV group configuration TSV.",
    )
    parser.add_argument(
        "--target-bed",
        required=True,
        type=Path,
        help="Target BED file in 0-based half-open format.",
    )
    parser.add_argument(
        "--reference-fasta",
        required=True,
        type=Path,
        help="Reference FASTA loaded by IGV.",
    )
    parser.add_argument(
        "--analysis-role",
        required=True,
        help="Analysis role to include, for example method_development.",
    )
    parser.add_argument(
        "--output-root",
        required=True,
        type=Path,
        help="Root directory for IGV batch files and snapshot directories.",
    )
    parser.add_argument(
        "--gene-annotation",
        type=Path,
        default=None,
        help=(
            "Optional gene-to-target annotation TSV. When supplied, gene "
            "symbols are included in snapshot filenames."
        ),
    )
    parser.add_argument(
        "--sleep-interval",
        type=int,
        default=1000,
        help="IGV batch sleep interval in milliseconds. Default: 1000.",
    )
    parser.add_argument(
        "--max-panel-height",
        type=int,
        default=700,
        help="IGV maximum panel height. Default: 700.",
    )
    return parser.parse_args()


def read_tsv(path: Path) -> list[dict[str, str]]:
    if not path.is_file():
        raise FileNotFoundError(f"TSV file not found: {path}")

    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames is None:
            raise ValueError(f"TSV file has no header: {path}")
        return list(reader)


def require_columns(
    rows: list[dict[str, str]],
    required: set[str],
    path: Path,
) -> None:
    if not rows:
        raise ValueError(f"TSV file contains no data rows: {path}")

    available = set(rows[0])
    missing = sorted(required - available)

    if missing:
        raise ValueError(
            f"Missing required column(s) in {path}: {', '.join(missing)}"
        )


def is_included(value: str | None) -> bool:
    return value is not None and value.strip().lower() in TRUE_VALUES


def safe_filename(value: str) -> str:
    cleaned = re.sub(r"[^A-Za-z0-9_.-]+", "_", value.strip())
    return cleaned.strip("_") or "unnamed"


def load_gene_labels(
    annotation_path: Path | None,
) -> dict[tuple[str, int, int], str]:
    if annotation_path is None:
        return {}

    rows = read_tsv(annotation_path)

    possible_gene_columns = ("gene_symbol", "gene", "Gene")
    possible_chrom_columns = ("chromosome", "chromosome")
    possible_start_columns = ("bed_start_0based", "start")
    possible_end_columns = ("bed_end_0based_exclusive", "end")

    def find_column(options: tuple[str, ...]) -> str | None:
        fields = set(rows[0])
        for option in options:
            if option in fields:
                return option
        return None

    gene_column = find_column(possible_gene_columns)
    chrom_column = find_column(possible_chrom_columns)
    start_column = find_column(possible_start_columns)
    end_column = find_column(possible_end_columns)

    if not all((gene_column, chrom_column, start_column, end_column)):
        print(
            "[WARNING] Gene annotation columns were not recognised. "
            "Snapshot filenames will use interval labels only.",
            file=sys.stderr,
        )
        return {}

    labels: dict[tuple[str, int, int], list[str]] = defaultdict(list)

    for row in rows:
        chrom = row[chrom_column].strip()
        start = int(row[start_column])
        end = int(row[end_column])
        gene = row[gene_column].strip()

        if gene and gene not in labels[(chrom, start, end)]:
            labels[(chrom, start, end)].append(gene)

    return {
        key: "_".join(sorted(genes))
        for key, genes in labels.items()
    }


def read_bed(
    bed_path: Path,
    gene_labels: dict[tuple[str, int, int], str],
) -> list[dict[str, str | int]]:
    if not bed_path.is_file():
        raise FileNotFoundError(f"Target BED not found: {bed_path}")

    intervals: list[dict[str, str | int]] = []

    with bed_path.open("r", encoding="utf-8") as handle:
        for line_number, raw_line in enumerate(handle, start=1):
            line = raw_line.strip()

            if not line or line.startswith("#"):
                continue

            fields = line.split("\t")

            if len(fields) < 3:
                raise ValueError(
                    f"Invalid BED row at line {line_number}: {raw_line.rstrip()}"
                )

            chrom = fields[0]
            start_0based = int(fields[1])
            end_halfopen = int(fields[2])

            if start_0based < 0 or end_halfopen <= start_0based:
                raise ValueError(
                    f"Invalid BED coordinates at line {line_number}: "
                    f"{chrom}:{start_0based}-{end_halfopen}"
                )

            gene_label = gene_labels.get(
                (chrom, start_0based, end_halfopen),
                "",
            )

            intervals.append(
                {
                    "interval_number": len(intervals) + 1,
                    "chromosome": chrom,
                    "bed_start_0based": start_0based,
                    "bed_end_0based_exclusive": end_halfopen,
                    "igv_start_1based": start_0based + 1,
                    "igv_end_1based": end_halfopen,
                    "gene_label": gene_label,
                }
            )

    if not intervals:
        raise ValueError(f"No target intervals found in BED file: {bed_path}")

    return intervals


def write_locus_manifest(
    intervals: list[dict[str, str | int]],
    output_path: Path,
) -> None:
    fieldnames = [
        "interval_number",
        "gene_label",
        "chromosome",
        "bed_start_0based",
        "bed_end_halfopen",
        "igv_start_1based",
        "igv_end_1based",
        "igv_locus",
    ]

    output_path.parent.mkdir(parents=True, exist_ok=True)

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

        for interval in intervals:
            chrom = str(interval["chromosome"])
            igv_start = int(interval["igv_start_1based"])
            igv_end = int(interval["igv_end_1based"])

            writer.writerow(
                {
                    "interval_number": interval["interval_number"],
                    "gene_label": interval["gene_label"],
                    "chromosome": chrom,
                    "bed_start_0based": interval["bed_start_0based"],
                    "bed_end_halfopen": interval["bed_end_0based_exclusive"],
                    "igv_start_1based": igv_start,
                    "igv_end_1based": igv_end,
                    "igv_locus": f"{chrom}:{igv_start}-{igv_end}",
                }
            )


def main() -> int:
    args = parse_arguments()

    sample_rows = read_tsv(args.samples)
    group_rows = read_tsv(args.groups)

    require_columns(
        sample_rows,
        {"sample_id", "analysis_role", "bam_path", "bai_path", "include"},
        args.samples,
    )
    require_columns(
        group_rows,
        {"igv_group", "analysis_role", "sample_id", "include"},
        args.groups,
    )

    if not args.reference_fasta.is_file():
        raise FileNotFoundError(
            f"Reference FASTA not found: {args.reference_fasta}"
        )

    samples_by_id = {
        row["sample_id"].strip(): row
        for row in sample_rows
        if is_included(row.get("include"))
    }

    requested_groups: dict[str, list[str]] = defaultdict(list)

    for row in group_rows:
        if not is_included(row.get("include")):
            continue

        if row["analysis_role"].strip() != args.analysis_role:
            continue

        group_name = row["igv_group"].strip()
        sample_id = row["sample_id"].strip()

        if not group_name:
            raise ValueError("An IGV group name is empty.")

        requested_groups[group_name].append(sample_id)

    if not requested_groups:
        raise ValueError(
            f"No IGV groups were found for analysis role: "
            f"{args.analysis_role}"
        )

    gene_labels = load_gene_labels(args.gene_annotation)
    intervals = read_bed(args.target_bed, gene_labels)

    batches_dir = args.output_root / "batch_scripts"
    snapshots_root = args.output_root / "snapshots"
    manifests_dir = args.output_root / "manifests"

    batches_dir.mkdir(parents=True, exist_ok=True)
    snapshots_root.mkdir(parents=True, exist_ok=True)
    manifests_dir.mkdir(parents=True, exist_ok=True)

    locus_manifest = manifests_dir / "target_region_loci.tsv"
    write_locus_manifest(intervals, locus_manifest)

    batch_manifest_rows: list[dict[str, str | int]] = []

    for group_name, sample_ids in sorted(requested_groups.items()):
        safe_group = safe_filename(group_name)
        snapshot_dir = snapshots_root / safe_group
        batch_path = batches_dir / f"{safe_group}.target_regions.batch.txt"

        snapshot_dir.mkdir(parents=True, exist_ok=True)

        bam_paths: list[Path] = []

        for sample_id in sample_ids:
            if sample_id not in samples_by_id:
                raise ValueError(
                    f"IGV group '{group_name}' references unavailable or "
                    f"excluded sample: {sample_id}"
                )

            sample = samples_by_id[sample_id]

            if sample["analysis_role"].strip() != args.analysis_role:
                raise ValueError(
                    f"Sample '{sample_id}' does not have analysis role "
                    f"'{args.analysis_role}'."
                )

            bam_path = Path(sample["bam_path"])
            bai_path = Path(sample["bai_path"])

            if not bam_path.is_file():
                raise FileNotFoundError(
                    f"BAM file not found for {sample_id}: {bam_path}"
                )

            if not bai_path.is_file():
                raise FileNotFoundError(
                    f"BAM index not found for {sample_id}: {bai_path}"
                )

            bam_paths.append(bam_path)

        with batch_path.open("w", encoding="utf-8") as handle:
            handle.write("new\n")
            handle.write(f"genome {args.reference_fasta}\n")
            handle.write(f"snapshotDirectory {snapshot_dir}\n")
            handle.write(f"maxPanelHeight {args.max_panel_height}\n")
            handle.write(f"setSleepInterval {args.sleep_interval}\n")

            for bam_path in bam_paths:
                handle.write(f"load {bam_path}\n")

            handle.write("collapse\n")

            for interval in intervals:
                number = int(interval["interval_number"])
                chrom = str(interval["chromosome"])
                igv_start = int(interval["igv_start_1based"])
                igv_end = int(interval["igv_end_1based"])
                gene_label = str(interval["gene_label"])

                if gene_label:
                    label = safe_filename(gene_label)
                else:
                    label = f"interval_{number:02d}"

                snapshot_name = (
                    f"{number:02d}_{label}_{chrom}_"
                    f"{igv_start}_{igv_end}.png"
                )

                handle.write(
                    f"goto {chrom}:{igv_start}-{igv_end}\n"
                )
                handle.write("sort base\n")
                handle.write(f"snapshot {snapshot_name}\n")

            handle.write("exit\n")

        batch_manifest_rows.append(
            {
                "igv_group": group_name,
                "analysis_role": args.analysis_role,
                "sample_count": len(sample_ids),
                "sample_ids": ",".join(sample_ids),
                "batch_file": str(batch_path.resolve()),
                "snapshot_directory": str(snapshot_dir.resolve()),
                "target_interval_count": len(intervals),
            }
        )

    batch_manifest_path = manifests_dir / "igv_batch_manifest.tsv"
    fieldnames = [
        "igv_group",
        "analysis_role",
        "sample_count",
        "sample_ids",
        "batch_file",
        "snapshot_directory",
        "target_interval_count",
    ]

    with batch_manifest_path.open(
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
        writer.writerows(batch_manifest_rows)

    print(f"[INFO] Analysis role: {args.analysis_role}")
    print(f"[INFO] IGV groups prepared: {len(batch_manifest_rows)}")
    print(f"[INFO] Target intervals per group: {len(intervals)}")
    print(f"[INFO] Batch directory: {batches_dir.resolve()}")
    print(f"[INFO] Snapshot root: {snapshots_root.resolve()}")
    print(f"[INFO] Locus manifest: {locus_manifest.resolve()}")
    print(f"[INFO] Batch manifest: {batch_manifest_path.resolve()}")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
