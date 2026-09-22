#!/usr/bin/env python3

"""
Select included samples by run ID and/or analysis role.

The original samples.tsv file is never modified. A new filtered TSV file is
written for the current pipeline execution.

Selection logic:
    - include must already be true in the source metadata;
    - --run-id restricts samples to one or more run IDs;
    - --role restricts samples to one or more analysis roles;
    - when both are supplied, both conditions must be satisfied.
"""

from __future__ import annotations

import argparse
import csv
import sys
from pathlib import Path


TRUE_VALUES = {"true", "yes", "1"}

ALLOWED_ROLES = {
    "method_development",
    "validation_control",
    "workflow_application",
}


def parse_arguments() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Select samples for one pipeline execution."
    )
    parser.add_argument(
        "--samples",
        required=True,
        type=Path,
        help="Source samples.tsv file",
    )
    parser.add_argument(
        "--output",
        required=True,
        type=Path,
        help="Filtered output samples.tsv file",
    )
    parser.add_argument(
        "--run-id",
        action="append",
        default=[],
        help=(
            "Run ID to include. May be supplied more than once."
        ),
    )
    parser.add_argument(
        "--role",
        action="append",
        default=[],
        help=(
            "Analysis role to include. May be supplied more than once."
        ),
    )
    return parser.parse_args()


def is_included(value: str) -> bool:
    return value.strip().lower() in TRUE_VALUES


def read_samples(
    path: Path,
) -> tuple[list[str], list[dict[str, str]]]:
    if not path.is_file():
        raise ValueError(f"Sample metadata not found: {path}")

    with path.open(
        "r",
        encoding="utf-8",
        newline="",
    ) as handle:
        reader = csv.DictReader(handle, delimiter="\t")

        if reader.fieldnames is None:
            raise ValueError("samples.tsv has no header.")

        required = {
            "sample_id",
            "analysis_role",
            "run_id",
            "include",
        }

        missing = required - set(reader.fieldnames)

        if missing:
            raise ValueError(
                "samples.tsv is missing required columns: "
                + ", ".join(sorted(missing))
            )

        rows = list(reader)

    return list(reader.fieldnames), rows


def validate_requested_roles(roles: list[str]) -> None:
    invalid = sorted(set(roles) - ALLOWED_ROLES)

    if invalid:
        raise ValueError(
            "Invalid requested analysis role(s): "
            + ", ".join(invalid)
        )


def select_rows(
    rows: list[dict[str, str]],
    run_ids: set[str],
    roles: set[str],
) -> list[dict[str, str]]:
    selected: list[dict[str, str]] = []

    for row in rows:
        if not is_included(row["include"]):
            continue

        if run_ids and row["run_id"].strip() not in run_ids:
            continue

        if roles and row["analysis_role"].strip() not in roles:
            continue

        selected.append(row)

    return selected


def write_samples(
    rows: list[dict[str, str]],
    fieldnames: list[str],
    output_path: Path,
) -> None:
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
        writer.writerows(rows)


def main() -> int:
    args = parse_arguments()

    try:
        validate_requested_roles(args.role)

        fieldnames, rows = read_samples(args.samples)

        selected = select_rows(
            rows=rows,
            run_ids=set(args.run_id),
            roles=set(args.role),
        )

        if not selected:
            raise ValueError(
                "No enabled samples matched the requested selection."
            )

        write_samples(
            rows=selected,
            fieldnames=fieldnames,
            output_path=args.output,
        )

    except ValueError as exc:
        print(f"[ERROR] {exc}", file=sys.stderr)
        return 1

    print(f"[INFO] Selected samples: {len(selected)}")
    print(
        "[INFO] Selected sample IDs: "
        + ", ".join(row["sample_id"] for row in selected)
    )
    print(f"[INFO] Output: {args.output}")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
