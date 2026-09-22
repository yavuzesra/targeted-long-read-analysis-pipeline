#!/usr/bin/env python3

"""Read one resource path from the project YAML configuration."""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import yaml


def parse_arguments() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Read one resource from project.yaml."
    )
    parser.add_argument(
        "--project-config",
        required=True,
        type=Path,
    )
    parser.add_argument(
        "--key",
        required=True,
    )
    return parser.parse_args()


def main() -> int:
    args = parse_arguments()

    if not args.project_config.is_file():
        print(
            f"[ERROR] Project configuration not found: "
            f"{args.project_config}",
            file=sys.stderr,
        )
        return 1

    with args.project_config.open(
        "r",
        encoding="utf-8",
    ) as handle:
        config = yaml.safe_load(handle)

    if not isinstance(config, dict):
        print(
            "[ERROR] Project configuration is not a YAML mapping.",
            file=sys.stderr,
        )
        return 1

    resources = config.get("resources")

    if not isinstance(resources, dict):
        print(
            "[ERROR] resources section is missing.",
            file=sys.stderr,
        )
        return 1

    value = resources.get(args.key)

    if not isinstance(value, str) or not value.strip():
        print(
            f"[ERROR] Resource is missing or empty: {args.key}",
            file=sys.stderr,
        )
        return 1

    print(value.strip())
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
