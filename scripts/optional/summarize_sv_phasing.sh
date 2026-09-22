#!/usr/bin/env bash

set -euo pipefail

# Thin wrapper for the exploratory SV and phasing Python summarization utility.

SCRIPT_DIR="$(
    cd -- "$(dirname -- "${BASH_SOURCE[0]}")"
    pwd
)"

exec \
    python \
    "$SCRIPT_DIR/summarize_sv_phasing.py" \
    "$@"
