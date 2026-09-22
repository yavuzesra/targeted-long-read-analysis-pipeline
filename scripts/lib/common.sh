#!/usr/bin/env bash

# Shared Bash functions for the QIAseq xHYB Long Read pipeline.

set -Eeuo pipefail

log_info() {
    printf '[INFO] %s\n' "$*"
}

log_warning() {
    printf '[WARNING] %s\n' "$*" >&2
}

log_error() {
    printf '[ERROR] %s\n' "$*" >&2
}

die() {
    log_error "$*"
    exit 1
}

require_file() {
    local file_path="$1"
    local label="$2"

    [[ -f "$file_path" ]] || die "$label not found: $file_path"
}

require_command() {
    local command_name="$1"

    command -v "$command_name" >/dev/null 2>&1 ||
        die "Required command not found: $command_name"
}

create_directory() {
    local directory_path="$1"
    mkdir -p "$directory_path"
}

timestamp() {
    date '+%Y-%m-%d %H:%M:%S'
}
