#!/usr/bin/env bash
set -euo pipefail

PROJECT_ROOT=$(
    cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd
)

PYTHON_ENVIRONMENT=${PYTHON_ENVIRONMENT:-human-variation-workflow}

OUTPUT_DIR="$PROJECT_ROOT/results/shared_qc/sequencing_run_target_metrics"

PER_SAMPLE_OUTPUT="$OUTPUT_DIR/target_metrics_all_samples.tsv"
RUN_OUTPUT="$OUTPUT_DIR/target_metrics_by_sequencing_run.tsv"

mkdir -p "$OUTPUT_DIR"

printf '[INFO] Calculating target metrics for all included samples.\n'
printf '[INFO] Role-based filtering is not applied in this dedicated run-level analysis.\n'

conda run -n "$PYTHON_ENVIRONMENT" \
    python "$PROJECT_ROOT/scripts/python/target_metrics.py" \
    --project-config "$PROJECT_ROOT/config/project.yaml" \
    --samples "$PROJECT_ROOT/config/samples.tsv" \
    --filtering-config "$PROJECT_ROOT/config/filtering.yaml" \
    --output "$PER_SAMPLE_OUTPUT" \
    --run-output "$RUN_OUTPUT"

printf '[INFO] Per-sample output: %s\n' "$PER_SAMPLE_OUTPUT"
printf '[INFO] Sequencing-run output: %s\n' "$RUN_OUTPUT"
