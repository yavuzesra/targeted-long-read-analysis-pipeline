#!/usr/bin/env bash

# Run optional GIAB truth-set sensitivity analyses separately from run_all.sh.

set -Eeuo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(cd "$SCRIPT_DIR/.." && pwd)"

# shellcheck source=scripts/lib/common.sh
source "$SCRIPT_DIR/lib/common.sh"

PROJECT_CONFIG="$PROJECT_ROOT/config/project.yaml"
SAMPLES_FILE="$PROJECT_ROOT/config/samples.tsv"
SENSITIVITY_CONFIG="$PROJECT_ROOT/config/giab_sensitivity.tsv"
PRIMARY_SUMMARY="$PROJECT_ROOT/results/validation_control/giab_benchmark/giab_small_variant_benchmark_summary.tsv"
PRIMARY_MANIFEST="$PROJECT_ROOT/work/intermediate/giab_validation/prepared_inputs/giab_input_manifest.tsv"

usage() {
    cat <<'EOF'
Usage:
  bash scripts/run_giab_truthset_sensitivity.sh [options]

Options:
  --config PATH              Path to project.yaml
  --samples PATH             Path to primary samples.tsv
  --sensitivity-config PATH  Path to giab_sensitivity.tsv
  --primary-summary PATH     Existing primary GIAB summary TSV
  --primary-manifest PATH    Existing primary prepared-input manifest TSV
  -h, --help                 Show this help message
EOF
}

while [[ $# -gt 0 ]]; do
    case "$1" in
        --config) PROJECT_CONFIG="$2"; shift 2 ;;
        --samples) SAMPLES_FILE="$2"; shift 2 ;;
        --sensitivity-config) SENSITIVITY_CONFIG="$2"; shift 2 ;;
        --primary-summary) PRIMARY_SUMMARY="$2"; shift 2 ;;
        --primary-manifest) PRIMARY_MANIFEST="$2"; shift 2 ;;
        -h|--help) usage; exit 0 ;;
        *) die "Unknown argument: $1" ;;
    esac
done

require_command conda
require_file "$PROJECT_CONFIG" "Project configuration"
require_file "$SAMPLES_FILE" "Primary sample metadata"
require_file "$SENSITIVITY_CONFIG" "GIAB sensitivity configuration"
require_file "$PRIMARY_SUMMARY" "Primary GIAB summary"
require_file "$PRIMARY_MANIFEST" "Primary GIAB input manifest"

PYTHON_ENVIRONMENT="human-variation-workflow"
WORK_ROOT="$PROJECT_ROOT/work/intermediate/giab_truthset_sensitivity"
PREPARED_ROOT="$WORK_ROOT/prepared_inputs"
TEMP_SAMPLES="$WORK_ROOT/sensitivity_samples.tsv"
RESULTS_ROOT="$PROJECT_ROOT/results/validation_control/giab_truthset_sensitivity"
SENSITIVITY_SUMMARY="$RESULTS_ROOT/giab_sensitivity_benchmark_summary.tsv"
COMPARISON_OUTPUT="$RESULTS_ROOT/giab_truthset_comparison.tsv"

TARGET_BED=$(conda run -n "$PYTHON_ENVIRONMENT" python \
    "$SCRIPT_DIR/python/get_project_resource.py" \
    --project-config "$PROJECT_CONFIG" --key target_bed)
REFERENCE_FASTA=$(conda run -n "$PYTHON_ENVIRONMENT" python \
    "$SCRIPT_DIR/python/get_project_resource.py" \
    --project-config "$PROJECT_CONFIG" --key reference_fasta)
RTG_SDF=$(conda run -n "$PYTHON_ENVIRONMENT" python \
    "$SCRIPT_DIR/python/get_project_resource.py" \
    --project-config "$PROJECT_CONFIG" --key rtg_sdf)

create_directory "$WORK_ROOT"

conda run -n "$PYTHON_ENVIRONMENT" python \
    "$SCRIPT_DIR/python/prepare_giab_sensitivity_metadata.py" \
    --samples "$SAMPLES_FILE" \
    --sensitivity-config "$SENSITIVITY_CONFIG" \
    --output "$TEMP_SAMPLES"

rm -rf "$PREPARED_ROOT" "$RESULTS_ROOT"
create_directory "$RESULTS_ROOT"

conda run -n giab-benchmark bash \
    "$SCRIPT_DIR/run_prepare_giab_inputs.sh" \
    --target-bed "$TARGET_BED" \
    --samples "$TEMP_SAMPLES" \
    --output-root "$PREPARED_ROOT"

conda run -n giab-benchmark bash \
    "$SCRIPT_DIR/run_giab_benchmark.sh" \
    --input-manifest "$PREPARED_ROOT/giab_input_manifest.tsv" \
    --reference-fasta "$REFERENCE_FASTA" \
    --rtg-sdf "$RTG_SDF" \
    --output-root "$RESULTS_ROOT"

conda run -n "$PYTHON_ENVIRONMENT" python \
    "$SCRIPT_DIR/python/summarize_giab_benchmark.py" \
    --run-table "$RESULTS_ROOT/giab_benchmark_runs.tsv" \
    --output "$SENSITIVITY_SUMMARY"

conda run -n "$PYTHON_ENVIRONMENT" python \
    "$SCRIPT_DIR/python/compare_giab_truthsets.py" \
    --sensitivity-config "$SENSITIVITY_CONFIG" \
    --primary-summary "$PRIMARY_SUMMARY" \
    --primary-manifest "$PRIMARY_MANIFEST" \
    --sensitivity-summary "$SENSITIVITY_SUMMARY" \
    --sensitivity-manifest "$PREPARED_ROOT/giab_input_manifest.tsv" \
    --output "$COMPARISON_OUTPUT"

log_info "GIAB truth-set sensitivity analysis completed."
log_info "Sensitivity summary: $SENSITIVITY_SUMMARY"
log_info "Primary-versus-sensitivity comparison: $COMPARISON_OUTPUT"
