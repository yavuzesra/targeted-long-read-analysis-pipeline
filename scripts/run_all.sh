#!/usr/bin/env bash

# Main entry point for the QIAseq xHYB Long Read analysis pipeline.

set -Eeuo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(cd "$SCRIPT_DIR/.." && pwd)"

# shellcheck source=scripts/lib/common.sh
source "$SCRIPT_DIR/lib/common.sh"

PROJECT_CONFIG="$PROJECT_ROOT/config/project.yaml"
SAMPLES_FILE="$PROJECT_ROOT/config/samples.tsv"
RUNS_FILE="$PROJECT_ROOT/config/runs.tsv"

VALIDATE_ONLY=false
PLAN_ONLY=false
OFF_TARGET_ANALYSIS_REQUESTED=false

REQUESTED_RUN_IDS=()
REQUESTED_ROLES=()

PYTHON_ENVIRONMENT="human-variation-workflow"

usage() {
    cat <<'EOF'
Usage:
  bash scripts/run_all.sh [options]

Options:
  --config PATH          Path to project.yaml
  --samples PATH         Path to samples.tsv
  --runs PATH            Path to runs.tsv
  --run-id ID            Restrict analysis to one run ID
                         May be supplied more than once
  --role ROLE            Restrict analysis to one analysis role
                         May be supplied more than once
  --validate-only        Validate inputs and stop
  --plan-only            Validate inputs, print execution plan, and stop
  --off-target-analysis  Run the optional exploratory base-level off-target analysis
  -h, --help             Show this help message

The selected role determines which implemented analysis modules are run.
EOF
}

while [[ $# -gt 0 ]]; do
    case "$1" in
        --config)
            [[ $# -ge 2 ]] || die "--config requires a path."
            PROJECT_CONFIG="$2"
            shift 2
            ;;
        --samples)
            [[ $# -ge 2 ]] || die "--samples requires a path."
            SAMPLES_FILE="$2"
            shift 2
            ;;
        --runs)
            [[ $# -ge 2 ]] || die "--runs requires a path."
            RUNS_FILE="$2"
            shift 2
            ;;
        --run-id)
            [[ $# -ge 2 ]] || die "--run-id requires a value."
            REQUESTED_RUN_IDS+=("$2")
            shift 2
            ;;
        --role)
            [[ $# -ge 2 ]] || die "--role requires a value."
            REQUESTED_ROLES+=("$2")
            shift 2
            ;;
        --validate-only)
            VALIDATE_ONLY=true
            shift
            ;;
        --plan-only)
            PLAN_ONLY=true
            shift
            ;;
        --off-target-analysis)
            OFF_TARGET_ANALYSIS_REQUESTED=true
            shift
            ;;

        -h|--help)
            usage
            exit 0
            ;;
        *)
            die "Unknown argument: $1"
            ;;
    esac
done

require_command conda
require_file "$PROJECT_CONFIG" "Project configuration"
require_file "$SAMPLES_FILE" "Sample metadata"
require_file "$RUNS_FILE" "Run metadata"

LOG_DIR="$PROJECT_ROOT/work/logs"
create_directory "$LOG_DIR"

RUN_TIMESTAMP="$(date '+%Y%m%d_%H%M%S')"
LOG_FILE="$LOG_DIR/pipeline_${RUN_TIMESTAMP}.log"

exec > >(tee -a "$LOG_FILE") 2>&1

log_info "Pipeline started: $(timestamp)"
log_info "Project root: $PROJECT_ROOT"
log_info "Project config: $PROJECT_CONFIG"
log_info "Source samples file: $SAMPLES_FILE"
log_info "Runs file: $RUNS_FILE"
log_info "Log file: $LOG_FILE"
echo

SELECTED_SAMPLES_FILE="$PROJECT_ROOT/work/tmp/selected_samples_${RUN_TIMESTAMP}.tsv"
create_directory "$PROJECT_ROOT/work/tmp"

SELECTION_ARGUMENTS=()

for run_id in "${REQUESTED_RUN_IDS[@]}"; do
    SELECTION_ARGUMENTS+=(--run-id "$run_id")
done

for role in "${REQUESTED_ROLES[@]}"; do
    SELECTION_ARGUMENTS+=(--role "$role")
done

if ((
    ${#REQUESTED_RUN_IDS[@]} > 0 ||
    ${#REQUESTED_ROLES[@]} > 0
)); then
    log_info "Selecting samples for this execution."

    conda run -n "$PYTHON_ENVIRONMENT" \
        python "$SCRIPT_DIR/python/select_samples.py" \
        --samples "$SAMPLES_FILE" \
        --output "$SELECTED_SAMPLES_FILE" \
        "${SELECTION_ARGUMENTS[@]}"

    ACTIVE_SAMPLES_FILE="$SELECTED_SAMPLES_FILE"
else
    ACTIVE_SAMPLES_FILE="$SAMPLES_FILE"
fi

log_info "Active samples file: $ACTIVE_SAMPLES_FILE"
echo

log_info "Validating project configuration and metadata."

conda run -n "$PYTHON_ENVIRONMENT" \
    python "$SCRIPT_DIR/python/validate_metadata.py" \
    --project-config "$PROJECT_CONFIG" \
    --samples "$ACTIVE_SAMPLES_FILE" \
    --runs "$RUNS_FILE"

echo
log_info "Metadata validation passed."

if [[ "$VALIDATE_ONLY" == true ]]; then
    log_info "Validation-only mode completed."
    exit 0
fi

echo
log_info "Preparing role-aware execution plan."

conda run -n "$PYTHON_ENVIRONMENT" \
    python "$SCRIPT_DIR/python/route_samples.py" \
    --samples "$ACTIVE_SAMPLES_FILE" \
    --project-config "$PROJECT_CONFIG" \
    --format text

ROUTING_FILE="$PROJECT_ROOT/work/tmp/routing_${RUN_TIMESTAMP}.sh"

conda run -n "$PYTHON_ENVIRONMENT" \
    python "$SCRIPT_DIR/python/route_samples.py" \
    --samples "$ACTIVE_SAMPLES_FILE" \
    --project-config "$PROJECT_CONFIG" \
    --format shell \
    > "$ROUTING_FILE"

if [[ ! -s "$ROUTING_FILE" ]]; then
    die "Routing decision file was not created: $ROUTING_FILE"
fi

# shellcheck disable=SC1090
source "$ROUTING_FILE"

if [[ "$OFF_TARGET_ANALYSIS_REQUESTED" == true ]]; then
    RUN_OFF_TARGET_ANALYSIS=true
    log_info "CLI override: off-target analysis enabled."
fi



if [[ "$RUN_OFF_TARGET_ANALYSIS" == true ]]; then
    log_info "Final module decision: off_target_analysis=RUN"
else
    log_info "Final module decision: off_target_analysis=SKIP"
fi


log_info "Routing decisions loaded from: $ROUTING_FILE"

if [[ "$PLAN_ONLY" == true ]]; then
    echo
    log_info "Plan-only mode completed."
    exit 0
fi

if [[ "$RUN_ALIGNMENT_QC" == true ]]; then
    echo
    log_info "Running alignment QC."

ALIGNMENT_QC_DIR="$PROJECT_ROOT/results/shared_qc/alignment_qc"
ALIGNMENT_QC_OUTPUT="$ALIGNMENT_QC_DIR/alignment_qc_per_sample.tsv"

create_directory "$ALIGNMENT_QC_DIR"

conda run -n "$PYTHON_ENVIRONMENT" \
    python "$SCRIPT_DIR/python/alignment_qc.py" \
    --samples "$ACTIVE_SAMPLES_FILE" \
    --filtering-config "$PROJECT_ROOT/config/filtering.yaml" \
    --output "$ALIGNMENT_QC_OUTPUT"

echo
log_info "Alignment QC completed."
log_info "Output: $ALIGNMENT_QC_OUTPUT"
else
    log_info "Skipping alignment QC."
fi

if [[ "$RUN_TARGET_METRICS" == true ]]; then
    echo
    log_info "Running target metrics."

TARGET_METRICS_DIR="$PROJECT_ROOT/results/shared_qc/target_metrics"
TARGET_METRICS_OUTPUT="$TARGET_METRICS_DIR/target_metrics_per_sample.tsv"
TARGET_METRICS_ROLE_RUN_OUTPUT="$TARGET_METRICS_DIR/target_metrics_by_selected_role_run.tsv"

create_directory "$TARGET_METRICS_DIR"

conda run -n "$PYTHON_ENVIRONMENT" \
    python "$SCRIPT_DIR/python/target_metrics.py" \
    --project-config "$PROJECT_CONFIG" \
    --samples "$ACTIVE_SAMPLES_FILE" \
    --filtering-config "$PROJECT_ROOT/config/filtering.yaml" \
    --output "$TARGET_METRICS_OUTPUT" \
    --run-output "$TARGET_METRICS_ROLE_RUN_OUTPUT"

echo
log_info "Target metrics completed."
log_info "Per-sample output: $TARGET_METRICS_OUTPUT"
log_info "Selected-role run summary: $TARGET_METRICS_ROLE_RUN_OUTPUT"
else
    log_info "Skipping target metrics."
fi

if [[ "$RUN_COVERAGE_ANALYSIS" == true ]]; then
    echo
    log_info "Running target coverage analysis."

COVERAGE_DIR="$PROJECT_ROOT/results/shared_qc/coverage"
REGION_COVERAGE_OUTPUT="$COVERAGE_DIR/coverage_per_region.tsv"
PANEL_COVERAGE_OUTPUT="$COVERAGE_DIR/coverage_panel_summary.tsv"

create_directory "$COVERAGE_DIR"

conda run -n "$PYTHON_ENVIRONMENT" \
    python "$SCRIPT_DIR/python/coverage_analysis.py" \
    --project-config "$PROJECT_CONFIG" \
    --samples "$ACTIVE_SAMPLES_FILE" \
    --filtering-config "$PROJECT_ROOT/config/filtering.yaml" \
    --region-output "$REGION_COVERAGE_OUTPUT" \
    --panel-output "$PANEL_COVERAGE_OUTPUT"

echo
log_info "Target coverage analysis completed."
log_info "Region output: $REGION_COVERAGE_OUTPUT"
log_info "Panel output: $PANEL_COVERAGE_OUTPUT"
else
    log_info "Skipping target coverage analysis."
fi

if [[ "$RUN_COVERAGE_SUMMARIES" == true ]]; then
    echo
    log_info "Preparing gene-to-target annotation."

ANNOTATION_DIR="$PROJECT_ROOT/work/intermediate/annotation"
GENE_TARGET_ANNOTATION="$ANNOTATION_DIR/gene_target_annotation.tsv"

create_directory "$ANNOTATION_DIR"

ANNOTATION_WORKBOOK=$(
    conda run -n "$PYTHON_ENVIRONMENT" \
        python -c '
import sys
import yaml

with open(sys.argv[1], encoding="utf-8") as handle:
    config = yaml.safe_load(handle)

print(config["resources"]["annotation_xlsx"])
' "$PROJECT_CONFIG"
)

ANNOTATION_SHEET=$(
    conda run -n "$PYTHON_ENVIRONMENT" \
        python -c '
import sys
import yaml

with open(sys.argv[1], encoding="utf-8") as handle:
    config = yaml.safe_load(handle)

print(config["resources"]["annotation_sheet"])
' "$PROJECT_CONFIG"
)

conda run -n "$PYTHON_ENVIRONMENT" \
    python "$SCRIPT_DIR/python/prepare_target_annotation.py" \
    --workbook "$ANNOTATION_WORKBOOK" \
    --sheet "$ANNOTATION_SHEET" \
    --output "$GENE_TARGET_ANNOTATION"

echo
log_info "Gene-to-target annotation prepared."
log_info "Output: $GENE_TARGET_ANNOTATION"

echo
log_info "Creating gene-level coverage summaries."

GENE_COVERAGE_PER_SAMPLE="$COVERAGE_DIR/coverage_by_gene_per_sample.tsv"
GENE_COVERAGE_SUMMARY="$COVERAGE_DIR/coverage_by_gene_summary.tsv"

conda run -n "$PYTHON_ENVIRONMENT" \
    python "$SCRIPT_DIR/python/coverage_summaries.py" \
    --region-coverage "$REGION_COVERAGE_OUTPUT" \
    --annotation "$GENE_TARGET_ANNOTATION" \
    --per-sample-output "$GENE_COVERAGE_PER_SAMPLE" \
    --cohort-output "$GENE_COVERAGE_SUMMARY"

echo
log_info "Gene-level coverage summaries completed."
log_info "Per-sample output: $GENE_COVERAGE_PER_SAMPLE"
log_info "Cohort output: $GENE_COVERAGE_SUMMARY"
else
    log_info "Skipping gene-level coverage summaries."
fi

if [[ "$RUN_METHOD_DEVELOPMENT_SUMMARY" == true ]]; then
    echo
    log_info "Creating method-development technical summaries."

    METHOD_DEVELOPMENT_DIR="$PROJECT_ROOT/results/method_development"
    METHOD_LIBRARY_SUMMARY="$METHOD_DEVELOPMENT_DIR/method_development_library_summary.tsv"
    METHOD_GROUP_SUMMARY="$METHOD_DEVELOPMENT_DIR/method_development_group_summary.tsv"

    create_directory "$METHOD_DEVELOPMENT_DIR"

    conda run -n "$PYTHON_ENVIRONMENT" \
        python "$SCRIPT_DIR/python/summarize_method_development.py" \
        --samples "$ACTIVE_SAMPLES_FILE" \
        --alignment-qc \
        "$PROJECT_ROOT/results/shared_qc/alignment_qc/alignment_qc_per_sample.tsv" \
        --target-metrics \
        "$PROJECT_ROOT/results/shared_qc/target_metrics/target_metrics_per_sample.tsv" \
        --panel-coverage \
        "$PROJECT_ROOT/results/shared_qc/coverage/coverage_panel_summary.tsv" \
        --library-output "$METHOD_LIBRARY_SUMMARY" \
        --group-output "$METHOD_GROUP_SUMMARY"

    echo
    log_info "Method-development summaries completed."
    log_info "Library output: $METHOD_LIBRARY_SUMMARY"
    log_info "Group output: $METHOD_GROUP_SUMMARY"
else
    log_info "Skipping method-development summaries."
fi

if [[ "$RUN_METHOD_DEVELOPMENT_REPRODUCIBILITY" == true ]]; then
    echo
    log_info "Running method-development coverage reproducibility analysis."

    METHOD_REPRODUCIBILITY_DIR="$PROJECT_ROOT/results/method_development/reproducibility"
    GENE_REPRODUCIBILITY_COMPARISON="$METHOD_REPRODUCIBILITY_DIR/gene_relative_coverage_experiment_comparison.tsv"

    create_directory "$METHOD_REPRODUCIBILITY_DIR"

    conda run -n "$PYTHON_ENVIRONMENT" \
        python "$SCRIPT_DIR/python/summarize_method_coverage_reproducibility.py" \
        --samples "$ACTIVE_SAMPLES_FILE" \
        --gene-coverage \
        "$PROJECT_ROOT/results/shared_qc/coverage/coverage_by_gene_per_sample.tsv" \
        --panel-coverage \
        "$PROJECT_ROOT/results/shared_qc/coverage/coverage_panel_summary.tsv" \
        --output-dir "$METHOD_REPRODUCIBILITY_DIR"

    conda run -n "$PYTHON_ENVIRONMENT" \
        python "$SCRIPT_DIR/python/plot_basewise_coverage_reproducibility.py" \
        --samples "$ACTIVE_SAMPLES_FILE" \
        --gene-coverage \
        "$PROJECT_ROOT/results/shared_qc/coverage/coverage_by_gene_per_sample.tsv" \
        --gene-comparison "$GENE_REPRODUCIBILITY_COMPARISON" \
        --filtering-config "$PROJECT_ROOT/config/filtering.yaml" \
        --output-dir "$METHOD_REPRODUCIBILITY_DIR" \
        --top-n 3

    echo
    log_info "Method-development coverage reproducibility completed."
    log_info "Output directory: $METHOD_REPRODUCIBILITY_DIR"
else
    log_info "Skipping method-development coverage reproducibility."
fi

if [[ "$RUN_GIAB_VALIDATION" == true ]]; then
    echo
    log_info "Running GIAB small-variant validation."

    GIAB_INTERMEDIATE_DIR="$PROJECT_ROOT/work/intermediate/giab_validation"
    GIAB_PREPARED_DIR="$GIAB_INTERMEDIATE_DIR/prepared_inputs"

    GIAB_RESULTS_DIR="$PROJECT_ROOT/results/validation_control/giab_benchmark"

    GIAB_INPUT_MANIFEST="$GIAB_PREPARED_DIR/giab_input_manifest.tsv"

    GIAB_RUN_TABLE="$GIAB_RESULTS_DIR/giab_benchmark_runs.tsv"

    GIAB_SUMMARY="$GIAB_RESULTS_DIR/giab_small_variant_benchmark_summary.tsv"

    create_directory "$GIAB_INTERMEDIATE_DIR"
    create_directory "$GIAB_RESULTS_DIR"

    TARGET_BED=$(
        conda run -n "$PYTHON_ENVIRONMENT" \
            python "$SCRIPT_DIR/python/get_project_resource.py" \
            --project-config "$PROJECT_CONFIG" \
            --key target_bed
    )

    REFERENCE_FASTA=$(
        conda run -n "$PYTHON_ENVIRONMENT" \
            python "$SCRIPT_DIR/python/get_project_resource.py" \
            --project-config "$PROJECT_CONFIG" \
            --key reference_fasta
    )

    RTG_SDF=$(
        conda run -n "$PYTHON_ENVIRONMENT" \
            python "$SCRIPT_DIR/python/get_project_resource.py" \
            --project-config "$PROJECT_CONFIG" \
            --key rtg_sdf
    )

    log_info "Preparing panel-restricted GIAB inputs."

    rm -rf "$GIAB_PREPARED_DIR"

    conda run -n giab-benchmark \
        bash "$SCRIPT_DIR/run_prepare_giab_inputs.sh" \
        --target-bed "$TARGET_BED" \
        --samples "$ACTIVE_SAMPLES_FILE" \
        --output-root "$GIAB_PREPARED_DIR"

    log_info "GIAB input preparation completed."
    log_info "Input manifest: $GIAB_INPUT_MANIFEST"

    log_info "Running hap.py with the RTG vcfeval engine."

    rm -rf "$GIAB_RESULTS_DIR"
    create_directory "$GIAB_RESULTS_DIR"

    conda run -n giab-benchmark \
        bash "$SCRIPT_DIR/run_giab_benchmark.sh" \
        --input-manifest "$GIAB_INPUT_MANIFEST" \
        --reference-fasta "$REFERENCE_FASTA" \
        --rtg-sdf "$RTG_SDF" \
        --output-root "$GIAB_RESULTS_DIR"

    log_info "hap.py benchmarking completed."
    log_info "Benchmark run table: $GIAB_RUN_TABLE"

    conda run -n "$PYTHON_ENVIRONMENT" \
        python "$SCRIPT_DIR/python/summarize_giab_benchmark.py" \
        --run-table "$GIAB_RUN_TABLE" \
        --output "$GIAB_SUMMARY"

    log_info "GIAB benchmark summary completed."
    log_info "Output: $GIAB_SUMMARY"
else
    log_info "Skipping GIAB validation."
fi

if [[ "$RUN_VARIANT_APPLICATION" == true ]]; then
    log_info "Running workflow-application variant analysis."

    VARIANT_APPLICATION_DIR="$PROJECT_ROOT/results/workflow_application"

    create_directory "$VARIANT_APPLICATION_DIR"

    conda run -n "$PYTHON_ENVIRONMENT" \
        python "$SCRIPT_DIR/python/run_variant_application.py" \
        --samples "$ACTIVE_SAMPLES_FILE" \
        --scripts-dir "$SCRIPT_DIR/python" \
        --output-root "$VARIANT_APPLICATION_DIR"

    log_info "Workflow-application variant analysis completed."
    log_info "Output: $VARIANT_APPLICATION_DIR"
else
    log_info "Skipping workflow-application variant analysis."
fi

if [[ "$RUN_IGV_REVIEW_BATCH_PREPARATION" == true ]]; then
    log_info "Preparing IGV review-candidate batch files."

    IGV_REVIEW_OUTPUT_ROOT="$PROJECT_ROOT/results/igv_review_candidates"

    IGV_REVIEW_REFERENCE_FASTA=$(
        conda run -n "$PYTHON_ENVIRONMENT" \
            python "$SCRIPT_DIR/python/get_project_resource.py" \
            --project-config "$PROJECT_CONFIG" \
            --key reference_fasta
    )

    create_directory "$IGV_REVIEW_OUTPUT_ROOT"

    conda run -n "$PYTHON_ENVIRONMENT" \
        python "$SCRIPT_DIR/python/prepare_igv_review_batches.py" \
        --samples "$ACTIVE_SAMPLES_FILE" \
        --project-root "$PROJECT_ROOT" \
        --reference-fasta "$IGV_REVIEW_REFERENCE_FASTA" \
        --analysis-role workflow_application \
        --output-root "$IGV_REVIEW_OUTPUT_ROOT" \
        --window-bp 100

    log_info "IGV review-candidate batch preparation completed."
    log_info "Output: $IGV_REVIEW_OUTPUT_ROOT"
else
    log_info "Skipping IGV review-candidate batch preparation."
fi



if [[ "$RUN_SV_PHASING_SUMMARY" == true ]]; then
    log_info "Running exploratory SV and phasing summary."

    SV_PHASING_SAMPLE_SHEET="$PROJECT_ROOT/config/samples_sv_phasing.tsv"
    SV_PHASING_RESULTS_DIR="$PROJECT_ROOT/results/exploratory/sv_phasing"
    SV_PHASING_OUTPUT="$SV_PHASING_RESULTS_DIR/sv_phasing_summary.tsv"

    if [[ ! -f "$SV_PHASING_SAMPLE_SHEET" ]]; then
        log_error "SV/phasing sample sheet not found: $SV_PHASING_SAMPLE_SHEET"
        exit 1
    fi

    SV_PHASING_REFERENCE_FASTA=$(
        conda run -n "$PYTHON_ENVIRONMENT" \
            python "$SCRIPT_DIR/python/get_project_resource.py" \
            --project-config "$PROJECT_CONFIG" \
            --key reference_fasta
    )

    create_directory "$SV_PHASING_RESULTS_DIR"

    conda run -n "$PYTHON_ENVIRONMENT" \
        bash "$SCRIPT_DIR/optional/summarize_sv_phasing.sh" \
        "$SV_PHASING_SAMPLE_SHEET" \
        "$SV_PHASING_REFERENCE_FASTA" \
        "$SV_PHASING_OUTPUT"

    log_info "Exploratory SV and phasing summary completed."
    log_info "Output: $SV_PHASING_OUTPUT"
else
    log_info "Skipping exploratory SV and phasing summary."
fi

if [[ "$RUN_OFF_TARGET_ANALYSIS" == true ]]; then
    log_info "Running exploratory base-level off-target analysis."

    TARGET_BED=$(
        conda run -n "$PYTHON_ENVIRONMENT" \
            python "$SCRIPT_DIR/python/get_project_resource.py" \
            --project-config "$PROJECT_CONFIG" \
            --key target_bed
    )

    REFERENCE_FAI=$(
        conda run -n "$PYTHON_ENVIRONMENT" \
            python "$SCRIPT_DIR/python/get_project_resource.py" \
            --project-config "$PROJECT_CONFIG" \
            --key reference_fai
    )

    OFF_TARGET_WORK_DIR="$PROJECT_ROOT/work/intermediate/off_target_analysis"
    OFF_TARGET_RESULTS_DIR="$PROJECT_ROOT/results/off_target_analysis"

    rm -rf \
        "$OFF_TARGET_WORK_DIR" \
        "$OFF_TARGET_RESULTS_DIR"

    create_directory "$OFF_TARGET_WORK_DIR"
    create_directory "$OFF_TARGET_RESULTS_DIR"

    conda run -n "$PYTHON_ENVIRONMENT" \
        python "$SCRIPT_DIR/python/off_target_analysis.py" \
        --samples "$ACTIVE_SAMPLES_FILE" \
        --target-bed "$TARGET_BED" \
        --reference-fai "$REFERENCE_FAI" \
        --work-root "$OFF_TARGET_WORK_DIR" \
        --output-root "$OFF_TARGET_RESULTS_DIR" \
        --window-size 1000000 \
        --top-windows 50

    log_info "Exploratory off-target analysis completed."
    log_info "Output: $OFF_TARGET_RESULTS_DIR"
else
    log_info "Skipping exploratory off-target analysis."
fi

if [[ "$RUN_IGV_BATCH_PREPARATION" == true ]]; then
    log_info "Preparing IGV target-region batch files."

    IGV_TARGET_BED=$(
        conda run -n "$PYTHON_ENVIRONMENT" \
            python "$SCRIPT_DIR/python/get_project_resource.py" \
            --project-config "$PROJECT_CONFIG" \
            --key target_bed
    )

    IGV_REFERENCE_FASTA=$(
        conda run -n "$PYTHON_ENVIRONMENT" \
            python "$SCRIPT_DIR/python/get_project_resource.py" \
            --project-config "$PROJECT_CONFIG" \
            --key reference_fasta
    )

    IGV_GROUPS_CONFIG="$PROJECT_ROOT/config/igv_groups.tsv"
    IGV_OUTPUT_ROOT="$PROJECT_ROOT/results/igv_target_inspection"
    GENE_TARGET_ANNOTATION="$PROJECT_ROOT/work/intermediate/annotation/gene_target_annotation.tsv"

    if [[ ! -f "$IGV_GROUPS_CONFIG" ]]; then
        log_error "IGV group configuration not found: $IGV_GROUPS_CONFIG"
        exit 1
    fi

    create_directory "$IGV_OUTPUT_ROOT"

    mapfile -t IGV_ANALYSIS_ROLES < <(
        tail -n +2 "$ACTIVE_SAMPLES_FILE" |
            cut -f2 |
            sort -u
    )

    for analysis_role in "${IGV_ANALYSIS_ROLES[@]}"; do
        case "$analysis_role" in
            method_development|validation_control|workflow_application)
                log_info "Preparing IGV batches for role: $analysis_role"

                ROLE_IGV_OUTPUT="$IGV_OUTPUT_ROOT/$analysis_role"
                create_directory "$ROLE_IGV_OUTPUT"

                conda run -n "$PYTHON_ENVIRONMENT" \
                    python "$SCRIPT_DIR/python/prepare_igv_target_batches.py" \
                    --samples "$ACTIVE_SAMPLES_FILE" \
                    --groups "$IGV_GROUPS_CONFIG" \
                    --target-bed "$IGV_TARGET_BED" \
                    --reference-fasta "$IGV_REFERENCE_FASTA" \
                    --analysis-role "$analysis_role" \
                    --output-root "$ROLE_IGV_OUTPUT" \
                    --gene-annotation "$GENE_TARGET_ANNOTATION"
                ;;
            *)
                log_warning "IGV batch preparation skipped for unsupported role: $analysis_role"
                ;;
        esac
    done

    log_info "IGV target-region batch preparation completed."
    log_info "Output: $IGV_OUTPUT_ROOT"
else
    log_info "Skipping IGV batch preparation."
fi

log_info "Current pipeline execution completed: $(timestamp)"
