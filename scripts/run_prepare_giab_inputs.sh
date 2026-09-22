#!/usr/bin/env bash

set -euo pipefail

# Prepare panel-restricted truth and query VCF files for GIAB benchmarking.
#
# Expected inputs:
#   --target-bed
#       Custom panel BED file using 0-based, half-open coordinates.
#
#   --samples
#       Active samples.tsv containing one or more included
#       validation_control samples.
#
#   --output-root
#       Directory in which prepared GIAB inputs will be written.
#
# Main outputs:
#   <output-root>/shared/target_panel.sorted.bed
#   <output-root>/shared/target_panel.sorted.merged.bed
#   <output-root>/<sample_id>/benchmarkable_panel.bed
#   <output-root>/<sample_id>/<sample_id>.truth.panel.vcf.gz
#   <output-root>/<sample_id>/<sample_id>.query.panel.vcf.gz
#   <output-root>/<sample_id>/preparation_summary.tsv
#   <output-root>/giab_prepared_inputs.tsv
#
# Failure behaviour:
#   The script stops immediately if an input is missing, a required tool
#   is unavailable, no validation_control sample is found, or the panel
#   does not overlap a sample's GIAB confident-region BED file.

usage() {
    cat <<'EOF'
Usage:
  run_prepare_giab_inputs.sh \
    --target-bed FILE \
    --samples FILE \
    --output-root DIRECTORY
EOF
}


TARGET_BED=""
SAMPLES_FILE=""
OUTPUT_ROOT=""


while [[ $# -gt 0 ]]; do
    case "$1" in
        --target-bed)
            [[ $# -ge 2 ]] || {
                echo "ERROR: --target-bed requires a path." >&2
                exit 2
            }

            TARGET_BED="$2"
            shift 2
            ;;

        --samples)
            [[ $# -ge 2 ]] || {
                echo "ERROR: --samples requires a path." >&2
                exit 2
            }

            SAMPLES_FILE="$2"
            shift 2
            ;;

        --output-root)
            [[ $# -ge 2 ]] || {
                echo "ERROR: --output-root requires a path." >&2
                exit 2
            }

            OUTPUT_ROOT="$2"
            shift 2
            ;;

        -h|--help)
            usage
            exit 0
            ;;

        *)
            echo "ERROR: Unknown argument: $1" >&2
            usage >&2
            exit 2
            ;;
    esac
done


if [[ -z "$TARGET_BED" ]]; then
    echo "ERROR: --target-bed is required." >&2
    exit 2
fi

if [[ -z "$SAMPLES_FILE" ]]; then
    echo "ERROR: --samples is required." >&2
    exit 2
fi

if [[ -z "$OUTPUT_ROOT" ]]; then
    echo "ERROR: --output-root is required." >&2
    exit 2
fi


for required_file in \
    "$TARGET_BED" \
    "$SAMPLES_FILE"
do
    if [[ ! -s "$required_file" ]]; then
        echo \
            "ERROR: Required input is missing or empty: $required_file" \
            >&2
        exit 1
    fi
done


for tool in \
    bedtools \
    bcftools
do
    if ! command -v "$tool" >/dev/null 2>&1; then
        echo "ERROR: Required tool is unavailable: $tool" >&2
        exit 1
    fi
done


mkdir -p "$OUTPUT_ROOT"


if [[ ! -s "$TARGET_BED" ]]; then
    echo "ERROR: Target BED is missing or empty: $TARGET_BED" >&2
    exit 1
fi


SHARED_DIR="$OUTPUT_ROOT/shared"

SORTED_PANEL="$SHARED_DIR/target_panel.sorted.bed"

MERGED_PANEL="$SHARED_DIR/target_panel.sorted.merged.bed"

MANIFEST="$OUTPUT_ROOT/giab_input_manifest.tsv"


mkdir -p "$SHARED_DIR"


bedtools sort \
    -i "$TARGET_BED" \
    > "$SORTED_PANEL"


bedtools merge \
    -i "$SORTED_PANEL" \
    > "$MERGED_PANEL"


if [[ ! -s "$MERGED_PANEL" ]]; then
    echo \
        "ERROR: Sorted and merged target BED was not created." \
        >&2
    exit 1
fi


printf \
    'sample_id\ttruth_version\tsex\tbenchmark_bed\ttruth_panel_vcf\tquery_panel_vcf\tbenchmarkable_bases\ttruth_variant_records\tquery_variant_records\n' \
    > "$MANIFEST"


VALIDATION_COUNT=0


while IFS=$'\t' read -r \
    sample_id \
    analysis_role \
    run_id \
    bam_path \
    bai_path \
    small_variant_vcf \
    small_variant_vcf_index \
    clinvar_vcf \
    clinvar_vcf_index \
    truth_vcf \
    truth_vcf_index \
    truth_bed \
    truth_version \
    sex \
    include \
    notes
do
    if [[ "$sample_id" == "sample_id" ]]; then
        continue
    fi

    if [[ "${include,,}" != "true" ]]; then
        continue
    fi

    if [[ "$analysis_role" != "validation_control" ]]; then
        continue
    fi

    VALIDATION_COUNT=$((VALIDATION_COUNT + 1))

    echo "[INFO] Preparing GIAB inputs: $sample_id"


    for required_file in \
        "$small_variant_vcf" \
        "$small_variant_vcf_index" \
        "$truth_vcf" \
        "$truth_vcf_index" \
        "$truth_bed"
    do
        if [[ ! -s "$required_file" ]]; then
            echo \
                "ERROR: Missing required GIAB input for ${sample_id}: ${required_file}" \
                >&2
            exit 1
        fi
    done


    SAMPLE_DIR="$OUTPUT_ROOT/$sample_id"

    BENCHMARK_BED="$SAMPLE_DIR/benchmarkable_panel.bed"

    TRUTH_PANEL="$SAMPLE_DIR/${sample_id}.truth.panel.vcf.gz"

    QUERY_PANEL="$SAMPLE_DIR/${sample_id}.query.panel.vcf.gz"

    SUMMARY="$SAMPLE_DIR/preparation_summary.tsv"


    mkdir -p "$SAMPLE_DIR"


    bedtools intersect \
        -a "$MERGED_PANEL" \
        -b "$truth_bed" \
        > "$BENCHMARK_BED"


    if [[ ! -s "$BENCHMARK_BED" ]]; then
        echo \
            "ERROR: No overlap between the panel and truth BED for ${sample_id}." \
            >&2
        exit 1
    fi


    bcftools view \
        -R "$BENCHMARK_BED" \
        -Oz \
        -o "$TRUTH_PANEL" \
        "$truth_vcf"


    bcftools index \
        --tbi \
        --force \
        "$TRUTH_PANEL"


    bcftools view \
        -R "$BENCHMARK_BED" \
        -Oz \
        -o "$QUERY_PANEL" \
        "$small_variant_vcf"


    bcftools index \
        --tbi \
        --force \
        "$QUERY_PANEL"


    for prepared_vcf in \
        "$TRUTH_PANEL" \
        "$QUERY_PANEL"
    do
        if [[ ! -s "$prepared_vcf" ]]; then
            echo \
                "ERROR: Prepared VCF is missing or empty: $prepared_vcf" \
                >&2
            exit 1
        fi

        if [[ ! -s "${prepared_vcf}.tbi" ]]; then
            echo \
                "ERROR: Prepared VCF index is missing: ${prepared_vcf}.tbi" \
                >&2
            exit 1
        fi

        bcftools view \
            -h \
            "$prepared_vcf" \
            >/dev/null
    done


    BENCHMARKABLE_BASES=$(
        awk '
            {
                total += $3 - $2
            }

            END {
                print total + 0
            }
        ' "$BENCHMARK_BED"
    )


    TRUTH_COUNT=$(
        bcftools view \
            -H \
            "$TRUTH_PANEL" |
        wc -l
    )


    QUERY_COUNT=$(
        bcftools view \
            -H \
            "$QUERY_PANEL" |
        wc -l
    )


    {
        printf 'metric\tvalue\n'
        printf 'sample_id\t%s\n' "$sample_id"
        printf 'analysis_role\t%s\n' "$analysis_role"
        printf 'run_id\t%s\n' "$run_id"
        printf 'truth_version\t%s\n' "$truth_version"
        printf 'sex\t%s\n' "$sex"
        printf 'benchmarkable_bases\t%s\n' "$BENCHMARKABLE_BASES"
        printf 'truth_variant_records\t%s\n' "$TRUTH_COUNT"
        printf 'query_variant_records\t%s\n' "$QUERY_COUNT"
    } > "$SUMMARY"


    printf \
        '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
        "$sample_id" \
        "$truth_version" \
        "$sex" \
        "$BENCHMARK_BED" \
        "$TRUTH_PANEL" \
        "$QUERY_PANEL" \
        "$BENCHMARKABLE_BASES" \
        "$TRUTH_COUNT" \
        "$QUERY_COUNT" \
        >> "$MANIFEST"
done < "$SAMPLES_FILE"


if (( VALIDATION_COUNT == 0 )); then
    echo \
        "ERROR: No included validation_control samples were found." \
        >&2
    exit 1
fi


if [[ ! -s "$MANIFEST" ]]; then
    echo "ERROR: GIAB input manifest was not created." >&2
    exit 1
fi


echo "[INFO] Prepared GIAB controls: $VALIDATION_COUNT"
echo "[INFO] GIAB input manifest: $MANIFEST"
