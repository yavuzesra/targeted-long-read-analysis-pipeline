# Targeted Long-Read Sequencing Analysis Pipeline

[![Tests](https://github.com/yavuzesra/targeted-long-read-analysis-pipeline/actions/workflows/tests.yml/badge.svg)](https://github.com/yavuzesra/targeted-long-read-analysis-pipeline/actions/workflows/tests.yml)

Metadata-driven Bash/Python workflow for quality control, target-enrichment analysis, coverage assessment, GIAB small-variant benchmarking and technical variant summarisation of targeted Oxford Nanopore Technologies (ONT) sequencing data.

The workflow was developed during a BSc biotechnology project evaluating a custom QIAseq xHYB Long Read hybrid-capture panel for genomic regions associated with autoinflammatory disorders.

> **Research use only.** This repository documents a technical-evaluation workflow and is not a clinically validated diagnostic pipeline.

## What the pipeline does

The workflow supports three analysis roles:

- `method_development` — evaluates technical-library and run reproducibility;
- `validation_control` — benchmarks small variants against GIAB truth sets;
- `workflow_application` — generates QC and technical variant summaries for application samples.

Core modules include:

```text
BAM / VCF / metadata
        |
        v
Input and metadata validation
        |
        +-- Alignment QC
        +-- Target coverage
        +-- Target-enrichment metrics
        +-- Gene-level coverage summaries
        +-- GIAB SNP/INDEL benchmarking
        +-- Technical variant summaries
        +-- Optional IGV / off-target / SV-phasing summaries
```

## Technical evaluation results

The thesis evaluation used three GIAB reference samples and five clinical workflow-application samples. Aggregate results reported in the thesis included:

| Metric | Result |
|---|---:|
| Target enrichment | 684–753× |
| Target bases covered at ≥20× | 97.2–99.6% |
| GIAB SNP F1 | 99.59% |
| GIAB INDEL F1 | 82.51% |

These values describe the evaluated dataset and should not be interpreted as expected performance on new data. Patient-level sequencing files and patient-level variant results are not distributed in this repository.

## Repository structure

```text
.
├── .github/workflows/tests.yml
├── config/
│   ├── filtering.yaml
│   ├── giab_sensitivity.tsv
│   ├── igv_groups.tsv
│   ├── project.yaml
│   ├── runs.tsv
│   ├── samples.tsv
│   └── samples_sv_phasing.tsv
├── resources/
│   ├── examples/example_targets.bed
│   ├── panel/autoinflammatory.bed
│   ├── panel/panel_annotation.xlsx
│   └── reference/GRCh38_no_alt_analysis_set.fna.fai
├── scripts/
│   ├── lib/
│   ├── optional/
│   ├── python/
│   ├── run_all.sh
│   ├── run_giab_benchmark.sh
│   ├── run_giab_truthset_sensitivity.sh
│   ├── run_prepare_giab_inputs.sh
│   └── run_sequencing_run_target_metrics.sh
├── tests/
├── environment.yml
├── environment-benchmark.yml
├── .gitignore
├── pytest.ini
└── README.md
```

`results/`, `work/`, sequencing files and patient-level outputs are intentionally excluded from version control.

## Panel and reference conventions

The analysis uses the GRCh38 no-alt analysis set:

```text
GCA_000001405.15_GRCh38_no_alt_analysis_set
```

The project panel is included as:

```text
resources/panel/autoinflammatory.bed
```

The configured panel contains 28 target intervals representing 29 gene-to-target relationships, with a merged target territory of 1,203,874 bp.

A minimal synthetic BED is also supplied:

```text
resources/examples/example_targets.bed
```

It exists only for examples and lightweight software tests; it is not the biological panel used in the thesis evaluation.

BED coordinates are interpreted as 0-based, half-open. IGV displays loci using 1-based coordinates.

The complete GRCh38 FASTA and RTG SDF are deliberately external because of their size. Their local paths must be configured before running analyses.

## Installation

The main environment contains the general analysis and test dependencies:

```bash
conda env create -f environment.yml
conda activate targeted-long-read-analysis
```

GIAB benchmarking uses a separate legacy-compatible environment because the Bioconda build of hap.py 0.3.15 uses Python 2.7:

```bash
conda env create -f environment-benchmark.yml
conda activate targeted-long-read-benchmark
```

The environment files document the principal dependencies rather than an OS-specific full package export.

Main command-line dependencies include:

- Python 3
- `samtools`
- `bcftools`
- `bedtools`
- PyYAML
- openpyxl
- NumPy
- Matplotlib
- pytest
- hap.py and RTG Tools for GIAB benchmarking

The thesis analysis used `samtools` 1.18 and `bcftools` 1.17. Upstream ONT workflow versions recorded for the evaluated dataset were `wf-alignment` 1.2.5/1.2.6 and `wf-human-variation` 2.8.0.

## Configuration

### `config/project.yaml`

Set local paths for large external resources before running the workflow:

```yaml
resources:
  reference_fasta: /path/to/GRCh38_no_alt_analysis_set.fna
  reference_fai: /path/to/GRCh38_no_alt_analysis_set.fna.fai
  target_bed: resources/panel/autoinflammatory.bed
  annotation_xlsx: resources/panel/panel_annotation.xlsx
  rtg_sdf: /path/to/GRCh38_no_alt.sdf
```

The repository does not contain personal workstation paths.

### `config/samples.tsv`

The supplied rows are non-identifying examples and are disabled with `include=false`. Replace paths and metadata locally, then set selected rows to `include=true`.

Public GIAB identifiers such as `HG002` may be used directly. Clinical examples use synthetic identifiers such as `CLINICAL_SAMPLE_01`.

Role-specific requirements are validated before analysis:

- method-development samples require BAM and BAI files;
- validation controls additionally require query VCF, truth VCF, indexes, confident-region BED, truth version and sex;
- workflow-application samples require BAM/BAI and technical plus ClinVar-annotated small-variant VCF pairs.

### Other configuration files

- `config/runs.tsv` records sequencing-run provenance;
- `config/filtering.yaml` defines SAM flag masks and analysis assumptions;
- `config/igv_groups.tsv` defines optional IGV sample groups;
- `config/samples_sv_phasing.tsv` maps application samples to optional EPI2ME output directories;
- `config/giab_sensitivity.tsv` defines an optional HG002 alternative-truth-set comparison.

## Running the pipeline

Run commands from the repository root.

Inspect a planned execution:

```bash
bash scripts/run_all.sh \
    --role validation_control \
    --plan-only
```

Run one analytical role after configuring its input files:

```bash
bash scripts/run_all.sh --role method_development
bash scripts/run_all.sh --role validation_control
bash scripts/run_all.sh --role workflow_application
```

Restrict method-development processing to one run:

```bash
bash scripts/run_all.sh \
    --role method_development \
    --run-id method_run_1
```

Enable exploratory off-target analysis for a specific execution:

```bash
bash scripts/run_all.sh \
    --role workflow_application \
    --off-target-analysis
```

## Main analyses

### Alignment QC

Reports mapping and read-level summary metrics from each BAM. Secondary and supplementary records are excluded from primary-read calculations according to `config/filtering.yaml`.

Main output:

```text
results/shared_qc/alignment_qc/alignment_qc_per_sample.tsv
```

### Target coverage

Coverage is calculated for every target interval and for the merged panel. Zero-depth target bases remain in the denominators. The configured technical QC criterion is ≥90% of target bases covered at ≥20×.

Main outputs:

```text
results/shared_qc/coverage/coverage_per_region.tsv
results/shared_qc/coverage/coverage_panel_summary.tsv
results/shared_qc/coverage/coverage_by_gene_per_sample.tsv
results/shared_qc/coverage/coverage_by_gene_summary.tsv
```

### Target enrichment

The analysis flag mask excludes unmapped, secondary, QC-failed, duplicate and supplementary records.

An on-target alignment record is an eligible record overlapping at least one merged target interval by at least one reference base.

```text
on-target record percentage
    = 100 × on-target records / eligible records
```

Target-base enrichment is calculated independently from depth sums:

```text
target-base enrichment factor
    = (target depth sum / merged target territory)
      / (whole-reference depth sum / reference territory)
```

### Method development

Method-development outputs combine shared QC into per-library and grouped summaries. When two run IDs are selected, the workflow also compares gene-level relative coverage and base-wise coverage reproducibility.

Relative gene coverage is:

```text
gene mean depth / panel mean depth
```

### GIAB validation controls

The GIAB module intersects the panel BED with each control's confident-region BED, restricts truth and query VCFs to the benchmarkable region, and runs hap.py with the RTG `vcfeval` engine.

The thesis evaluation used:

- HG002: GIAB v5.0q
- HG003: GIAB v4.2.1
- HG004: GIAB v4.2.1

The summary reports truth/query totals, TP, FP, FN, unknown calls, precision, recall and F1 separately for SNPs and INDELs.

### Workflow application

Application samples have no independent truth set in this workflow. The module summarises existing small-variant and ClinVar-annotated VCFs and prepares technical review tables.

A ClinVar annotation or VCF `PASS` status is not an independent clinical classification.

## Testing

Run the Python test suite:

```bash
pytest -q
```

Two plan-level integration tests are automatically skipped when `conda` is unavailable; in a configured Ubuntu/WSL environment they run as part of the same command.

Additional syntax checks:

```bash
python -m compileall -q scripts tests

for shell_file in scripts/*.sh scripts/lib/*.sh scripts/optional/*.sh; do
    bash -n "$shell_file" || exit 1
done
```

GitHub Actions runs the portable Python tests plus Python and shell syntax checks on every push and pull request. A green `Tests` badge therefore means that these defined automated checks passed; it does not establish biological or clinical validity.

## Data-protection boundary

This public repository contains source code, non-identifying example metadata, public GIAB identifiers and the target-region definition used in the technical evaluation.

It intentionally excludes:

- FASTQ, BAM, CRAM and patient-level VCF files;
- patient-level variant tables;
- IGV screenshots derived from clinical samples;
- real clinical sample identifiers;
- EPI2ME instance identifiers;
- local workstation paths;
- generated `results/`, `work/` and log files.

## Limitations

- Performance values above apply only to the evaluated dataset and panel.
- GIAB benchmarking applies to configured small variants within benchmarkable confident regions.
- Workflow-application variants, structural variants, phasing, off-target regions and IGV observations are not independently validated by this repository.
- Low coverage, repetitive sequence, strand imbalance and mapping artefacts require cautious interpretation.

## Citation / project context

Developed as part of the BSc thesis **“Long-read-based diagnostic of autoinflammatory disorders”** at HES-SO Valais-Wallis, with clinical workflow application performed in collaboration with CHUV.

## License

No software license is assigned in this public-preparation version. Add a license only after confirming the appropriate institutional and project-level redistribution terms.
