from __future__ import annotations

import json
import subprocess
import sys
from pathlib import Path

import yaml


PROJECT_ROOT = Path(__file__).resolve().parents[1]
ROUTER = PROJECT_ROOT / "scripts" / "python" / "route_samples.py"


def write_samples(tmp_path: Path) -> Path:
    """Create minimal workflow-application metadata for routing tests."""
    samples_path = tmp_path / "samples.tsv"
    samples_path.write_text(
        "\t".join(
            [
                "sample_id",
                "analysis_role",
                "small_variant_vcf",
                "include",
            ]
        )
        + "\n"
        + "\t".join(
            [
                "CLINICAL_TEST",
                "workflow_application",
                "CLINICAL_TEST.vcf.gz",
                "true",
            ]
        )
        + "\n",
        encoding="utf-8",
    )
    return samples_path


def write_project_config(
    tmp_path: Path,
    *,
    igv_batch_preparation: bool,
) -> Path:
    """Create a project configuration with a controlled IGV setting."""
    config_path = tmp_path / "project.yaml"
    config = {
        "modules": {
            "off_target_analysis": False,
            "igv_batch_preparation": igv_batch_preparation,
        }
    }

    with config_path.open("w", encoding="utf-8") as handle:
        yaml.safe_dump(config, handle, sort_keys=False)

    return config_path


def run_router(
    samples_path: Path,
    config_path: Path,
) -> dict:
    """Run the routing script and return its JSON execution plan."""
    result = subprocess.run(
        [
            sys.executable,
            str(ROUTER),
            "--samples",
            str(samples_path),
            "--project-config",
            str(config_path),
            "--format",
            "json",
        ],
        cwd=PROJECT_ROOT,
        capture_output=True,
        text=True,
        check=False,
    )

    assert result.returncode == 0, result.stderr
    return json.loads(result.stdout)


def test_igv_batch_preparation_runs_when_enabled(
    tmp_path: Path,
) -> None:
    samples_path = write_samples(tmp_path)
    config_path = write_project_config(
        tmp_path,
        igv_batch_preparation=True,
    )

    plan = run_router(samples_path, config_path)

    assert plan["modules"]["igv_batch_preparation"] is True
    assert plan["modules"]["off_target_analysis"] is False


def test_igv_batch_preparation_skips_when_disabled(
    tmp_path: Path,
) -> None:
    samples_path = write_samples(tmp_path)
    config_path = write_project_config(
        tmp_path,
        igv_batch_preparation=False,
    )

    plan = run_router(samples_path, config_path)

    assert plan["modules"]["igv_batch_preparation"] is False
    assert plan["modules"]["off_target_analysis"] is False
