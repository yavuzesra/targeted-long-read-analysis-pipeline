from __future__ import annotations

import shutil
import subprocess
from pathlib import Path

import pytest


PROJECT_ROOT = Path(__file__).resolve().parents[1]
RUN_ALL = PROJECT_ROOT / "scripts" / "run_all.sh"


def run_pipeline_plan(*extra_args: str) -> subprocess.CompletedProcess[str]:
    """Run the workflow-application execution plan without analysing BAM files."""
    return subprocess.run(
        [
            "bash",
            str(RUN_ALL),
            "--role",
            "workflow_application",
            *extra_args,
            "--plan-only",
        ],
        cwd=PROJECT_ROOT,
        capture_output=True,
        text=True,
        check=False,
    )


@pytest.mark.skipif(shutil.which("conda") is None, reason="requires conda")
def test_off_target_option_is_accepted() -> None:
    result = run_pipeline_plan("--off-target-analysis")

    combined_output = result.stdout + result.stderr

    assert result.returncode == 0, combined_output
    assert "Unknown argument: --off-target-analysis" not in combined_output
    assert "CLI override: off-target analysis enabled." in combined_output
    assert (
        "Final module decision: off_target_analysis=RUN"
        in combined_output
    )


@pytest.mark.skipif(shutil.which("conda") is None, reason="requires conda")
def test_off_target_is_not_forced_without_cli_option() -> None:
    result = run_pipeline_plan()

    combined_output = result.stdout + result.stderr

    assert result.returncode == 0, combined_output
    assert "CLI override: off-target analysis enabled." not in combined_output
    assert (
        "Final module decision: off_target_analysis=SKIP"
        in combined_output
    )
