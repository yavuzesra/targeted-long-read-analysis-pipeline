import shutil
import subprocess
from pathlib import Path

import pytest


REPO_ROOT = Path(__file__).resolve().parents[1]
RUN_ALL = REPO_ROOT / "scripts" / "run_all.sh"

CONDA_ENV = Path("/usr/share/miniconda/envs/human-variation-workflow")

HAS_PIPELINE_ENV = (
    shutil.which("conda") is not None
    and CONDA_ENV.is_dir()
)

requires_pipeline_env = pytest.mark.skipif(
    not HAS_PIPELINE_ENV,
    reason=f"requires conda environment at {CONDA_ENV}",
)


def run_pipeline_plan(*extra_args: str) -> subprocess.CompletedProcess[str]:
    return subprocess.run(
        [
            "bash",
            str(RUN_ALL),
            "--plan-only",
            *extra_args,
        ],
        cwd=REPO_ROOT,
        text=True,
        capture_output=True,
        check=False,
    )


@requires_pipeline_env
def test_off_target_option_is_accepted() -> None:
    result = run_pipeline_plan("--off-target-analysis")

    combined_output = result.stdout + result.stderr

    assert result.returncode == 0, combined_output
    assert "Unknown argument: --off-target-analysis" not in combined_output
    assert "CLI override: off-target analysis enabled." in combined_output
    assert "Final module decision: off_target_analysis=RUN" in combined_output


@requires_pipeline_env
def test_off_target_is_not_forced_without_cli_option() -> None:
    result = run_pipeline_plan()

    combined_output = result.stdout + result.stderr

    assert result.returncode == 0, combined_output
    assert "CLI override: off-target analysis enabled." not in combined_output
    assert "Final module decision: off_target_analysis=SKIP" in combined_output
