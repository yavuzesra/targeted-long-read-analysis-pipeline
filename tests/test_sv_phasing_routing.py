from __future__ import annotations

from pathlib import Path
import importlib.util


PROJECT_ROOT = Path(__file__).resolve().parents[1]
MODULE_PATH = PROJECT_ROOT / "scripts" / "python" / "route_samples.py"

SPEC = importlib.util.spec_from_file_location(
    "route_samples",
    MODULE_PATH,
)

assert SPEC is not None
assert SPEC.loader is not None

route_samples = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(route_samples)


def make_sample(sample_id: str, role: str) -> dict[str, str]:
    return {
        "sample_id": sample_id,
        "analysis_role": role,
        "include": "true",
        "small_variant_vcf": "NA",
    }


def test_sv_phasing_runs_for_workflow_application_when_enabled() -> None:
    samples = [
        make_sample(
            sample_id="CLINICAL_TEST",
            role="workflow_application",
        )
    ]

    project_config = {
        "modules": {
            "sv_phasing_summary": True,
        }
    }

    plan = route_samples.build_execution_plan(
        samples,
        project_config,
    )

    assert plan["modules"]["sv_phasing_summary"] is True


def test_sv_phasing_skips_without_workflow_application_samples() -> None:
    samples = [
        make_sample(
            sample_id="HG_TEST",
            role="validation_control",
        )
    ]

    project_config = {
        "modules": {
            "sv_phasing_summary": True,
        }
    }

    plan = route_samples.build_execution_plan(
        samples,
        project_config,
    )

    assert plan["modules"]["sv_phasing_summary"] is False
