from __future__ import annotations

import importlib.util
from pathlib import Path


PROJECT_ROOT = Path(__file__).resolve().parents[1]
MODULE_PATH = (
    PROJECT_ROOT
    / "scripts/python/route_samples.py"
)

SPEC = importlib.util.spec_from_file_location(
    "route_samples",
    MODULE_PATH,
)

if SPEC is None or SPEC.loader is None:
    raise RuntimeError(
        f"Unable to load module: {MODULE_PATH}"
    )

route_samples = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(route_samples)


def make_sample(
    sample_id: str,
    role: str,
    small_variant_vcf: str = "NA",
    run_id: str = "run1",
) -> dict[str, str]:
    return {
        "sample_id": sample_id,
        "analysis_role": role,
        "run_id": run_id,
        "small_variant_vcf": small_variant_vcf,
        "include": "true",
    }


def test_method_development_routing() -> None:
    samples = [
        make_sample(
            sample_id="Method01",
            role="method_development",
        )
    ]

    plan = route_samples.build_execution_plan(samples, {})

    assert plan["modules"]["alignment_qc"] is True
    assert (
        plan["modules"]["method_development_summary"]
        is True
    )
    assert (
        plan["modules"]["method_development_reproducibility"]
        is False
    )
    assert plan["modules"]["giab_validation"] is False
    assert plan["modules"]["variant_application"] is False


def test_method_reproducibility_requires_two_runs() -> None:
    samples = [
        make_sample("E1_A_1", "method_development", run_id="experiment1"),
        make_sample("E2_A_1", "method_development", run_id="experiment2"),
    ]

    plan = route_samples.build_execution_plan(samples, {})

    assert plan["modules"]["method_development_reproducibility"] is True


def test_validation_control_routing() -> None:
    samples = [
        make_sample(
            sample_id="Control01",
            role="validation_control",
            small_variant_vcf="control.vcf.gz",
        )
    ]

    plan = route_samples.build_execution_plan(samples, {})

    assert plan["modules"]["alignment_qc"] is True
    assert plan["modules"]["giab_validation"] is True
    assert (
        plan["modules"]["method_development_summary"]
        is False
    )
    assert plan["modules"]["variant_application"] is False


def test_workflow_application_with_vcf() -> None:
    samples = [
        make_sample(
            sample_id="Clinical01",
            role="workflow_application",
            small_variant_vcf="clinical.vcf.gz",
        )
    ]

    plan = route_samples.build_execution_plan(samples, {})

    assert plan["modules"]["variant_application"] is True
    assert (
        plan["workflow_application_with_vcf"]
        == ["Clinical01"]
    )
    assert (
        plan["workflow_application_without_vcf"]
        == []
    )


def test_workflow_application_without_vcf() -> None:
    samples = [
        make_sample(
            sample_id="Clinical01",
            role="workflow_application",
            small_variant_vcf="NA",
        )
    ]

    plan = route_samples.build_execution_plan(samples, {})

    assert plan["modules"]["alignment_qc"] is True
    assert plan["modules"]["variant_application"] is False
    assert (
        plan["workflow_application_without_vcf"]
        == ["Clinical01"]
    )


def test_shell_plan_output(capsys) -> None:
    plan = {
        "modules": {
            "input_validation": True,
            "giab_validation": False,
        }
    }

    route_samples.print_shell_plan(plan)

    captured = capsys.readouterr()

    assert (
        captured.out
        == "RUN_INPUT_VALIDATION=true\n"
        "RUN_GIAB_VALIDATION=false\n"
    )
