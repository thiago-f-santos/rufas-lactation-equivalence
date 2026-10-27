import json
from pathlib import Path
import sys
from unittest.mock import patch
import pytest

from rufas_lactation_equivalence.atomic_io import atomic_write_json, atomic_write_text
from rufas_lactation_equivalence.cli import main, run_pipeline

PROJECT_ROOT = Path(__file__).resolve().parent.parent


def test_atomic_write_text_and_json(tmp_path: Path) -> None:
    """Verify atomic write creates and replaces files without leftover temp files."""
    txt_file = tmp_path / "sample.txt"
    atomic_write_text(txt_file, "initial")
    assert txt_file.read_text(encoding="utf-8") == "initial"
    atomic_write_text(txt_file, "updated")
    assert txt_file.read_text(encoding="utf-8") == "updated"

    json_file = tmp_path / "sample.json"
    atomic_write_json(json_file, {"key": "value", "num": 1})
    assert json.loads(json_file.read_text(encoding="utf-8")) == {"key": "value", "num": 1}
    assert list(tmp_path.glob("*.tmp*")) == []


def test_run_pipeline_generates_reports_and_patches_rufas_json(tmp_path: Path) -> None:
    """Verify end-to-end atomic pipeline generates audit reports and exports brazil_state_to_region_mapping."""
    mock_rufas_json = tmp_path / "lactation_curve_adjustment_inputs.json"
    initial_rufas_payload = {
        "adjustments": {"region": {"Southeast": {"l": -2.0, "m": 2.59, "n": 2.6}}},
        "state_to_region_mapping": {"31": "Northern Plains", "55": "Wisconsin"},
    }
    mock_rufas_json.write_text(json.dumps(initial_rufas_payload, indent=2), encoding="utf-8")

    out_dir = tmp_path / "output"
    report = run_pipeline(
        data_dir=PROJECT_ROOT / "data",
        output_dir=out_dir,
        export_rufas_json=mock_rufas_json,
    )

    assert (out_dir / "equivalence_audit_report.json").is_file()
    assert (out_dir / "equivalence_audit_report.md").is_file()
    assert len(report["brazil_state_to_region_mapping"]) == 27

    updated_rufas = json.loads(mock_rufas_json.read_text(encoding="utf-8"))
    # US FIPS mapping must be untouched
    assert updated_rufas["state_to_region_mapping"] == {"31": "Northern Plains", "55": "Wisconsin"}
    # Brazil IBGE mapping must be injected with all 27 UFs
    assert "brazil_state_to_region_mapping" in updated_rufas
    assert len(updated_rufas["brazil_state_to_region_mapping"]) == 27
    assert updated_rufas["brazil_state_to_region_mapping"]["31"] != "Northern Plains"


def test_atomic_write_cleans_up_on_failure(tmp_path: Path) -> None:
    """Verify that failure during atomic write deletes the temporary file."""
    target_file = tmp_path / "failed.txt"
    with patch("os.fsync", side_effect=OSError("Disk write failed")):
        with pytest.raises(OSError, match="Disk write failed"):
            atomic_write_text(target_file, "content")
    assert not target_file.exists()
    assert list(tmp_path.glob("*.tmp*")) == []


def test_cli_main_execution(tmp_path: Path) -> None:
    """Verify CLI main entry point runs with arguments."""
    out_dir = tmp_path / "cli_out"
    test_args = [
        "cli.py",
        "--data-dir",
        str(PROJECT_ROOT / "data"),
        "--output-dir",
        str(out_dir),
    ]
    with patch.object(sys, "argv", test_args):
        main()

    assert (out_dir / "equivalence_audit_report.json").is_file()
    assert (out_dir / "equivalence_audit_report.md").is_file()

