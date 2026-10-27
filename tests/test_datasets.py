import json
from pathlib import Path
import pytest

from rufas_lactation_equivalence.biophysics import normalize_profile
from rufas_lactation_equivalence.models import RawEdaphoclimaticProfile
from rufas_lactation_equivalence.similarity import compute_all_equivalences

DATA_DIR = Path(__file__).resolve().parent.parent / "data"


def test_reference_coordinates_completeness() -> None:
    """Verify reference_coordinates.json exists and contains 48 US states and 27 Brazilian UFs."""
    coords_path = DATA_DIR / "reference_coordinates.json"
    assert coords_path.exists(), f"Missing {coords_path}"

    with open(coords_path, "r", encoding="utf-8") as f:
        data = json.load(f)

    assert "us_states" in data
    assert "brazil_ufs" in data

    us_states = data["us_states"]
    brazil_ufs = data["brazil_ufs"]

    assert len(us_states) == 48, f"Expected 48 US states, got {len(us_states)}"
    assert len(brazil_ufs) == 27, f"Expected 27 Brazilian UFs, got {len(brazil_ufs)}"

    # Check US states structure
    for s in us_states:
        assert "code" in s and s["code"].isdigit()
        assert "abbrev" in s and len(s["abbrev"]) == 2
        assert "name" in s and len(s["name"]) > 0
        assert "rufas_region" in s and s["rufas_region"] not in (None, "None")
        assert -90.0 <= s["lat"] <= 90.0
        assert -180.0 <= s["lon"] <= 180.0

    # Check Brazil UFs structure & exact IBGE codes 11..53
    expected_ibge_codes = {
        "11", "12", "13", "14", "15", "16", "17",
        "21", "22", "23", "24", "25", "26", "27", "28", "29",
        "31", "32", "33", "35",
        "41", "42", "43",
        "50", "51", "52", "53",
    }
    br_codes = {uf["code"] for uf in brazil_ufs}
    assert br_codes == expected_ibge_codes

    for uf in brazil_ufs:
        assert "code" in uf
        assert "abbrev" in uf and len(uf["abbrev"]) == 2
        assert "name" in uf and len(uf["name"]) > 0
        assert -90.0 <= uf["lat"] <= 90.0
        assert -180.0 <= uf["lon"] <= 180.0


def test_rufas_region_adjustments_match_rufas_inputs() -> None:
    """Verify rufas_region_adjustments.json contains region adjustments and state mapping."""
    adj_path = DATA_DIR / "rufas_region_adjustments.json"
    assert adj_path.exists(), f"Missing {adj_path}"

    with open(adj_path, "r", encoding="utf-8") as f:
        data = json.load(f)

    assert "region_adjustments" in data
    assert "state_to_region_mapping" in data

    region_adjustments = data["region_adjustments"]
    state_to_region = data["state_to_region_mapping"]

    expected_regions = {
        "Appalachian", "Corn Belt", "Delta", "Lake", "Mountain",
        "Northeast", "Northern Plains", "New York", "Pennsylvania",
        "Southeast", "Southern Plains", "West Coast", "Wisconsin", "None",
    }
    assert set(region_adjustments.keys()) == expected_regions

    for reg, params in region_adjustments.items():
        assert "l" in params and "m" in params and "n" in params


def test_generated_datasets_structure_and_provenance() -> None:
    """Verify us_states_edaphoclimatic.json and brazil_ufs_edaphoclimatic.json structure."""
    us_path = DATA_DIR / "us_states_edaphoclimatic.json"
    br_path = DATA_DIR / "brazil_ufs_edaphoclimatic.json"

    assert us_path.exists(), f"Missing {us_path}"
    assert br_path.exists(), f"Missing {br_path}"

    with open(us_path, "r", encoding="utf-8") as f:
        us_data = json.load(f)
    with open(br_path, "r", encoding="utf-8") as f:
        br_data = json.load(f)

    # Check provenance metadata
    for ds in (us_data, br_data):
        assert "metadata" in ds
        assert "sources" in ds["metadata"]
        meta_str = str(ds["metadata"]).lower()
        assert "nasa" in meta_str
        assert "soilgrids" in meta_str or "isric" in meta_str

    assert len(us_data["profiles"]) == 48
    assert len(br_data["profiles"]) == 27


def test_datasets_load_and_normalize_cleanly() -> None:
    """Verify all raw profiles load cleanly and normalize into valid [0.0, 1.0]^6 vectors."""
    us_path = DATA_DIR / "us_states_edaphoclimatic.json"
    br_path = DATA_DIR / "brazil_ufs_edaphoclimatic.json"

    with open(us_path, "r", encoding="utf-8") as f:
        us_data = json.load(f)
    with open(br_path, "r", encoding="utf-8") as f:
        br_data = json.load(f)

    us_profiles = [RawEdaphoclimaticProfile.from_dict(p) for p in us_data["profiles"]]
    br_profiles = [RawEdaphoclimaticProfile.from_dict(p) for p in br_data["profiles"]]

    assert len(us_profiles) == 48
    assert len(br_profiles) == 27

    for p in us_profiles + br_profiles:
        vec = normalize_profile(p)
        for val in vec.as_tuple():
            assert 0.0 <= val <= 1.0


def test_compute_all_equivalences_with_real_datasets() -> None:
    """Verify full equivalence pipeline across all 27 Brazilian UFs and 48 US states."""
    us_path = DATA_DIR / "us_states_edaphoclimatic.json"
    br_path = DATA_DIR / "brazil_ufs_edaphoclimatic.json"
    adj_path = DATA_DIR / "rufas_region_adjustments.json"

    with open(us_path, "r", encoding="utf-8") as f:
        us_profiles = [RawEdaphoclimaticProfile.from_dict(p) for p in json.load(f)["profiles"]]
    with open(br_path, "r", encoding="utf-8") as f:
        br_profiles = [RawEdaphoclimaticProfile.from_dict(p) for p in json.load(f)["profiles"]]
    with open(adj_path, "r", encoding="utf-8") as f:
        region_adjustments = json.load(f)["region_adjustments"]

    results = compute_all_equivalences(
        brazil_profiles=br_profiles,
        us_profiles=us_profiles,
        region_adjustments=region_adjustments,
        top_k=3,
    )

    assert len(results) == 27
    for res in results:
        assert res.brazil_profile.country == "BRA"
        assert res.best_match.us_code.isdigit()
        assert res.best_match.rufas_region in region_adjustments
        assert 0.0 <= res.best_match.gower_similarity <= 1.0
        assert 0.0 <= res.best_match.euclidean_similarity <= 1.0
        assert len(res.top_candidates) == 3
        assert "l" in res.wood_adjustments
        assert "m" in res.wood_adjustments
        assert "n" in res.wood_adjustments
