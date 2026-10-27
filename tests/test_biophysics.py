import pytest

from rufas_lactation_equivalence.biophysics import (
    calculate_pawc,
    calculate_thi,
    min_max_normalize,
    normalize_profile,
)
from rufas_lactation_equivalence.models import (
    MatchCandidate,
    NormalizedVector,
    RawEdaphoclimaticProfile,
)


def test_calculate_thi_thermoneutral_and_heat_stress() -> None:
    """Verify NRC (1971) / Bohmanova et al. (2007) THI formula."""
    # At T = 20°C (68°F), RH = 50% -> THI = 68 - (0.55 - 0.275)*(68 - 58) = 68 - 2.75 = 65.25
    thi_mild = calculate_thi(temp_c=20.0, rh_pct=50.0)
    assert pytest.approx(thi_mild, rel=1e-4) == 65.25

    # At T = 30°C (86°F), RH = 80% -> THI = 86 - (0.55 - 0.44)*(86 - 58) = 86 - 3.08 = 82.92
    thi_hot = calculate_thi(temp_c=30.0, rh_pct=80.0)
    assert pytest.approx(thi_hot, rel=1e-4) == 82.92


def test_calculate_pawc_loam_vs_sand() -> None:
    """Verify Saxton & Rawls (2006) PAWC increases from coarse sand to silt-clay loam."""
    pawc_sand = calculate_pawc(sand_pct=85.0, clay_pct=5.0, soc_g_kg=8.0)
    pawc_loam = calculate_pawc(sand_pct=30.0, clay_pct=25.0, soc_g_kg=25.0)
    assert 0.03 <= pawc_sand < pawc_loam <= 0.25


def test_min_max_normalize_clamping_and_scaling() -> None:
    """Verify strict [0.0, 1.0] clamping and linear midpoint."""
    assert min_max_normalize(10.0, 0.0, 20.0) == 0.5
    assert min_max_normalize(-5.0, 0.0, 20.0) == 0.0
    assert min_max_normalize(25.0, 0.0, 20.0) == 1.0

    with pytest.raises(ValueError, match="upper .* must be greater than lower"):
        min_max_normalize(5.0, 10.0, 5.0)

    with pytest.raises(ValueError, match="upper .* must be greater than lower"):
        min_max_normalize(5.0, 10.0, 4.0)


def test_models_dataclass_methods() -> None:
    """Verify RawEdaphoclimaticProfile.from_dict and dictionary serialization."""
    data = {
        "code": "31",
        "abbrev": "MG",
        "name": "Minas Gerais",
        "country": "BRA",
        "rufas_region": None,
        "temp_mean_c": 21.2,
        "rh_mean_pct": 70.0,
        "temp_seasonality_c": 6.5,
        "precip_annual_mm": 1450.0,
        "pet_annual_mm": 1150.0,
        "soil_clay_pct": 46.0,
        "soil_sand_pct": 34.0,
        "soil_ph": 5.2,
        "soil_soc_g_kg": 19.0,
    }
    profile = RawEdaphoclimaticProfile.from_dict(data)
    assert profile.code == "31"
    assert profile.country == "BRA"

    vec = normalize_profile(profile)
    d = vec.to_dict()
    assert set(d.keys()) == {
        "thi_index",
        "temp_seasonality",
        "moisture_index",
        "soil_clay_fraction",
        "soil_pawc_index",
        "soil_fertility_index",
    }
    assert all(0.0 <= v <= 1.0 for v in d.values())

    cand = MatchCandidate(
        us_code="12",
        us_abbrev="FL",
        us_name="Florida",
        rufas_region="Southeast",
        gower_similarity=0.9123456,
        euclidean_similarity=0.8876543,
    )
    cand_d = cand.to_dict()
    assert cand_d["us_code"] == "12"
    assert cand_d["gower_similarity"] == 0.912346


def test_normalize_profile_produces_strictly_bounded_vector() -> None:
    """Verify full profile decomposition into 6D vector in [0.0, 1.0]."""
    raw = RawEdaphoclimaticProfile(
        code="31",
        abbrev="MG",
        name="Minas Gerais",
        country="BRA",
        rufas_region=None,
        temp_mean_c=21.2,
        rh_mean_pct=70.0,
        temp_seasonality_c=6.5,
        precip_annual_mm=1450.0,
        pet_annual_mm=1150.0,
        soil_clay_pct=46.0,
        soil_sand_pct=34.0,
        soil_ph=5.2,
        soil_soc_g_kg=19.0,
    )
    vec = normalize_profile(raw)
    assert isinstance(vec, NormalizedVector)
    for val in vec.as_tuple():
        assert 0.0 <= val <= 1.0
    assert len(vec.as_tuple()) == 6


def test_normalized_vector_bounds_validation() -> None:
    """Verify NormalizedVector raises ValueError on out-of-bounds fields."""
    # Valid vector
    NormalizedVector(0.0, 0.5, 1.0, 0.2, 0.4, 0.8)

    # Test each field < 0.0 or > 1.0
    with pytest.raises(ValueError, match="thi_index"):
        NormalizedVector(-0.01, 0.5, 1.0, 0.2, 0.4, 0.8)
    with pytest.raises(ValueError, match="temp_seasonality"):
        NormalizedVector(0.5, 1.05, 1.0, 0.2, 0.4, 0.8)
    with pytest.raises(ValueError, match="moisture_index"):
        NormalizedVector(0.5, 0.5, -0.1, 0.2, 0.4, 0.8)
    with pytest.raises(ValueError, match="soil_clay_fraction"):
        NormalizedVector(0.5, 0.5, 0.5, 1.5, 0.4, 0.8)
    with pytest.raises(ValueError, match="soil_pawc_index"):
        NormalizedVector(0.5, 0.5, 0.5, 0.5, -0.5, 0.8)
    with pytest.raises(ValueError, match="soil_fertility_index"):
        NormalizedVector(0.5, 0.5, 0.5, 0.5, 0.5, 1.001)


def test_match_candidate_bounds_validation() -> None:
    """Verify MatchCandidate raises ValueError when similarity is outside [0.0, 1.0]."""
    with pytest.raises(ValueError, match="gower_similarity"):
        MatchCandidate("12", "FL", "Florida", "Southeast", -0.1, 0.5)
    with pytest.raises(ValueError, match="gower_similarity"):
        MatchCandidate("12", "FL", "Florida", "Southeast", 1.01, 0.5)
    with pytest.raises(ValueError, match="euclidean_similarity"):
        MatchCandidate("12", "FL", "Florida", "Southeast", 0.5, -0.01)
    with pytest.raises(ValueError, match="euclidean_similarity"):
        MatchCandidate("12", "FL", "Florida", "Southeast", 0.5, 1.1)


def test_raw_profile_validation() -> None:
    """Verify RawEdaphoclimaticProfile raises ValueError on invalid physical parameters."""
    valid_args = dict(
        code="31",
        abbrev="MG",
        name="Minas Gerais",
        country="BRA",
        rufas_region=None,
        temp_mean_c=21.2,
        rh_mean_pct=70.0,
        temp_seasonality_c=6.5,
        precip_annual_mm=1450.0,
        pet_annual_mm=1150.0,
        soil_clay_pct=46.0,
        soil_sand_pct=34.0,
        soil_ph=5.2,
        soil_soc_g_kg=19.0,
    )

    # Valid profile should succeed
    RawEdaphoclimaticProfile(**valid_args)

    # RH out of bounds
    with pytest.raises(ValueError, match="rh_mean_pct"):
        RawEdaphoclimaticProfile(**{**valid_args, "rh_mean_pct": -1.0})
    with pytest.raises(ValueError, match="rh_mean_pct"):
        RawEdaphoclimaticProfile(**{**valid_args, "rh_mean_pct": 100.1})

    # Clay out of bounds
    with pytest.raises(ValueError, match="soil_clay_pct"):
        RawEdaphoclimaticProfile(**{**valid_args, "soil_clay_pct": -5.0})
    with pytest.raises(ValueError, match="soil_clay_pct"):
        RawEdaphoclimaticProfile(**{**valid_args, "soil_clay_pct": 101.0})

    # Sand out of bounds
    with pytest.raises(ValueError, match="soil_sand_pct"):
        RawEdaphoclimaticProfile(**{**valid_args, "soil_sand_pct": -0.1})
    with pytest.raises(ValueError, match="soil_sand_pct"):
        RawEdaphoclimaticProfile(**{**valid_args, "soil_sand_pct": 105.0})

    # Precip negative
    with pytest.raises(ValueError, match="precip_annual_mm"):
        RawEdaphoclimaticProfile(**{**valid_args, "precip_annual_mm": -10.0})

    # PET negative
    with pytest.raises(ValueError, match="pet_annual_mm"):
        RawEdaphoclimaticProfile(**{**valid_args, "pet_annual_mm": -5.0})

    # pH <= 0
    with pytest.raises(ValueError, match="soil_ph"):
        RawEdaphoclimaticProfile(**{**valid_args, "soil_ph": 0.0})
    with pytest.raises(ValueError, match="soil_ph"):
        RawEdaphoclimaticProfile(**{**valid_args, "soil_ph": -1.0})

    # SOC < 0
    with pytest.raises(ValueError, match="soil_soc_g_kg"):
        RawEdaphoclimaticProfile(**{**valid_args, "soil_soc_g_kg": -0.5})
