import pytest

from rufas_lactation_equivalence.models import NormalizedVector, RawEdaphoclimaticProfile
from rufas_lactation_equivalence.similarity import (
    compute_all_equivalences,
    euclidean_similarity,
    gower_similarity,
)


def _create_sample_profile(
    code: str,
    abbrev: str,
    name: str,
    country: str,
    rufas_region: str | None,
    temp_mean_c: float = 20.0,
    rh_mean_pct: float = 70.0,
    temp_seasonality_c: float = 10.0,
    precip_annual_mm: float = 1200.0,
    pet_annual_mm: float = 1000.0,
    soil_clay_pct: float = 30.0,
    soil_sand_pct: float = 50.0,
    soil_ph: float = 6.0,
    soil_soc_g_kg: float = 20.0,
) -> RawEdaphoclimaticProfile:
    return RawEdaphoclimaticProfile(
        code=code,
        abbrev=abbrev,
        name=name,
        country=country,
        rufas_region=rufas_region,
        temp_mean_c=temp_mean_c,
        rh_mean_pct=rh_mean_pct,
        temp_seasonality_c=temp_seasonality_c,
        precip_annual_mm=precip_annual_mm,
        pet_annual_mm=pet_annual_mm,
        soil_clay_pct=soil_clay_pct,
        soil_sand_pct=soil_sand_pct,
        soil_ph=soil_ph,
        soil_soc_g_kg=soil_soc_g_kg,
    )


def test_gower_and_euclidean_metric_properties() -> None:
    """Verify identity, symmetry, and [0.0, 1.0] range of Gower and Euclidean similarities."""
    zeros = NormalizedVector(0.0, 0.0, 0.0, 0.0, 0.0, 0.0)
    ones = NormalizedVector(1.0, 1.0, 1.0, 1.0, 1.0, 1.0)
    mid = NormalizedVector(0.5, 0.5, 0.5, 0.5, 0.5, 0.5)

    assert gower_similarity(zeros, zeros) == 1.0
    assert euclidean_similarity(ones, ones) == 1.0
    assert gower_similarity(zeros, ones) == 0.0
    assert euclidean_similarity(zeros, ones) == 0.0
    assert gower_similarity(zeros, mid) == pytest.approx(0.5)
    assert euclidean_similarity(zeros, mid) == pytest.approx(0.5)
    assert gower_similarity(zeros, mid) == gower_similarity(mid, zeros)
    assert euclidean_similarity(zeros, mid) == euclidean_similarity(mid, zeros)


def test_compute_all_equivalences_selects_closest_homologue() -> None:
    """Verify a warm/humid UF matches Southeast over cold/continental Wisconsin."""
    br_uf = _create_sample_profile(
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
    us_fl = _create_sample_profile(
        code="12",
        abbrev="FL",
        name="Florida",
        country="USA",
        rufas_region="Southeast",
        temp_mean_c=21.5,
        rh_mean_pct=74.0,
        temp_seasonality_c=11.0,
        precip_annual_mm=1380.0,
        pet_annual_mm=1200.0,
        soil_clay_pct=18.0,
        soil_sand_pct=65.0,
        soil_ph=5.5,
        soil_soc_g_kg=16.0,
    )
    us_wi = _create_sample_profile(
        code="55",
        abbrev="WI",
        name="Wisconsin",
        country="USA",
        rufas_region="Wisconsin",
        temp_mean_c=6.5,
        rh_mean_pct=71.0,
        temp_seasonality_c=29.0,
        precip_annual_mm=820.0,
        pet_annual_mm=650.0,
        soil_clay_pct=22.0,
        soil_sand_pct=38.0,
        soil_ph=6.4,
        soil_soc_g_kg=32.0,
    )
    adjustments = {
        "Southeast": {"l": -2.0, "m": 2.59, "n": 2.6},
        "Wisconsin": {"l": 1.4, "m": -0.83, "n": -0.2},
    }

    results = compute_all_equivalences([br_uf], [us_fl, us_wi], adjustments, top_k=2)
    assert len(results) == 1
    res = results[0]
    assert res.best_match.us_code == "12"
    assert res.best_match.rufas_region == "Southeast"
    assert res.wood_adjustments == {"l": -2.0, "m": 2.59, "n": 2.6}
    assert res.best_match.gower_similarity > res.top_candidates[1].gower_similarity
    assert len(res.top_candidates) == 2


def test_compute_all_equivalences_validation_and_errors() -> None:
    """Verify validation triggers on unknown region or missing candidate profiles."""
    br_uf = _create_sample_profile("35", "SP", "São Paulo", "BRA", None)
    us_unknown = _create_sample_profile("99", "XX", "Unknown State", "USA", "Atlantis")
    us_none_region = _create_sample_profile("98", "YY", "Unmapped State", "USA", None)
    us_literal_none = _create_sample_profile("97", "ZZ", "None State", "USA", "None")

    valid_adjustments = {"Southeast": {"l": 1.0, "m": 1.0, "n": 1.0}}

    # Unknown region not present in dictionary
    with pytest.raises(ValueError, match="references unknown RuFaS region 'Atlantis'"):
        compute_all_equivalences([br_uf], [us_unknown], valid_adjustments)

    # Empty candidate list
    with pytest.raises(ValueError, match="No valid US candidate states provided."):
        compute_all_equivalences([br_uf], [], valid_adjustments)

    # Candidates only with None or "None" rufas_region
    with pytest.raises(ValueError, match="No valid US candidate states provided."):
        compute_all_equivalences([br_uf], [us_none_region, us_literal_none], valid_adjustments)


def test_compute_all_equivalences_sorting_and_tie_breaking() -> None:
    """Verify Brazil UFs are ordered by integer code and tie-breaking ranks deterministically."""
    br_mg = _create_sample_profile("31", "MG", "Minas Gerais", "BRA", None)
    br_ro = _create_sample_profile("11", "RO", "Rondônia", "BRA", None)

    # Two US profiles that have identical values to trigger tie-breaking by code
    us_a = _create_sample_profile("20", "AA", "State A", "USA", "Southeast")
    us_b = _create_sample_profile("10", "BB", "State B", "USA", "Southeast")

    adjustments = {"Southeast": {"l": 0.0, "m": 0.0, "n": 0.0}}

    # Passing out of order [br_mg (31), br_ro (11)]
    results = compute_all_equivalences([br_mg, br_ro], [us_a, us_b], adjustments, top_k=2)

    assert len(results) == 2
    # Result should be sorted by int(code): 11 first, then 31
    assert results[0].brazil_profile.code == "11"
    assert results[1].brazil_profile.code == "31"

    # For candidates with identical similarities, lower numeric code (10 vs 20) should come first
    assert results[0].best_match.us_code == "10"
    assert results[0].top_candidates[1].us_code == "20"


def test_compute_all_equivalences_empty_brazil_list() -> None:
    """Verify empty Brazil profile list yields empty results list."""
    us_fl = _create_sample_profile("12", "FL", "Florida", "USA", "Southeast")
    adjustments = {"Southeast": {"l": 0.0, "m": 0.0, "n": 0.0}}

    results = compute_all_equivalences([], [us_fl], adjustments)
    assert results == []
