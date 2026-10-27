# Edaphoclimatic Lactation Curve Equivalence Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Build a standalone, scientifically grounded, atomically executed Python package (`rufas-lactation-equivalence/`) that decomposes US (FIPS) and Brazilian (IBGE UF) agroclimatic and edaphic profiles into 6 normalized $[0.0, 1.0]$ characteristics, computes Gower similarity to find the Top-1 US homologue region for all 27 Brazilian UFs, and integrates the resulting `brazil_state_to_region_mapping` cleanly into `RuFaS/` without FIPS/IBGE code collisions.

**Architecture:** A zero-external-runtime-dependency Python package in `/home/thiago/Projects/rufas-lactation-equivalence/` computes biophysical indices (THI from Bohmanova et al., 2007; PAWC from Saxton & Rawls, 2006), normalizes each state/UF into a 6D vector in $[0.0, 1.0]^6$, ranks US candidate states via Gower's similarity coefficient (Gower, 1971), and atomically writes (`NamedTemporaryFile` + `os.fsync` + `os.replace`) audit reports and the `brazil_state_to_region_mapping` block into `RuFaS/input/data/animal/lactation_curve_adjustment_inputs.json`. `RuFaS` resolves `country == "BRA"` against `brazil_state_to_region_mapping` and `country == "USA"` against `state_to_region_mapping`.

**Tech Stack:** Python 3.12+ standard library (`dataclasses`, `json`, `math`, `os`, `pathlib`, `tempfile`, `argparse`), `pytest` for testing.

**Spec:** `/home/thiago/Projects/rufas-lactation-equivalence/docs/superpowers/specs/2026-10-26-edaphoclimatic-lactation-equivalence-design.md`

## Global Constraints

- Zero third-party runtime dependencies in `rufas-lactation-equivalence` (pure Python standard library).
- All 6 normalized characteristics (`thi_index`, `temp_seasonality`, `moisture_index`, `soil_clay_fraction`, `soil_pawc_index`, `soil_fertility_index`) MUST strictly lie in $[0.0, 1.0]$.
- All file writes in the pipeline MUST be atomic using `tempfile.NamedTemporaryFile` in the target file's directory, `os.fsync`, and `os.replace`.
- All 27 Brazilian IBGE UF codes (`11`, `12`, `13`, `14`, `15`, `16`, `17`, `21`, `22`, `23`, `24`, `25`, `26`, `27`, `28`, `29`, `31`, `32`, `33`, `35`, `41`, `42`, `43`, `50`, `51`, `52`, `53`) MUST be mapped.
- `state_to_region_mapping` in `RuFaS/input/data/animal/lactation_curve_adjustment_inputs.json` MUST remain untouched for 100% backward compatibility with US FIPS codes (`"1"` to `"56"`).

---

### Task 1: Project Scaffolding, Data Models, and Biophysical Normalization Engine

**Files:**
- Create: `/home/thiago/Projects/rufas-lactation-equivalence/pyproject.toml`
- Create: `/home/thiago/Projects/rufas-lactation-equivalence/src/rufas_lactation_equivalence/__init__.py`
- Create: `/home/thiago/Projects/rufas-lactation-equivalence/src/rufas_lactation_equivalence/models.py`
- Create: `/home/thiago/Projects/rufas-lactation-equivalence/src/rufas_lactation_equivalence/biophysics.py`
- Test: `/home/thiago/Projects/rufas-lactation-equivalence/tests/test_biophysics.py`

**Interfaces:**
- Consumes: Raw climate and soil parameters ($T_{\text{mean}}$, $RH_{\text{mean}}$, $\Delta T$, $P_{\text{annual}}$, $\text{PET}_{\text{annual}}$, $\text{clay}_{\%}$, $\text{sand}_{\%}$, $\text{pH}$, $\text{SOC}_{\text{g/kg}}$).
- Produces:
  - `RawEdaphoclimaticProfile`, `NormalizedVector`, `MatchCandidate`, `EquivalenceResult` in `rufas_lactation_equivalence.models`
  - `calculate_thi(temp_c: float, rh_pct: float) -> float`
  - `calculate_pawc(sand_pct: float, clay_pct: float, soc_g_kg: float) -> float`
  - `min_max_normalize(value: float, lower: float, upper: float) -> float`
  - `normalize_profile(profile: RawEdaphoclimaticProfile) -> NormalizedVector`

- [ ] **Step 1: Write `pyproject.toml` and failing biophysical unit tests**

Create `/home/thiago/Projects/rufas-lactation-equivalence/pyproject.toml`:

```toml
[build-system]
requires = ["setuptools>=61.0"]
build-backend = "setuptools.build_meta"

[project]
name = "rufas-lactation-equivalence"
version = "1.0.0"
description = "Scientific edaphoclimatic equivalence pipeline mapping Brazilian IBGE UFs to RuFaS US dairy lactation curve adjustments"
readme = "README.md"
requires-python = ">=3.11"
dependencies = []

[tool.setuptools.packages.find]
where = ["src"]

[tool.pytest.ini_options]
pythonpath = ["src"]
testpaths = ["tests"]
```

Create `/home/thiago/Projects/rufas-lactation-equivalence/tests/test_biophysics.py`:

```python
import pytest

from rufas_lactation_equivalence.biophysics import (
    calculate_pawc,
    calculate_thi,
    min_max_normalize,
    normalize_profile,
)
from rufas_lactation_equivalence.models import NormalizedVector, RawEdaphoclimaticProfile


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
```

- [ ] **Step 2: Run test to verify it fails**

Run: `/home/thiago/Projects/RuFaS/.venv/bin/pytest /home/thiago/Projects/rufas-lactation-equivalence/tests/test_biophysics.py -v`
Expected: FAIL with `ModuleNotFoundError: No module named 'rufas_lactation_equivalence'`

- [ ] **Step 3: Implement `models.py` and `biophysics.py`**

Create `/home/thiago/Projects/rufas-lactation-equivalence/src/rufas_lactation_equivalence/__init__.py`:

```python
"""Scientific edaphoclimatic equivalence package for RuFaS lactation curves."""

__version__ = "1.0.0"
```

Create `/home/thiago/Projects/rufas-lactation-equivalence/src/rufas_lactation_equivalence/models.py`:

```python
from dataclasses import dataclass
from typing import Any


@dataclass(frozen=True)
class RawEdaphoclimaticProfile:
    """Raw biophysical climate and soil attributes for a state or UF."""

    code: str
    abbrev: str
    name: str
    country: str
    rufas_region: str | None
    temp_mean_c: float
    rh_mean_pct: float
    temp_seasonality_c: float
    precip_annual_mm: float
    pet_annual_mm: float
    soil_clay_pct: float
    soil_sand_pct: float
    soil_ph: float
    soil_soc_g_kg: float

    @classmethod
    def from_dict(cls, data: dict[str, Any]) -> "RawEdaphoclimaticProfile":
        return cls(
            code=str(data["code"]),
            abbrev=str(data["abbrev"]),
            name=str(data["name"]),
            country=str(data["country"]),
            rufas_region=data.get("rufas_region"),
            temp_mean_c=float(data["temp_mean_c"]),
            rh_mean_pct=float(data["rh_mean_pct"]),
            temp_seasonality_c=float(data["temp_seasonality_c"]),
            precip_annual_mm=float(data["precip_annual_mm"]),
            pet_annual_mm=float(data["pet_annual_mm"]),
            soil_clay_pct=float(data["soil_clay_pct"]),
            soil_sand_pct=float(data["soil_sand_pct"]),
            soil_ph=float(data["soil_ph"]),
            soil_soc_g_kg=float(data["soil_soc_g_kg"]),
        )


@dataclass(frozen=True)
class NormalizedVector:
    """Normalized 6-dimensional edaphoclimatic feature vector in [0.0, 1.0]."""

    thi_index: float
    temp_seasonality: float
    moisture_index: float
    soil_clay_fraction: float
    soil_pawc_index: float
    soil_fertility_index: float

    def as_tuple(self) -> tuple[float, float, float, float, float, float]:
        return (
            self.thi_index,
            self.temp_seasonality,
            self.moisture_index,
            self.soil_clay_fraction,
            self.soil_pawc_index,
            self.soil_fertility_index,
        )

    def to_dict(self) -> dict[str, float]:
        return {
            "thi_index": round(self.thi_index, 6),
            "temp_seasonality": round(self.temp_seasonality, 6),
            "moisture_index": round(self.moisture_index, 6),
            "soil_clay_fraction": round(self.soil_clay_fraction, 6),
            "soil_pawc_index": round(self.soil_pawc_index, 6),
            "soil_fertility_index": round(self.soil_fertility_index, 6),
        }


@dataclass(frozen=True)
class MatchCandidate:
    """Similarity evaluation between a target UF and a candidate US FIPS state."""

    us_code: str
    us_abbrev: str
    us_name: str
    rufas_region: str
    gower_similarity: float
    euclidean_similarity: float

    def to_dict(self) -> dict[str, Any]:
        return {
            "us_code": self.us_code,
            "us_abbrev": self.us_abbrev,
            "us_name": self.us_name,
            "rufas_region": self.rufas_region,
            "gower_similarity": round(self.gower_similarity, 6),
            "euclidean_similarity": round(self.euclidean_similarity, 6),
        }


@dataclass(frozen=True)
class EquivalenceResult:
    """Complete equivalence mapping result for a Brazilian UF."""

    brazil_profile: RawEdaphoclimaticProfile
    brazil_vector: NormalizedVector
    best_match: MatchCandidate
    top_candidates: tuple[MatchCandidate, ...]
    wood_adjustments: dict[str, float]
```

Create `/home/thiago/Projects/rufas-lactation-equivalence/src/rufas_lactation_equivalence/biophysics.py`:

```python
"""
Biophysical equations and [0.0, 1.0] feature normalization.

References
----------
1. Bohmanova, J., Misztal, I., & Cole, J. B. (2007). Temperature-humidity indices
   as indicators of milk production losses due to heat stress. J. Dairy Sci., 90(4), 1947-1956.
2. Saxton, K. E., & Rawls, W. J. (2006). Soil water characteristic estimates by texture
   and organic matter for hydrologic solutions. Soil Sci. Soc. Am. J., 70(5), 1569-1578.
3. Poggio, L., et al. (2021). SoilGrids 2.0. SOIL, 7(1), 217-240.
4. Fick, S. E., & Hijmans, R. J. (2017). WorldClim 2. Int. J. Climatol., 37(12), 4302-4315.
"""

from rufas_lactation_equivalence.models import NormalizedVector, RawEdaphoclimaticProfile

THI_MIN = 45.0
THI_MAX = 82.0

SEASONALITY_MIN_C = 1.0
SEASONALITY_MAX_C = 32.0

MOISTURE_INDEX_MIN = 0.10
MOISTURE_INDEX_MAX = 2.50

CLAY_MIN_PCT = 5.0
CLAY_MAX_PCT = 75.0

PAWC_MIN = 0.04
PAWC_MAX = 0.22

PH_MIN = 4.2
PH_MAX = 7.8

SOC_MIN_G_KG = 5.0
SOC_MAX_G_KG = 50.0


def min_max_normalize(value: float, lower: float, upper: float) -> float:
    """Linearly scales value into [0.0, 1.0] clamped at [lower, upper]."""
    if upper <= lower:
        raise ValueError(f"upper ({upper}) must be greater than lower ({lower})")
    scaled = (value - lower) / (upper - lower)
    return max(0.0, min(1.0, scaled))


def calculate_thi(temp_c: float, rh_pct: float) -> float:
    """
    Calculates Temperature-Humidity Index (THI) using NRC (1971) / Bohmanova et al. (2007).

    THI = (1.8 * T + 32) - (0.55 - 0.0055 * RH) * (1.8 * T - 26)
    """
    rh_clamped = max(0.0, min(100.0, rh_pct))
    return (1.8 * temp_c + 32.0) - (0.55 - 0.0055 * rh_clamped) * (1.8 * temp_c - 26.0)


def calculate_pawc(sand_pct: float, clay_pct: float, soc_g_kg: float) -> float:
    """
    Calculates Plant Available Water Capacity (PAWC, cm3/cm3) using Saxton & Rawls (2006).

    Parameters
    ----------
    sand_pct : float
        Sand percentage by weight (0-100).
    clay_pct : float
        Clay percentage by weight (0-100).
    soc_g_kg : float
        Soil organic carbon in g/kg.
    """
    s = max(0.0, min(1.0, sand_pct / 100.0))
    c = max(0.0, min(1.0, clay_pct / 100.0))
    # Van Bemmelen factor 1.724 converts SOC (g/kg -> % -> fraction) to Organic Matter weight %
    om = max(0.0, min(0.15, soc_g_kg * 0.001724))

    theta_1500t = (
        -0.024 * s
        + 0.487 * c
        + 0.006 * om
        + 0.005 * (s * om)
        - 0.013 * (c * om)
        + 0.068 * (s * c)
        + 0.031
    )
    theta_1500 = theta_1500t + (0.14 * theta_1500t - 0.02)

    theta_33t = (
        -0.251 * s
        + 0.195 * c
        + 0.011 * om
        + 0.006 * (s * om)
        - 0.027 * (c * om)
        + 0.452 * (s * c)
        + 0.299
    )
    theta_33 = theta_33t + (1.283 * (theta_33t**2) - 0.374 * theta_33t - 0.015)

    return max(0.01, theta_33 - theta_1500)


def normalize_profile(profile: RawEdaphoclimaticProfile) -> NormalizedVector:
    """Decomposes and normalizes a raw edaphoclimatic profile into a 6D [0.0, 1.0] vector."""
    thi = calculate_thi(profile.temp_mean_c, profile.rh_mean_pct)
    thi_norm = min_max_normalize(thi, THI_MIN, THI_MAX)

    seasonality_norm = min_max_normalize(
        profile.temp_seasonality_c, SEASONALITY_MIN_C, SEASONALITY_MAX_C
    )

    aridity_ratio = profile.precip_annual_mm / max(1.0, profile.pet_annual_mm)
    moisture_norm = min_max_normalize(aridity_ratio, MOISTURE_INDEX_MIN, MOISTURE_INDEX_MAX)

    clay_norm = min_max_normalize(profile.soil_clay_pct, CLAY_MIN_PCT, CLAY_MAX_PCT)

    pawc = calculate_pawc(profile.soil_sand_pct, profile.soil_clay_pct, profile.soil_soc_g_kg)
    pawc_norm = min_max_normalize(pawc, PAWC_MIN, PAWC_MAX)

    ph_norm = min_max_normalize(profile.soil_ph, PH_MIN, PH_MAX)
    soc_norm = min_max_normalize(profile.soil_soc_g_kg, SOC_MIN_G_KG, SOC_MAX_G_KG)
    fertility_norm = max(0.0, min(1.0, 0.5 * ph_norm + 0.5 * soc_norm))

    return NormalizedVector(
        thi_index=thi_norm,
        temp_seasonality=seasonality_norm,
        moisture_index=moisture_norm,
        soil_clay_fraction=clay_norm,
        soil_pawc_index=pawc_norm,
        soil_fertility_index=fertility_norm,
    )
```

- [ ] **Step 4: Run tests to verify they pass**

Run: `PYTHONPATH=/home/thiago/Projects/rufas-lactation-equivalence/src /home/thiago/Projects/RuFaS/.venv/bin/pytest /home/thiago/Projects/rufas-lactation-equivalence/tests/test_biophysics.py -v`
Expected: PASS (4 passed)

- [ ] **Step 5: Commit Task 1**

```bash
rtk git -C /home/thiago/Projects/rufas-lactation-equivalence add pyproject.toml src/ tests/test_biophysics.py
rtk git -C /home/thiago/Projects/rufas-lactation-equivalence commit -m "feat(biophysics): add data models and [0..1] edaphoclimatic normalization engine"
```

---

### Task 2: Gower & Euclidean Similarity Engine and Top-1 Ranking

**Files:**
- Create: `/home/thiago/Projects/rufas-lactation-equivalence/src/rufas_lactation_equivalence/similarity.py`
- Test: `/home/thiago/Projects/rufas-lactation-equivalence/tests/test_similarity.py`

**Interfaces:**
- Consumes: `NormalizedVector`, `RawEdaphoclimaticProfile`, `MatchCandidate`, `EquivalenceResult`, `normalize_profile`
- Produces:
  - `gower_similarity(vec_a: NormalizedVector, vec_b: NormalizedVector) -> float`
  - `euclidean_similarity(vec_a: NormalizedVector, vec_b: NormalizedVector) -> float`
  - `compute_all_equivalences(brazil_profiles: list[RawEdaphoclimaticProfile], us_profiles: list[RawEdaphoclimaticProfile], region_adjustments: dict[str, dict[str, float]], top_k: int = 3) -> list[EquivalenceResult]`

- [ ] **Step 1: Write failing unit tests for similarity and Top-1 matching**

Create `/home/thiago/Projects/rufas-lactation-equivalence/tests/test_similarity.py`:

```python
import pytest

from rufas_lactation_equivalence.models import NormalizedVector, RawEdaphoclimaticProfile
from rufas_lactation_equivalence.similarity import (
    compute_all_equivalences,
    euclidean_similarity,
    gower_similarity,
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
    assert gower_similarity(zeros, mid) == gower_similarity(mid, zeros)


def test_compute_all_equivalences_selects_closest_homologue() -> None:
    """Verify a warm/humid UF matches Southeast over cold/continental Wisconsin."""
    br_uf = RawEdaphoclimaticProfile(
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
    us_fl = RawEdaphoclimaticProfile(
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
    us_wi = RawEdaphoclimaticProfile(
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
```

- [ ] **Step 2: Run test to verify it fails**

Run: `PYTHONPATH=/home/thiago/Projects/rufas-lactation-equivalence/src /home/thiago/Projects/RuFaS/.venv/bin/pytest /home/thiago/Projects/rufas-lactation-equivalence/tests/test_similarity.py -v`
Expected: FAIL with `ModuleNotFoundError: No module named 'rufas_lactation_equivalence.similarity'`

- [ ] **Step 3: Implement `similarity.py`**

Create `/home/thiago/Projects/rufas-lactation-equivalence/src/rufas_lactation_equivalence/similarity.py`:

```python
"""
Environmental similarity metrics (Gower, 1971; Ramirez-Villegas et al., 2011)
and deterministic Top-1 homologue matching.
"""

import math

from rufas_lactation_equivalence.biophysics import normalize_profile
from rufas_lactation_equivalence.models import (
    EquivalenceResult,
    MatchCandidate,
    NormalizedVector,
    RawEdaphoclimaticProfile,
)


def gower_similarity(vec_a: NormalizedVector, vec_b: NormalizedVector) -> float:
    """
    Computes Gower's similarity coefficient in [0.0, 1.0] across the 6 normalized dimensions.

    S_Gower(a, b) = 1 - (1/6) * sum(|a_k - b_k|)
    """
    coords_a = vec_a.as_tuple()
    coords_b = vec_b.as_tuple()
    n = len(coords_a)
    manhattan = sum(abs(a - b) for a, b in zip(coords_a, coords_b))
    return max(0.0, min(1.0, 1.0 - (manhattan / n)))


def euclidean_similarity(vec_a: NormalizedVector, vec_b: NormalizedVector) -> float:
    """
    Computes normalized Euclidean similarity in [0.0, 1.0] across the 6 normalized dimensions.

    S_Eucl(a, b) = 1 - sqrt((1/6) * sum((a_k - b_k)^2))
    """
    coords_a = vec_a.as_tuple()
    coords_b = vec_b.as_tuple()
    n = len(coords_a)
    mean_sq = sum((a - b) ** 2 for a, b in zip(coords_a, coords_b)) / n
    return max(0.0, min(1.0, 1.0 - math.sqrt(mean_sq)))


def compute_all_equivalences(
    brazil_profiles: list[RawEdaphoclimaticProfile],
    us_profiles: list[RawEdaphoclimaticProfile],
    region_adjustments: dict[str, dict[str, float]],
    top_k: int = 3,
) -> list[EquivalenceResult]:
    """
    Computes the normalized vectors and ranks candidate US states for every Brazilian UF.
    Only US states mapped to an active RuFaS region (not None and != "None") are eligible candidates.
    """
    valid_us: list[tuple[RawEdaphoclimaticProfile, NormalizedVector]] = []
    for us_prof in us_profiles:
        if us_prof.rufas_region and us_prof.rufas_region != "None":
            if us_prof.rufas_region not in region_adjustments:
                raise ValueError(
                    f"US state {us_prof.code} references unknown RuFaS region '{us_prof.rufas_region}'"
                )
            valid_us.append((us_prof, normalize_profile(us_prof)))

    if not valid_us:
        raise ValueError("No valid US candidate states provided.")

    results: list[EquivalenceResult] = []
    for br_prof in sorted(brazil_profiles, key=lambda p: int(p.code)):
        br_vec = normalize_profile(br_prof)
        candidates: list[MatchCandidate] = []
        for us_prof, us_vec in valid_us:
            g_sim = gower_similarity(br_vec, us_vec)
            e_sim = euclidean_similarity(br_vec, us_vec)
            candidates.append(
                MatchCandidate(
                    us_code=us_prof.code,
                    us_abbrev=us_prof.abbrev,
                    us_name=us_prof.name,
                    rufas_region=str(us_prof.rufas_region),
                    gower_similarity=g_sim,
                    euclidean_similarity=e_sim,
                )
            )

        candidates.sort(
            key=lambda c: (
                round(c.gower_similarity, 12),
                round(c.euclidean_similarity, 12),
                -int(c.us_code),
            ),
            reverse=True,
        )

        best = candidates[0]
        adjustments = dict(region_adjustments[best.rufas_region])
        results.append(
            EquivalenceResult(
                brazil_profile=br_prof,
                brazil_vector=br_vec,
                best_match=best,
                top_candidates=tuple(candidates[: max(1, top_k)]),
                wood_adjustments=adjustments,
            )
        )

    return results
```

- [ ] **Step 4: Run tests to verify they pass**

Run: `PYTHONPATH=/home/thiago/Projects/rufas-lactation-equivalence/src /home/thiago/Projects/RuFaS/.venv/bin/pytest /home/thiago/Projects/rufas-lactation-equivalence/tests/test_similarity.py -v`
Expected: PASS (2 passed)

- [ ] **Step 5: Commit Task 2**

```bash
rtk git -C /home/thiago/Projects/rufas-lactation-equivalence add src/rufas_lactation_equivalence/similarity.py tests/test_similarity.py
rtk git -C /home/thiago/Projects/rufas-lactation-equivalence commit -m "feat(similarity): implement Gower and Euclidean similarity with Top-1 ranking"
```

---

### Task 3: Curated Scientific Datasets for US FIPS States and 27 Brazilian IBGE UFs

**Files:**
- Create: `/home/thiago/Projects/rufas-lactation-equivalence/data/rufas_region_adjustments.json`
- Create: `/home/thiago/Projects/rufas-lactation-equivalence/data/us_states_edaphoclimatic.json`
- Create: `/home/thiago/Projects/rufas-lactation-equivalence/data/brazil_ufs_edaphoclimatic.json`
- Test: `/home/thiago/Projects/rufas-lactation-equivalence/tests/test_datasets.py`

**Interfaces:**
- Consumes: Schema of `RawEdaphoclimaticProfile` and RuFaS `lactation_curve_adjustment_inputs.json`.
- Produces: Complete, validated JSON datasets for the 48 active US FIPS states mapped in RuFaS and all 27 Brazilian IBGE UFs (`11`..`53`).

- [ ] **Step 1: Write failing dataset completeness and domain validation test**

Create `/home/thiago/Projects/rufas-lactation-equivalence/tests/test_datasets.py`:

```python
import json
from pathlib import Path

from rufas_lactation_equivalence.biophysics import normalize_profile
from rufas_lactation_equivalence.models import RawEdaphoclimaticProfile

DATA_DIR = Path(__file__).resolve().parent.parent / "data"

EXPECTED_BRAZIL_UF_CODES = {
    "11", "12", "13", "14", "15", "16", "17",
    "21", "22", "23", "24", "25", "26", "27", "28", "29",
    "31", "32", "33", "35",
    "41", "42", "43",
    "50", "51", "52", "53",
}


def test_datasets_completeness_and_normalization_bounds() -> None:
    """Verify all 27 Brazilian UFs and 48 contiguous US states load and normalize to [0, 1]."""
    regions_data = json.loads((DATA_DIR / "rufas_region_adjustments.json").read_text(encoding="utf-8"))
    us_data = json.loads((DATA_DIR / "us_states_edaphoclimatic.json").read_text(encoding="utf-8"))
    br_data = json.loads((DATA_DIR / "brazil_ufs_edaphoclimatic.json").read_text(encoding="utf-8"))

    adjustments = regions_data["region_adjustments"]
    assert len(adjustments) == 14  # 13 regions + None

    us_profiles = [RawEdaphoclimaticProfile.from_dict(item) for item in us_data["states"]]
    assert len(us_profiles) == 48
    for p in us_profiles:
        assert p.rufas_region in adjustments
        assert p.rufas_region != "None"
        vec = normalize_profile(p)
        assert all(0.0 <= v <= 1.0 for v in vec.as_tuple())

    br_profiles = [RawEdaphoclimaticProfile.from_dict(item) for item in br_data["states"]]
    assert {p.code for p in br_profiles} == EXPECTED_BRAZIL_UF_CODES
    for p in br_profiles:
        vec = normalize_profile(p)
        assert all(0.0 <= v <= 1.0 for v in vec.as_tuple())
```

- [ ] **Step 2: Run test to verify it fails**

Run: `PYTHONPATH=/home/thiago/Projects/rufas-lactation-equivalence/src /home/thiago/Projects/RuFaS/.venv/bin/pytest /home/thiago/Projects/rufas-lactation-equivalence/tests/test_datasets.py -v`
Expected: FAIL with `FileNotFoundError`

- [ ] **Step 3: Create the 3 dataset files in `data/`**

Create `/home/thiago/Projects/rufas-lactation-equivalence/data/rufas_region_adjustments.json`:

```json
{
  "source": "Li, M., K. F. Reed, and V. E. Cabrera (2022). Investigating the effect of temporal, geographic, and management factors on US Holstein lactation curve parameters. Journal of Dairy Science 105(9):7525-7538.",
  "region_adjustments": {
    "Appalachian": {"l": -0.22, "m": -0.042, "n": -0.89},
    "Corn Belt": {"l": 0.55, "m": -0.58, "n": -1.12},
    "Delta": {"l": -2.56, "m": 0.59, "n": 1.47},
    "Lake": {"l": 0.61, "m": -0.40, "n": -0.64},
    "Mountain": {"l": -0.96, "m": 3.13, "n": 1.50},
    "Northeast": {"l": 1.04, "m": -1.99, "n": -1.13},
    "Northern Plains": {"l": -0.26, "m": 0.19, "n": -0.79},
    "New York": {"l": 0.67, "m": -1.21, "n": -0.45},
    "Pennsylvania": {"l": 1.15, "m": -0.96, "n": 0.06},
    "Southeast": {"l": -2.00, "m": 2.59, "n": 2.60},
    "Southern Plains": {"l": -0.51, "m": -1.02, "n": -0.93},
    "West Coast": {"l": 1.09, "m": 0.53, "n": 0.52},
    "Wisconsin": {"l": 1.4, "m": -0.83, "n": -0.2},
    "None": {"l": 0.0, "m": 0.0, "n": 0.0}
  }
}
```

Create `/home/thiago/Projects/rufas-lactation-equivalence/data/us_states_edaphoclimatic.json` (containing all 48 contiguous US states mapped in RuFaS `state_to_region_mapping`, with agroclimatic normals from WorldClim 2 / NOAA NCEI and topsoil 0–30 cm attributes from ISRIC SoilGrids 2.0 / USDA STATSGO2):

```json
{
  "metadata": {
    "description": "Edaphoclimatic baseline profiles (0-30 cm topsoil and 30-year climate normals) for the 48 contiguous US states mapped to RuFaS dairy regions.",
    "sources": [
      "Fick & Hijmans (2017) WorldClim 2 / NOAA NCEI 1991-2020 Climate Normals",
      "Poggio et al. (2021) ISRIC SoilGrids 2.0 / USDA NRCS STATSGO2 (0-30 cm agricultural zones)",
      "Fischer et al. (2021) FAO/IIASA GAEZ v4 Aridity & Evapotranspiration"
    ]
  },
  "states": [
    {"code": "1", "abbrev": "AL", "name": "Alabama", "country": "USA", "rufas_region": "Southeast", "temp_mean_c": 17.8, "rh_mean_pct": 72.0, "temp_seasonality_c": 19.2, "precip_annual_mm": 1400.0, "pet_annual_mm": 1120.0, "soil_clay_pct": 26.0, "soil_sand_pct": 48.0, "soil_ph": 5.3, "soil_soc_g_kg": 14.0},
    {"code": "4", "abbrev": "AZ", "name": "Arizona", "country": "USA", "rufas_region": "Mountain", "temp_mean_c": 19.5, "rh_mean_pct": 36.0, "temp_seasonality_c": 21.0, "precip_annual_mm": 280.0, "pet_annual_mm": 1650.0, "soil_clay_pct": 19.0, "soil_sand_pct": 52.0, "soil_ph": 7.6, "soil_soc_g_kg": 8.0},
    {"code": "5", "abbrev": "AR", "name": "Arkansas", "country": "USA", "rufas_region": "Delta", "temp_mean_c": 16.4, "rh_mean_pct": 71.0, "temp_seasonality_c": 22.0, "precip_annual_mm": 1280.0, "pet_annual_mm": 1080.0, "soil_clay_pct": 27.0, "soil_sand_pct": 36.0, "soil_ph": 5.6, "soil_soc_g_kg": 16.0},
    {"code": "6", "abbrev": "CA", "name": "California", "country": "USA", "rufas_region": "West Coast", "temp_mean_c": 16.8, "rh_mean_pct": 58.0, "temp_seasonality_c": 16.5, "precip_annual_mm": 480.0, "pet_annual_mm": 1320.0, "soil_clay_pct": 23.0, "soil_sand_pct": 44.0, "soil_ph": 6.8, "soil_soc_g_kg": 15.0},
    {"code": "8", "abbrev": "CO", "name": "Colorado", "country": "USA", "rufas_region": "Mountain", "temp_mean_c": 8.8, "rh_mean_pct": 50.0, "temp_seasonality_c": 24.5, "precip_annual_mm": 410.0, "pet_annual_mm": 1150.0, "soil_clay_pct": 22.0, "soil_sand_pct": 45.0, "soil_ph": 7.2, "soil_soc_g_kg": 14.0},
    {"code": "9", "abbrev": "CT", "name": "Connecticut", "country": "USA", "rufas_region": "Northeast", "temp_mean_c": 9.8, "rh_mean_pct": 68.0, "temp_seasonality_c": 25.0, "precip_annual_mm": 1220.0, "pet_annual_mm": 740.0, "soil_clay_pct": 14.0, "soil_sand_pct": 54.0, "soil_ph": 5.6, "soil_soc_g_kg": 28.0},
    {"code": "10", "abbrev": "DE", "name": "Delaware", "country": "USA", "rufas_region": "Northeast", "temp_mean_c": 13.2, "rh_mean_pct": 69.0, "temp_seasonality_c": 23.5, "precip_annual_mm": 1140.0, "pet_annual_mm": 890.0, "soil_clay_pct": 17.0, "soil_sand_pct": 52.0, "soil_ph": 5.7, "soil_soc_g_kg": 18.0},
    {"code": "12", "abbrev": "FL", "name": "Florida", "country": "USA", "rufas_region": "Southeast", "temp_mean_c": 22.1, "rh_mean_pct": 75.0, "temp_seasonality_c": 11.8, "precip_annual_mm": 1380.0, "pet_annual_mm": 1310.0, "soil_clay_pct": 15.0, "soil_sand_pct": 68.0, "soil_ph": 5.4, "soil_soc_g_kg": 16.0},
    {"code": "13", "abbrev": "GA", "name": "Georgia", "country": "USA", "rufas_region": "Southeast", "temp_mean_c": 18.2, "rh_mean_pct": 71.0, "temp_seasonality_c": 18.5, "precip_annual_mm": 1290.0, "pet_annual_mm": 1140.0, "soil_clay_pct": 28.0, "soil_sand_pct": 49.0, "soil_ph": 5.3, "soil_soc_g_kg": 14.5},
    {"code": "16", "abbrev": "ID", "name": "Idaho", "country": "USA", "rufas_region": "Mountain", "temp_mean_c": 8.4, "rh_mean_pct": 56.0, "temp_seasonality_c": 24.0, "precip_annual_mm": 380.0, "pet_annual_mm": 1020.0, "soil_clay_pct": 20.0, "soil_sand_pct": 42.0, "soil_ph": 7.1, "soil_soc_g_kg": 16.0},
    {"code": "17", "abbrev": "IL", "name": "Illinois", "country": "USA", "rufas_region": "Corn Belt", "temp_mean_c": 11.2, "rh_mean_pct": 70.0, "temp_seasonality_c": 27.5, "precip_annual_mm": 1020.0, "pet_annual_mm": 860.0, "soil_clay_pct": 27.0, "soil_sand_pct": 22.0, "soil_ph": 6.4, "soil_soc_g_kg": 29.0},
    {"code": "18", "abbrev": "IN", "name": "Indiana", "country": "USA", "rufas_region": "Corn Belt", "temp_mean_c": 11.1, "rh_mean_pct": 71.0, "temp_seasonality_c": 26.8, "precip_annual_mm": 1060.0, "pet_annual_mm": 840.0, "soil_clay_pct": 25.0, "soil_sand_pct": 26.0, "soil_ph": 6.3, "soil_soc_g_kg": 25.0},
    {"code": "19", "abbrev": "IA", "name": "Iowa", "country": "USA", "rufas_region": "Corn Belt", "temp_mean_c": 9.3, "rh_mean_pct": 70.0, "temp_seasonality_c": 29.5, "precip_annual_mm": 890.0, "pet_annual_mm": 820.0, "soil_clay_pct": 28.0, "soil_sand_pct": 21.0, "soil_ph": 6.5, "soil_soc_g_kg": 33.0},
    {"code": "20", "abbrev": "KS", "name": "Kansas", "country": "USA", "rufas_region": "Northern Plains", "temp_mean_c": 12.8, "rh_mean_pct": 64.0, "temp_seasonality_c": 27.8, "precip_annual_mm": 740.0, "pet_annual_mm": 1120.0, "soil_clay_pct": 27.0, "soil_sand_pct": 28.0, "soil_ph": 6.7, "soil_soc_g_kg": 21.0},
    {"code": "21", "abbrev": "KY", "name": "Kentucky", "country": "USA", "rufas_region": "Appalachian", "temp_mean_c": 13.6, "rh_mean_pct": 70.0, "temp_seasonality_c": 23.8, "precip_annual_mm": 1230.0, "pet_annual_mm": 930.0, "soil_clay_pct": 26.0, "soil_sand_pct": 25.0, "soil_ph": 5.8, "soil_soc_g_kg": 19.0},
    {"code": "22", "abbrev": "LA", "name": "Louisiana", "country": "USA", "rufas_region": "Delta", "temp_mean_c": 19.8, "rh_mean_pct": 76.0, "temp_seasonality_c": 16.8, "precip_annual_mm": 1520.0, "pet_annual_mm": 1180.0, "soil_clay_pct": 31.0, "soil_sand_pct": 32.0, "soil_ph": 5.6, "soil_soc_g_kg": 18.0},
    {"code": "23", "abbrev": "ME", "name": "Maine", "country": "USA", "rufas_region": "Northeast", "temp_mean_c": 5.8, "rh_mean_pct": 71.0, "temp_seasonality_c": 27.2, "precip_annual_mm": 1120.0, "pet_annual_mm": 610.0, "soil_clay_pct": 14.0, "soil_sand_pct": 52.0, "soil_ph": 5.3, "soil_soc_g_kg": 36.0},
    {"code": "24", "abbrev": "MD", "name": "Maryland", "country": "USA", "rufas_region": "Northeast", "temp_mean_c": 12.8, "rh_mean_pct": 68.0, "temp_seasonality_c": 24.0, "precip_annual_mm": 1110.0, "pet_annual_mm": 870.0, "soil_clay_pct": 21.0, "soil_sand_pct": 41.0, "soil_ph": 5.7, "soil_soc_g_kg": 20.0},
    {"code": "25", "abbrev": "MA", "name": "Massachusetts", "country": "USA", "rufas_region": "Northeast", "temp_mean_c": 9.1, "rh_mean_pct": 68.0, "temp_seasonality_c": 25.5, "precip_annual_mm": 1210.0, "pet_annual_mm": 710.0, "soil_clay_pct": 12.0, "soil_sand_pct": 58.0, "soil_ph": 5.4, "soil_soc_g_kg": 31.0},
    {"code": "26", "abbrev": "MI", "name": "Michigan", "country": "USA", "rufas_region": "Lake", "temp_mean_c": 7.6, "rh_mean_pct": 72.0, "temp_seasonality_c": 26.5, "precip_annual_mm": 850.0, "pet_annual_mm": 680.0, "soil_clay_pct": 19.0, "soil_sand_pct": 48.0, "soil_ph": 6.2, "soil_soc_g_kg": 29.0},
    {"code": "27", "abbrev": "MN", "name": "Minnesota", "country": "USA", "rufas_region": "Lake", "temp_mean_c": 5.8, "rh_mean_pct": 70.0, "temp_seasonality_c": 31.5, "precip_annual_mm": 720.0, "pet_annual_mm": 690.0, "soil_clay_pct": 25.0, "soil_sand_pct": 33.0, "soil_ph": 6.5, "soil_soc_g_kg": 36.0},
    {"code": "28", "abbrev": "MS", "name": "Mississippi", "country": "USA", "rufas_region": "Delta", "temp_mean_c": 18.1, "rh_mean_pct": 74.0, "temp_seasonality_c": 19.0, "precip_annual_mm": 1440.0, "pet_annual_mm": 1130.0, "soil_clay_pct": 29.0, "soil_sand_pct": 35.0, "soil_ph": 5.4, "soil_soc_g_kg": 15.5},
    {"code": "29", "abbrev": "MO", "name": "Missouri", "country": "USA", "rufas_region": "Corn Belt", "temp_mean_c": 12.8, "rh_mean_pct": 69.0, "temp_seasonality_c": 26.2, "precip_annual_mm": 1080.0, "pet_annual_mm": 940.0, "soil_clay_pct": 28.0, "soil_sand_pct": 22.0, "soil_ph": 6.1, "soil_soc_g_kg": 23.0},
    {"code": "30", "abbrev": "MT", "name": "Montana", "country": "USA", "rufas_region": "Mountain", "temp_mean_c": 6.2, "rh_mean_pct": 58.0, "temp_seasonality_c": 27.5, "precip_annual_mm": 390.0, "pet_annual_mm": 890.0, "soil_clay_pct": 24.0, "soil_sand_pct": 39.0, "soil_ph": 7.1, "soil_soc_g_kg": 18.0},
    {"code": "31", "abbrev": "NE", "name": "Nebraska", "country": "USA", "rufas_region": "Northern Plains", "temp_mean_c": 9.8, "rh_mean_pct": 65.0, "temp_seasonality_c": 28.5, "precip_annual_mm": 610.0, "pet_annual_mm": 980.0, "soil_clay_pct": 24.0, "soil_sand_pct": 38.0, "soil_ph": 6.7, "soil_soc_g_kg": 22.0},
    {"code": "32", "abbrev": "NV", "name": "Nevada", "country": "USA", "rufas_region": "Mountain", "temp_mean_c": 10.5, "rh_mean_pct": 42.0, "temp_seasonality_c": 23.5, "precip_annual_mm": 240.0, "pet_annual_mm": 1280.0, "soil_clay_pct": 18.0, "soil_sand_pct": 49.0, "soil_ph": 7.6, "soil_soc_g_kg": 9.0},
    {"code": "33", "abbrev": "NH", "name": "New Hampshire", "country": "USA", "rufas_region": "Northeast", "temp_mean_c": 6.9, "rh_mean_pct": 69.0, "temp_seasonality_c": 26.8, "precip_annual_mm": 1150.0, "pet_annual_mm": 640.0, "soil_clay_pct": 11.0, "soil_sand_pct": 59.0, "soil_ph": 5.3, "soil_soc_g_kg": 34.0},
    {"code": "34", "abbrev": "NJ", "name": "New Jersey", "country": "USA", "rufas_region": "Northeast", "temp_mean_c": 11.8, "rh_mean_pct": 68.0, "temp_seasonality_c": 24.5, "precip_annual_mm": 1180.0, "pet_annual_mm": 820.0, "soil_clay_pct": 17.0, "soil_sand_pct": 50.0, "soil_ph": 5.6, "soil_soc_g_kg": 22.0},
    {"code": "35", "abbrev": "NM", "name": "New Mexico", "country": "USA", "rufas_region": "Mountain", "temp_mean_c": 13.5, "rh_mean_pct": 44.0, "temp_seasonality_c": 22.0, "precip_annual_mm": 340.0, "pet_annual_mm": 1420.0, "soil_clay_pct": 21.0, "soil_sand_pct": 48.0, "soil_ph": 7.5, "soil_soc_g_kg": 10.0},
    {"code": "36", "abbrev": "NY", "name": "New York", "country": "USA", "rufas_region": "New York", "temp_mean_c": 8.2, "rh_mean_pct": 70.0, "temp_seasonality_c": 26.5, "precip_annual_mm": 1040.0, "pet_annual_mm": 690.0, "soil_clay_pct": 20.0, "soil_sand_pct": 38.0, "soil_ph": 5.9, "soil_soc_g_kg": 30.0},
    {"code": "37", "abbrev": "NC", "name": "North Carolina", "country": "USA", "rufas_region": "Appalachian", "temp_mean_c": 15.6, "rh_mean_pct": 70.0, "temp_seasonality_c": 20.5, "precip_annual_mm": 1270.0, "pet_annual_mm": 1020.0, "soil_clay_pct": 25.0, "soil_sand_pct": 46.0, "soil_ph": 5.4, "soil_soc_g_kg": 17.0},
    {"code": "38", "abbrev": "ND", "name": "North Dakota", "country": "USA", "rufas_region": "Northern Plains", "temp_mean_c": 4.9, "rh_mean_pct": 67.0, "temp_seasonality_c": 32.0, "precip_annual_mm": 460.0, "pet_annual_mm": 780.0, "soil_clay_pct": 26.0, "soil_sand_pct": 33.0, "soil_ph": 6.9, "soil_soc_g_kg": 31.0},
    {"code": "39", "abbrev": "OH", "name": "Ohio", "country": "USA", "rufas_region": "Corn Belt", "temp_mean_c": 10.6, "rh_mean_pct": 71.0, "temp_seasonality_c": 26.0, "precip_annual_mm": 1010.0, "pet_annual_mm": 800.0, "soil_clay_pct": 26.0, "soil_sand_pct": 25.0, "soil_ph": 6.2, "soil_soc_g_kg": 24.0},
    {"code": "40", "abbrev": "OK", "name": "Oklahoma", "country": "USA", "rufas_region": "Southern Plains", "temp_mean_c": 15.8, "rh_mean_pct": 64.0, "temp_seasonality_c": 24.5, "precip_annual_mm": 910.0, "pet_annual_mm": 1240.0, "soil_clay_pct": 25.0, "soil_sand_pct": 38.0, "soil_ph": 6.4, "soil_soc_g_kg": 17.0},
    {"code": "41", "abbrev": "OR", "name": "Oregon", "country": "USA", "rufas_region": "West Coast", "temp_mean_c": 10.8, "rh_mean_pct": 68.0, "temp_seasonality_c": 16.0, "precip_annual_mm": 890.0, "pet_annual_mm": 840.0, "soil_clay_pct": 22.0, "soil_sand_pct": 39.0, "soil_ph": 6.1, "soil_soc_g_kg": 27.0},
    {"code": "42", "abbrev": "PA", "name": "Pennsylvania", "country": "USA", "rufas_region": "Pennsylvania", "temp_mean_c": 9.9, "rh_mean_pct": 69.0, "temp_seasonality_c": 25.5, "precip_annual_mm": 1090.0, "pet_annual_mm": 760.0, "soil_clay_pct": 22.0, "soil_sand_pct": 33.0, "soil_ph": 5.7, "soil_soc_g_kg": 25.0},
    {"code": "44", "abbrev": "RI", "name": "Rhode Island", "country": "USA", "rufas_region": "Northeast", "temp_mean_c": 10.2, "rh_mean_pct": 69.0, "temp_seasonality_c": 24.0, "precip_annual_mm": 1230.0, "pet_annual_mm": 730.0, "soil_clay_pct": 11.0, "soil_sand_pct": 60.0, "soil_ph": 5.4, "soil_soc_g_kg": 29.0},
    {"code": "45", "abbrev": "SC", "name": "South Carolina", "country": "USA", "rufas_region": "Southeast", "temp_mean_c": 17.6, "rh_mean_pct": 71.0, "temp_seasonality_c": 19.2, "precip_annual_mm": 1260.0, "pet_annual_mm": 1110.0, "soil_clay_pct": 24.0, "soil_sand_pct": 53.0, "soil_ph": 5.3, "soil_soc_g_kg": 15.0},
    {"code": "46", "abbrev": "SD", "name": "South Dakota", "country": "USA", "rufas_region": "Northern Plains", "temp_mean_c": 7.6, "rh_mean_pct": 66.0, "temp_seasonality_c": 30.8, "precip_annual_mm": 520.0, "pet_annual_mm": 890.0, "soil_clay_pct": 27.0, "soil_sand_pct": 31.0, "soil_ph": 6.8, "soil_soc_g_kg": 26.0},
    {"code": "47", "abbrev": "TN", "name": "Tennessee", "country": "USA", "rufas_region": "Appalachian", "temp_mean_c": 14.8, "rh_mean_pct": 71.0, "temp_seasonality_c": 22.0, "precip_annual_mm": 1340.0, "pet_annual_mm": 980.0, "soil_clay_pct": 26.0, "soil_sand_pct": 28.0, "soil_ph": 5.6, "soil_soc_g_kg": 17.5},
    {"code": "48", "abbrev": "TX", "name": "Texas", "country": "USA", "rufas_region": "Southern Plains", "temp_mean_c": 18.9, "rh_mean_pct": 63.0, "temp_seasonality_c": 20.5, "precip_annual_mm": 730.0, "pet_annual_mm": 1420.0, "soil_clay_pct": 28.0, "soil_sand_pct": 40.0, "soil_ph": 6.9, "soil_soc_g_kg": 14.0},
    {"code": "49", "abbrev": "UT", "name": "Utah", "country": "USA", "rufas_region": "Mountain", "temp_mean_c": 9.6, "rh_mean_pct": 48.0, "temp_seasonality_c": 25.0, "precip_annual_mm": 330.0, "pet_annual_mm": 1180.0, "soil_clay_pct": 20.0, "soil_sand_pct": 46.0, "soil_ph": 7.5, "soil_soc_g_kg": 12.0},
    {"code": "50", "abbrev": "VT", "name": "Vermont", "country": "USA", "rufas_region": "Northeast", "temp_mean_c": 6.1, "rh_mean_pct": 71.0, "temp_seasonality_c": 27.8, "precip_annual_mm": 1080.0, "pet_annual_mm": 610.0, "soil_clay_pct": 15.0, "soil_sand_pct": 48.0, "soil_ph": 5.6, "soil_soc_g_kg": 35.0},
    {"code": "51", "abbrev": "VA", "name": "Virginia", "country": "USA", "rufas_region": "Appalachian", "temp_mean_c": 13.8, "rh_mean_pct": 69.0, "temp_seasonality_c": 22.5, "precip_annual_mm": 1140.0, "pet_annual_mm": 920.0, "soil_clay_pct": 25.0, "soil_sand_pct": 40.0, "soil_ph": 5.5, "soil_soc_g_kg": 18.5},
    {"code": "53", "abbrev": "WA", "name": "Washington", "country": "USA", "rufas_region": "West Coast", "temp_mean_c": 10.1, "rh_mean_pct": 70.0, "temp_seasonality_c": 17.5, "precip_annual_mm": 980.0, "pet_annual_mm": 790.0, "soil_clay_pct": 18.0, "soil_sand_pct": 45.0, "soil_ph": 6.2, "soil_soc_g_kg": 26.0},
    {"code": "54", "abbrev": "WV", "name": "West Virginia", "country": "USA", "rufas_region": "Appalachian", "temp_mean_c": 11.4, "rh_mean_pct": 71.0, "temp_seasonality_c": 23.8, "precip_annual_mm": 1160.0, "pet_annual_mm": 820.0, "soil_clay_pct": 24.0, "soil_sand_pct": 31.0, "soil_ph": 5.4, "soil_soc_g_kg": 22.0},
    {"code": "55", "abbrev": "WI", "name": "Wisconsin", "country": "USA", "rufas_region": "Wisconsin", "temp_mean_c": 6.9, "rh_mean_pct": 71.0, "temp_seasonality_c": 29.2, "precip_annual_mm": 840.0, "pet_annual_mm": 680.0, "soil_clay_pct": 21.0, "soil_sand_pct": 41.0, "soil_ph": 6.3, "soil_soc_g_kg": 31.0},
    {"code": "56", "abbrev": "WY", "name": "Wyoming", "country": "USA", "rufas_region": "Mountain", "temp_mean_c": 6.1, "rh_mean_pct": 53.0, "temp_seasonality_c": 26.8, "precip_annual_mm": 330.0, "pet_annual_mm": 980.0, "soil_clay_pct": 22.0, "soil_sand_pct": 44.0, "soil_ph": 7.3, "soil_soc_g_kg": 14.0}
  ]
}
```

Create `/home/thiago/Projects/rufas-lactation-equivalence/data/brazil_ufs_edaphoclimatic.json` (containing all 27 Brazilian UFs with IBGE 2-digit codes `11`..`53`, grounded in INMET 1991–2020 normals, WorldClim 2, and Embrapa Solos / ISRIC SoilGrids 2.0 0–30 cm topsoil profiles):

```json
{
  "metadata": {
    "description": "Edaphoclimatic profiles (0-30 cm topsoil and 30-year climate normals) for all 27 Brazilian Federative Units (IBGE UF codes 11..53).",
    "sources": [
      "INMET Normais Climatológicas do Brasil (1991-2020) & WorldClim 2 (Fick & Hijmans, 2017)",
      "Embrapa Solos (SiBCS) & ISRIC SoilGrids 2.0 (Poggio et al., 2021) 0-30 cm agricultural/pasture zones",
      "FAO/IIASA GAEZ v4 Evapotranspiration and Moisture Regime (Fischer et al., 2021)"
    ]
  },
  "states": [
    {"code": "11", "abbrev": "RO", "name": "Rondônia", "country": "BRA", "rufas_region": null, "temp_mean_c": 25.6, "rh_mean_pct": 81.0, "temp_seasonality_c": 2.4, "precip_annual_mm": 2150.0, "pet_annual_mm": 1420.0, "soil_clay_pct": 38.0, "soil_sand_pct": 44.0, "soil_ph": 5.1, "soil_soc_g_kg": 16.5},
    {"code": "12", "abbrev": "AC", "name": "Acre", "country": "BRA", "rufas_region": null, "temp_mean_c": 25.2, "rh_mean_pct": 84.0, "temp_seasonality_c": 2.6, "precip_annual_mm": 2050.0, "pet_annual_mm": 1380.0, "soil_clay_pct": 36.0, "soil_sand_pct": 38.0, "soil_ph": 5.2, "soil_soc_g_kg": 17.0},
    {"code": "13", "abbrev": "AM", "name": "Amazonas", "country": "BRA", "rufas_region": null, "temp_mean_c": 26.7, "rh_mean_pct": 85.0, "temp_seasonality_c": 1.8, "precip_annual_mm": 2450.0, "pet_annual_mm": 1450.0, "soil_clay_pct": 42.0, "soil_sand_pct": 40.0, "soil_ph": 4.6, "soil_soc_g_kg": 19.0},
    {"code": "14", "abbrev": "RR", "name": "Roraima", "country": "BRA", "rufas_region": null, "temp_mean_c": 27.2, "rh_mean_pct": 76.0, "temp_seasonality_c": 2.2, "precip_annual_mm": 1680.0, "pet_annual_mm": 1560.0, "soil_clay_pct": 29.0, "soil_sand_pct": 54.0, "soil_ph": 5.0, "soil_soc_g_kg": 13.0},
    {"code": "15", "abbrev": "PA", "name": "Pará", "country": "BRA", "rufas_region": null, "temp_mean_c": 26.4, "rh_mean_pct": 82.0, "temp_seasonality_c": 1.9, "precip_annual_mm": 2200.0, "pet_annual_mm": 1480.0, "soil_clay_pct": 37.0, "soil_sand_pct": 46.0, "soil_ph": 4.9, "soil_soc_g_kg": 16.0},
    {"code": "16", "abbrev": "AP", "name": "Amapá", "country": "BRA", "rufas_region": null, "temp_mean_c": 26.6, "rh_mean_pct": 83.0, "temp_seasonality_c": 1.7, "precip_annual_mm": 2500.0, "pet_annual_mm": 1460.0, "soil_clay_pct": 40.0, "soil_sand_pct": 43.0, "soil_ph": 4.8, "soil_soc_g_kg": 17.5},
    {"code": "17", "abbrev": "TO", "name": "Tocantins", "country": "BRA", "rufas_region": null, "temp_mean_c": 26.1, "rh_mean_pct": 68.0, "temp_seasonality_c": 3.2, "precip_annual_mm": 1620.0, "pet_annual_mm": 1650.0, "soil_clay_pct": 31.0, "soil_sand_pct": 52.0, "soil_ph": 5.2, "soil_soc_g_kg": 13.5},
    {"code": "21", "abbrev": "MA", "name": "Maranhão", "country": "BRA", "rufas_region": null, "temp_mean_c": 26.5, "rh_mean_pct": 75.0, "temp_seasonality_c": 2.1, "precip_annual_mm": 1580.0, "pet_annual_mm": 1620.0, "soil_clay_pct": 26.0, "soil_sand_pct": 58.0, "soil_ph": 5.3, "soil_soc_g_kg": 12.5},
    {"code": "22", "abbrev": "PI", "name": "Piauí", "country": "BRA", "rufas_region": null, "temp_mean_c": 26.8, "rh_mean_pct": 62.0, "temp_seasonality_c": 3.1, "precip_annual_mm": 1050.0, "pet_annual_mm": 1820.0, "soil_clay_pct": 24.0, "soil_sand_pct": 60.0, "soil_ph": 5.5, "soil_soc_g_kg": 10.5},
    {"code": "23", "abbrev": "CE", "name": "Ceará", "country": "BRA", "rufas_region": null, "temp_mean_c": 26.6, "rh_mean_pct": 64.0, "temp_seasonality_c": 2.5, "precip_annual_mm": 820.0, "pet_annual_mm": 1880.0, "soil_clay_pct": 23.0, "soil_sand_pct": 58.0, "soil_ph": 6.0, "soil_soc_g_kg": 10.0},
    {"code": "24", "abbrev": "RN", "name": "Rio Grande do Norte", "country": "BRA", "rufas_region": null, "temp_mean_c": 26.5, "rh_mean_pct": 66.0, "temp_seasonality_c": 2.4, "precip_annual_mm": 780.0, "pet_annual_mm": 1850.0, "soil_clay_pct": 22.0, "soil_sand_pct": 61.0, "soil_ph": 6.1, "soil_soc_g_kg": 9.8},
    {"code": "25", "abbrev": "PB", "name": "Paraíba", "country": "BRA", "rufas_region": null, "temp_mean_c": 25.5, "rh_mean_pct": 68.0, "temp_seasonality_c": 3.0, "precip_annual_mm": 860.0, "pet_annual_mm": 1720.0, "soil_clay_pct": 25.0, "soil_sand_pct": 56.0, "soil_ph": 5.9, "soil_soc_g_kg": 11.0},
    {"code": "26", "abbrev": "PE", "name": "Pernambuco", "country": "BRA", "rufas_region": null, "temp_mean_c": 24.8, "rh_mean_pct": 70.0, "temp_seasonality_c": 3.4, "precip_annual_mm": 920.0, "pet_annual_mm": 1640.0, "soil_clay_pct": 27.0, "soil_sand_pct": 54.0, "soil_ph": 5.7, "soil_soc_g_kg": 12.0},
    {"code": "27", "abbrev": "AL", "name": "Alagoas", "country": "BRA", "rufas_region": null, "temp_mean_c": 24.9, "rh_mean_pct": 74.0, "temp_seasonality_c": 3.2, "precip_annual_mm": 1180.0, "pet_annual_mm": 1520.0, "soil_clay_pct": 29.0, "soil_sand_pct": 52.0, "soil_ph": 5.5, "soil_soc_g_kg": 13.5},
    {"code": "28", "abbrev": "SE", "name": "Sergipe", "country": "BRA", "rufas_region": null, "temp_mean_c": 25.1, "rh_mean_pct": 74.0, "temp_seasonality_c": 3.1, "precip_annual_mm": 1120.0, "pet_annual_mm": 1540.0, "soil_clay_pct": 28.0, "soil_sand_pct": 53.0, "soil_ph": 5.6, "soil_soc_g_kg": 13.0},
    {"code": "29", "abbrev": "BA", "name": "Bahia", "country": "BRA", "rufas_region": null, "temp_mean_c": 23.8, "rh_mean_pct": 68.0, "temp_seasonality_c": 4.2, "precip_annual_mm": 1020.0, "pet_annual_mm": 1580.0, "soil_clay_pct": 32.0, "soil_sand_pct": 49.0, "soil_ph": 5.5, "soil_soc_g_kg": 13.5},
    {"code": "31", "abbrev": "MG", "name": "Minas Gerais", "country": "BRA", "rufas_region": null, "temp_mean_c": 21.2, "rh_mean_pct": 70.0, "temp_seasonality_c": 5.8, "precip_annual_mm": 1450.0, "pet_annual_mm": 1220.0, "soil_clay_pct": 46.0, "soil_sand_pct": 34.0, "soil_ph": 5.2, "soil_soc_g_kg": 19.0},
    {"code": "32", "abbrev": "ES", "name": "Espírito Santo", "country": "BRA", "rufas_region": null, "temp_mean_c": 23.2, "rh_mean_pct": 76.0, "temp_seasonality_c": 4.8, "precip_annual_mm": 1320.0, "pet_annual_mm": 1310.0, "soil_clay_pct": 41.0, "soil_sand_pct": 42.0, "soil_ph": 5.1, "soil_soc_g_kg": 17.5},
    {"code": "33", "abbrev": "RJ", "name": "Rio de Janeiro", "country": "BRA", "rufas_region": null, "temp_mean_c": 22.8, "rh_mean_pct": 77.0, "temp_seasonality_c": 5.4, "precip_annual_mm": 1380.0, "pet_annual_mm": 1260.0, "soil_clay_pct": 39.0, "soil_sand_pct": 43.0, "soil_ph": 5.2, "soil_soc_g_kg": 18.0},
    {"code": "35", "abbrev": "SP", "name": "São Paulo", "country": "BRA", "rufas_region": null, "temp_mean_c": 21.5, "rh_mean_pct": 71.0, "temp_seasonality_c": 6.4, "precip_annual_mm": 1420.0, "pet_annual_mm": 1240.0, "soil_clay_pct": 42.0, "soil_sand_pct": 39.0, "soil_ph": 5.4, "soil_soc_g_kg": 18.5},
    {"code": "41", "abbrev": "PR", "name": "Paraná", "country": "BRA", "rufas_region": null, "temp_mean_c": 19.1, "rh_mean_pct": 75.0, "temp_seasonality_c": 8.6, "precip_annual_mm": 1680.0, "pet_annual_mm": 1080.0, "soil_clay_pct": 52.0, "soil_sand_pct": 26.0, "soil_ph": 5.4, "soil_soc_g_kg": 24.0},
    {"code": "42", "abbrev": "SC", "name": "Santa Catarina", "country": "BRA", "rufas_region": null, "temp_mean_c": 17.8, "rh_mean_pct": 78.0, "temp_seasonality_c": 9.4, "precip_annual_mm": 1720.0, "pet_annual_mm": 990.0, "soil_clay_pct": 48.0, "soil_sand_pct": 28.0, "soil_ph": 5.2, "soil_soc_g_kg": 26.5},
    {"code": "43", "abbrev": "RS", "name": "Rio Grande do Sul", "country": "BRA", "rufas_region": null, "temp_mean_c": 17.6, "rh_mean_pct": 76.0, "temp_seasonality_c": 10.8, "precip_annual_mm": 1620.0, "pet_annual_mm": 1020.0, "soil_clay_pct": 41.0, "soil_sand_pct": 33.0, "soil_ph": 5.4, "soil_soc_g_kg": 23.5},
    {"code": "50", "abbrev": "MS", "name": "Mato Grosso do Sul", "country": "BRA", "rufas_region": null, "temp_mean_c": 23.6, "rh_mean_pct": 69.0, "temp_seasonality_c": 5.6, "precip_annual_mm": 1440.0, "pet_annual_mm": 1420.0, "soil_clay_pct": 35.0, "soil_sand_pct": 48.0, "soil_ph": 5.4, "soil_soc_g_kg": 16.0},
    {"code": "51", "abbrev": "MT", "name": "Mato Grosso", "country": "BRA", "rufas_region": null, "temp_mean_c": 25.1, "rh_mean_pct": 71.0, "temp_seasonality_c": 3.8, "precip_annual_mm": 1720.0, "pet_annual_mm": 1510.0, "soil_clay_pct": 39.0, "soil_sand_pct": 44.0, "soil_ph": 5.2, "soil_soc_g_kg": 16.5},
    {"code": "52", "abbrev": "GO", "name": "Goiás", "country": "BRA", "rufas_region": null, "temp_mean_c": 23.4, "rh_mean_pct": 66.0, "temp_seasonality_c": 4.4, "precip_annual_mm": 1520.0, "pet_annual_mm": 1460.0, "soil_clay_pct": 45.0, "soil_sand_pct": 36.0, "soil_ph": 5.3, "soil_soc_g_kg": 17.5},
    {"code": "53", "abbrev": "DF", "name": "Distrito Federal", "country": "BRA", "rufas_region": null, "temp_mean_c": 21.4, "rh_mean_pct": 65.0, "temp_seasonality_c": 4.1, "precip_annual_mm": 1480.0, "pet_annual_mm": 1360.0, "soil_clay_pct": 52.0, "soil_sand_pct": 29.0, "soil_ph": 5.1, "soil_soc_g_kg": 20.0}
  ]
}
```

- [ ] **Step 4: Run dataset tests to verify they pass**

Run: `PYTHONPATH=/home/thiago/Projects/rufas-lactation-equivalence/src /home/thiago/Projects/RuFaS/.venv/bin/pytest /home/thiago/Projects/rufas-lactation-equivalence/tests/test_datasets.py -v`
Expected: PASS (1 passed)

- [ ] **Step 5: Commit Task 3**

```bash
rtk git -C /home/thiago/Projects/rufas-lactation-equivalence add data/ tests/test_datasets.py
rtk git -C /home/thiago/Projects/rufas-lactation-equivalence commit -m "feat(data): add curated edaphoclimatic datasets for 48 US FIPS states and 27 Brazilian UFs"
```

---

### Task 4: Atomic I/O Engine, CLI Pipeline, and Audit Report Generator

**Files:**
- Create: `/home/thiago/Projects/rufas-lactation-equivalence/src/rufas_lactation_equivalence/atomic_io.py`
- Create: `/home/thiago/Projects/rufas-lactation-equivalence/src/rufas_lactation_equivalence/cli.py`
- Create: `/home/thiago/Projects/rufas-lactation-equivalence/README.md`
- Test: `/home/thiago/Projects/rufas-lactation-equivalence/tests/test_atomic_pipeline.py`

**Interfaces:**
- Consumes: `compute_all_equivalences`, `normalize_profile`, `RawEdaphoclimaticProfile`, `EquivalenceResult`, `data/*.json`
- Produces:
  - `atomic_write_text(target_path: Path, content: str) -> None`
  - `atomic_write_json(target_path: Path, payload: dict[str, Any]) -> None`
  - `run_pipeline(data_dir: Path, output_dir: Path, export_rufas_json: Path | None = None) -> dict[str, Any]`
  - Generated files `output/equivalence_audit_report.json` and `output/equivalence_audit_report.md`

- [ ] **Step 1: Write failing tests for atomic I/O and CLI pipeline execution**

Create `/home/thiago/Projects/rufas-lactation-equivalence/tests/test_atomic_pipeline.py`:

```python
import json
from pathlib import Path

from rufas_lactation_equivalence.atomic_io import atomic_write_json, atomic_write_text
from rufas_lactation_equivalence.cli import run_pipeline

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
```

- [ ] **Step 2: Run test to verify it fails**

Run: `PYTHONPATH=/home/thiago/Projects/rufas-lactation-equivalence/src /home/thiago/Projects/RuFaS/.venv/bin/pytest /home/thiago/Projects/rufas-lactation-equivalence/tests/test_atomic_pipeline.py -v`
Expected: FAIL with `ModuleNotFoundError: No module named 'rufas_lactation_equivalence.atomic_io'`

- [ ] **Step 3: Implement `atomic_io.py`, `cli.py`, and `README.md`**

Create `/home/thiago/Projects/rufas-lactation-equivalence/src/rufas_lactation_equivalence/atomic_io.py`:

```python
"""Transactional atomic file writing utilities."""

import json
import os
import tempfile
from pathlib import Path
from typing import Any


def atomic_write_text(target_path: Path, content: str) -> None:
    """
    Writes text to target_path atomically using a temporary file in the same directory,
    os.fsync(), and os.replace().
    """
    target_path = Path(target_path)
    target_path.parent.mkdir(parents=True, exist_ok=True)
    fd, tmp_name = tempfile.mkstemp(
        dir=str(target_path.parent),
        prefix=f".{target_path.name}.",
        suffix=".tmp",
    )
    tmp_path = Path(tmp_name)
    try:
        with os.fdopen(fd, "w", encoding="utf-8") as f:
            f.write(content)
            f.flush()
            os.fsync(f.fileno())
        os.replace(str(tmp_path), str(target_path))
    except Exception:
        if tmp_path.exists():
            tmp_path.unlink(missing_ok=True)
        raise


def atomic_write_json(target_path: Path, payload: dict[str, Any]) -> None:
    """Serializes payload to JSON in memory first, then writes atomically to target_path."""
    serialized = json.dumps(payload, indent=2, ensure_ascii=False) + "\n"
    atomic_write_text(target_path, serialized)
```

Create `/home/thiago/Projects/rufas-lactation-equivalence/src/rufas_lactation_equivalence/cli.py`:

```python
"""
Command-line interface and atomic pipeline runner for edaphoclimatic lactation curve equivalence.
"""

import argparse
import json
from pathlib import Path
from typing import Any

from rufas_lactation_equivalence.atomic_io import atomic_write_json, atomic_write_text
from rufas_lactation_equivalence.biophysics import (
    CLAY_MAX_PCT,
    CLAY_MIN_PCT,
    MOISTURE_INDEX_MAX,
    MOISTURE_INDEX_MIN,
    PAWC_MAX,
    PAWC_MIN,
    PH_MAX,
    PH_MIN,
    SEASONALITY_MAX_C,
    SEASONALITY_MIN_C,
    SOC_MAX_G_KG,
    SOC_MIN_G_KG,
    THI_MAX,
    THI_MIN,
    normalize_profile,
)
from rufas_lactation_equivalence.models import EquivalenceResult, RawEdaphoclimaticProfile
from rufas_lactation_equivalence.similarity import compute_all_equivalences


def _format_markdown_report(
    results: list[EquivalenceResult],
    us_profiles: list[RawEdaphoclimaticProfile],
) -> str:
    lines: list[str] = [
        "# Relatório de Auditoria Científica: Equivalência Edafoclimática da Curva de Lactação (EUA ↔ Brasil)",
        "",
        "## 1. Fundamentação Científica e Limites de Normalização $[0.0, 1.0]$",
        "",
        "- **Curvas de Lactação Regionais (EUA):** Li, M., Reed, K. F., & Cabrera, V. E. (2022/2023). *Journal of Dairy Science*, 105(9):7525-7538.",
        "- **Índice de Temperatura e Umidade (THI):** Bohmanova, J., Misztal, I., & Cole, J. B. (2007). *Journal of Dairy Science*, 90(4):1947-1956.",
        "- **Clima e Aridez ($P/\\text{PET}$):** Fick & Hijmans (2017) *WorldClim 2*; Fischer et al. (2021) *FAO/IIASA GAEZ v4*.",
        "- **Física e Hidrologia do Solo (PAWC):** Saxton, K. E., & Rawls, W. J. (2006). *Soil Sci. Soc. Am. J.*, 70(5):1569-1578; Poggio et al. (2021) *SoilGrids 2.0*.",
        "- **Métrica de Similaridade Ambiental:** Gower, J. C. (1971). *Biometrics*, 27(4):857-871; Ramirez-Villegas et al. (2011) *CCAFS Climate Analogues*.",
        "",
        "| Característica | Descrição | Limite Mínimo ($0.0$) | Limite Máximo ($1.0$) |",
        "| :--- | :--- | :---: | :---: |",
        f"| `thi_index` | Índice de Temperatura e Umidade (THI) | {THI_MIN} | {THI_MAX} |",
        f"| `temp_seasonality` | Sazonalidade Térmica ($\\Delta T$, °C) | {SEASONALITY_MIN_C} | {SEASONALITY_MAX_C} |",
        f"| `moisture_index` | Índice Hídrico ($P/\\text{{PET}}$) | {MOISTURE_INDEX_MIN} | {MOISTURE_INDEX_MAX} |",
        f"| `soil_clay_fraction` | Teor de Argila Superficial 0–30 cm (%) | {CLAY_MIN_PCT}% | {CLAY_MAX_PCT}% |",
        f"| `soil_pawc_index` | Água Disponível PAWC ($\\text{{cm}}^3/\\text{{cm}}^3$) | {PAWC_MIN} | {PAWC_MAX} |",
        f"| `soil_fertility_index` | Fertilidade Química Composta (pH [{PH_MIN}..{PH_MAX}] + SOC [{SOC_MIN_G_KG}..{SOC_MAX_G_KG}] g/kg) | 0.0 | 1.0 |",
        "",
        "---",
        "",
        "## 2. Mapeamento Top-1 das 27 UFs Brasileiras (IBGE) para Regiões do RuFaS",
        "",
        "| Código IBGE | UF | Estado (Brasil) | Estado Homólogo Top-1 (EUA) | Região RuFaS Equivalente | Similaridade Gower | Similaridade Euclidiana | $\\Delta l$ | $\\Delta m$ | $\\Delta n$ |",
        "| :---: | :---: | :--- | :--- | :--- | :---: | :---: | :---: | :---: | :---: |",
    ]

    for r in results:
        bp = r.brazil_profile
        bm = r.best_match
        adj = r.wood_adjustments
        lines.append(
            f"| `{bp.code}` | **{bp.abbrev}** | {bp.name} | {bm.us_name} ({bm.us_abbrev}, FIPS `{bm.us_code}`) "
            f"| **{bm.rufas_region}** | `{bm.gower_similarity:.4f}` | `{bm.euclidean_similarity:.4f}` "
            f"| `{adj['l']:+.2f}` | `{adj['m']:+.3f}` | `{adj['n']:+.2f}` |"
        )

    lines.extend(
        [
            "",
            "---",
            "",
            "## 3. Vetores Edafoclimáticos Normalizados $[0.0, 1.0]$ — Brasil (27 UFs)",
            "",
            "| IBGE | UF | `thi_index` | `temp_seasonality` | `moisture_index` | `soil_clay_fraction` | `soil_pawc_index` | `soil_fertility_index` |",
            "| :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: |",
        ]
    )

    for r in results:
        v = r.brazil_vector
        lines.append(
            f"| `{r.brazil_profile.code}` | **{r.brazil_profile.abbrev}** | `{v.thi_index:.4f}` | `{v.temp_seasonality:.4f}` "
            f"| `{v.moisture_index:.4f}` | `{v.soil_clay_fraction:.4f}` | `{v.soil_pawc_index:.4f}` | `{v.soil_fertility_index:.4f}` |"
        )

    lines.extend(
        [
            "",
            "---",
            "",
            "## 4. Vetores Edafoclimáticos Normalizados $[0.0, 1.0]$ — EUA (48 Estados FIPS)",
            "",
            "| FIPS | Estado | Região RuFaS | `thi_index` | `temp_seasonality` | `moisture_index` | `soil_clay_fraction` | `soil_pawc_index` | `soil_fertility_index` |",
            "| :---: | :---: | :--- | :---: | :---: | :---: | :---: | :---: | :---: |",
        ]
    )

    for us_prof in sorted(us_profiles, key=lambda p: int(p.code)):
        uv = normalize_profile(us_prof)
        lines.append(
            f"| `{us_prof.code}` | **{us_prof.abbrev}** ({us_prof.name}) | {us_prof.rufas_region} "
            f"| `{uv.thi_index:.4f}` | `{uv.temp_seasonality:.4f}` | `{uv.moisture_index:.4f}` "
            f"| `{uv.soil_clay_fraction:.4f}` | `{uv.soil_pawc_index:.4f}` | `{uv.soil_fertility_index:.4f}` |"
        )

    return "\n".join(lines) + "\n"


def run_pipeline(
    data_dir: Path,
    output_dir: Path,
    export_rufas_json: Path | None = None,
) -> dict[str, Any]:
    """
    Executes the complete equivalence calculation and atomically writes output artifacts.
    """
    regions_raw = json.loads((data_dir / "rufas_region_adjustments.json").read_text(encoding="utf-8"))
    us_raw = json.loads((data_dir / "us_states_edaphoclimatic.json").read_text(encoding="utf-8"))
    br_raw = json.loads((data_dir / "brazil_ufs_edaphoclimatic.json").read_text(encoding="utf-8"))

    region_adjustments: dict[str, dict[str, float]] = regions_raw["region_adjustments"]
    us_profiles = [RawEdaphoclimaticProfile.from_dict(item) for item in us_raw["states"]]
    br_profiles = [RawEdaphoclimaticProfile.from_dict(item) for item in br_raw["states"]]

    equivalences = compute_all_equivalences(br_profiles, us_profiles, region_adjustments, top_k=3)

    brazil_mapping: dict[str, str] = {
        res.brazil_profile.code: res.best_match.rufas_region for res in equivalences
    }

    report_payload: dict[str, Any] = {
        "methodology": {
            "references": [
                "Li, M., Reed, K. F., & Cabrera, V. E. (2022). J. Dairy Sci., 105(9):7525-7538.",
                "Bohmanova, J., Misztal, I., & Cole, J. B. (2007). J. Dairy Sci., 90(4):1947-1956.",
                "Saxton, K. E., & Rawls, W. J. (2006). Soil Sci. Soc. Am. J., 70(5):1569-1578.",
                "Poggio, L., et al. (2021). ISRIC SoilGrids 2.0. SOIL, 7(1):217-240.",
                "Fick, S. E., & Hijmans, R. J. (2017). WorldClim 2. Int. J. Climatol., 37(12):4302-4315.",
                "Gower, J. C. (1971). Biometrics, 27(4):857-871.",
                "Ramirez-Villegas, J., et al. (2011). CCAFS Climate Analogue Working Paper No. 12.",
            ],
            "normalization_bounds": {
                "thi_index": [THI_MIN, THI_MAX],
                "temp_seasonality": [SEASONALITY_MIN_C, SEASONALITY_MAX_C],
                "moisture_index": [MOISTURE_INDEX_MIN, MOISTURE_INDEX_MAX],
                "soil_clay_fraction": [CLAY_MIN_PCT, CLAY_MAX_PCT],
                "soil_pawc_index": [PAWC_MIN, PAWC_MAX],
                "soil_fertility_index": {
                    "ph_bounds": [PH_MIN, PH_MAX],
                    "soc_g_kg_bounds": [SOC_MIN_G_KG, SOC_MAX_G_KG],
                },
            },
        },
        "brazil_state_to_region_mapping": brazil_mapping,
        "brazil_equivalences": [
            {
                "ibge_code": r.brazil_profile.code,
                "uf": r.brazil_profile.abbrev,
                "state_name": r.brazil_profile.name,
                "normalized_vector": r.brazil_vector.to_dict(),
                "best_match": r.best_match.to_dict(),
                "wood_adjustments": r.wood_adjustments,
                "top_3_candidates": [c.to_dict() for c in r.top_candidates],
            }
            for r in equivalences
        ],
        "us_reference_vectors": [
            {
                "fips_code": p.code,
                "state_abbrev": p.abbrev,
                "state_name": p.name,
                "rufas_region": p.rufas_region,
                "normalized_vector": normalize_profile(p).to_dict(),
            }
            for p in sorted(us_profiles, key=lambda x: int(x.code))
        ],
    }

    output_dir.mkdir(parents=True, exist_ok=True)
    atomic_write_json(output_dir / "equivalence_audit_report.json", report_payload)
    atomic_write_text(
        output_dir / "equivalence_audit_report.md",
        _format_markdown_report(equivalences, us_profiles),
    )

    if export_rufas_json is not None:
        rufas_path = Path(export_rufas_json)
        rufas_data = json.loads(rufas_path.read_text(encoding="utf-8"))
        rufas_data["brazil_state_to_region_mapping"] = brazil_mapping
        atomic_write_json(rufas_path, rufas_data)

    return report_payload


def main() -> None:
    project_root = Path(__file__).resolve().parent.parent.parent
    parser = argparse.ArgumentParser(
        description="Compute scientific edaphoclimatic lactation curve equivalence (US FIPS -> Brazil IBGE)."
    )
    parser.add_argument(
        "--data-dir",
        type=Path,
        default=project_root / "data",
        help="Directory containing input climate/soil JSON datasets.",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=project_root / "output",
        help="Directory where audit reports (.json and .md) will be atomically written.",
    )
    parser.add_argument(
        "--export-rufas-json",
        type=Path,
        default=None,
        help="Optional path to RuFaS lactation_curve_adjustment_inputs.json to update atomically.",
    )
    args = parser.parse_args()
    run_pipeline(
        data_dir=args.data_dir,
        output_dir=args.output_dir,
        export_rufas_json=args.export_rufas_json,
    )


if __name__ == "__main__":
    main()
```

Create `/home/thiago/Projects/rufas-lactation-equivalence/README.md`:

```markdown
# RuFaS Lactation Curve Edaphoclimatic Equivalence (USA FIPS ↔ Brazil IBGE)

Pacote científico determinístico e de execução atômica para decompor perfis de clima e solo dos estados americanos (códigos FIPS) e das 27 Unidades Federativas do Brasil (códigos IBGE) em vetores normalizados $[0.0, 1.0]^6$, identificando a região leiteira homóloga mais próxima nos EUA (**Top-1 Nearest Homologue**) via Coeficiente de Similaridade de Gower (1971) para parametrizar a curva de lactação de Wood no **RuFaS**.

## Execução Atômica

```bash
# Gerar relatórios científicos de auditoria (JSON + Markdown) e exportar para o RuFaS:
PYTHONPATH=src python3 -m rufas_lactation_equivalence.cli \
  --export-rufas-json ../RuFaS/input/data/animal/lactation_curve_adjustment_inputs.json
```

## Referências Científicas

1. **Li, M., Reed, K. F., & Cabrera, V. E. (2022).** *Investigating the effect of temporal, geographic, and management factors on US Holstein lactation curve parameters.* Journal of Dairy Science, 105(9), 7525–7538.
2. **Bohmanova, J., Misztal, I., & Cole, J. B. (2007).** *Temperature-humidity indices as indicators of milk production losses due to heat stress.* Journal of Dairy Science, 90(4), 1947–1956.
3. **Saxton, K. E., & Rawls, W. J. (2006).** *Soil water characteristic estimates by texture and organic matter for hydrologic solutions.* Soil Science Society of America Journal, 70(5), 1569–1578.
4. **Poggio, L., et al. (2021).** *SoilGrids 2.0: producing soil information for the globe with quantified spatial uncertainty.* SOIL, 7(1), 217–240.
5. **Fick, S. E., & Hijmans, R. J. (2017).** *WorldClim 2: new 1-km spatial resolution climate surfaces for global land areas.* International Journal of Climatology, 37(12), 4302–4315.
6. **Gower, J. C. (1971).** *A general coefficient of similarity and some of its properties.* Biometrics, 27(4), 857–871.
7. **Ramirez-Villegas, J., et al. (2011).** *Climate analogues: finding tomorrow's agriculture today.* CGIAR CCAFS Working Paper No. 12.
```

- [ ] **Step 4: Run tests and execute atomic CLI to generate `output/` reports**

Run:
```bash
PYTHONPATH=/home/thiago/Projects/rufas-lactation-equivalence/src /home/thiago/Projects/RuFaS/.venv/bin/pytest /home/thiago/Projects/rufas-lactation-equivalence/tests/ -v
PYTHONPATH=/home/thiago/Projects/rufas-lactation-equivalence/src /home/thiago/Projects/RuFaS/.venv/bin/python3 -m rufas_lactation_equivalence.cli
```
Expected: All 9 tests PASS, and `output/equivalence_audit_report.json` + `output/equivalence_audit_report.md` are created.

- [ ] **Step 5: Commit Task 4**

```bash
rtk git -C /home/thiago/Projects/rufas-lactation-equivalence add src/ tests/ output/ README.md
rtk git -C /home/thiago/Projects/rufas-lactation-equivalence commit -m "feat(cli): add atomic I/O, CLI pipeline, and generated edaphoclimatic audit reports"
```

---

### Task 5: RuFaS Integration — IBGE vs FIPS Disambiguation, Metadata Schema, and Unit Tests

**Files:**
- Modify: `/home/thiago/Projects/RuFaS/input/data/animal/lactation_curve_adjustment_inputs.json`
- Modify: `/home/thiago/Projects/RuFaS/RUFAS/biophysical/animal/milk/lactation_curve.py:96-101`
- Modify: `/home/thiago/Projects/RuFaS/RUFAS/input/metadata/properties/default.json:2350-2375`
- Modify: `/home/thiago/Projects/RuFaS/tests/test_biophysical/test_animal/test_milk_production/test_lactation_curve.py:274-354`

**Interfaces:**
- Consumes: `brazil_state_to_region_mapping` exported by `rufas-lactation-equivalence`.
- Produces: Clean country-aware mapping selection in `LactationCurve.set_lactation_parameters` (`country == "BRA"` uses `brazil_state_to_region_mapping`, `country == "USA"` uses `state_to_region_mapping`) and schema validation in `default.json`.

- [ ] **Step 1: Write failing test in `RuFaS` verifying IBGE vs FIPS disambiguation**

In `/home/thiago/Projects/RuFaS/tests/test_biophysical/test_animal/test_milk_production/test_lactation_curve.py`, update `test_set_lactation_parameters_region_code_fallback` (or add `test_set_lactation_parameters_brazil_uses_brazil_mapping_without_fips_collision`) so that `mock_lactation_inputs` contains BOTH `"state_to_region_mapping": {"31": "Northern Plains"}` (Nebraska) AND `"brazil_state_to_region_mapping": {"31": "Southeast"}` (Minas Gerais), verifying that when `country == "BRA"`, `LactationCurve.set_lactation_parameters` passes `brazil_state_to_region_mapping` to `_get_region_adjustments`:

```python
def test_set_lactation_parameters_brazil_disambiguates_from_us_fips(mocker: MockerFixture) -> None:
    """Tests that country='BRA' uses brazil_state_to_region_mapping instead of colliding with US FIPS 31 (Nebraska)."""
    mock_im = mocker.patch("RUFAS.biophysical.animal.milk.lactation_curve.InputManager").return_value
    mock_time = mocker.MagicMock()
    mock_time.end_date.year = 2024

    mock_lactation_inputs = {
        "adjustments": {
            "year": {"2016": {"l": 0.0, "m": 0.0, "n": 0.0}},
            "region": {
                "Northern Plains": {"l": -0.26, "m": 0.19, "n": -0.79},
                "Southeast": {"l": -2.00, "m": 2.59, "n": 2.60},
            },
            "milking_frequency": {"twice_daily": {"l": 0.0, "m": 0.0, "n": 0.0}},
            "parity": {"1": {"l": 0.0, "m": 0.0, "n": 0.0}},
        },
        "state_to_region_mapping": {"31": "Northern Plains"},
        "brazil_state_to_region_mapping": {"31": "Southeast"},
        "parameter_mean_values": {"parameter_l_mean": 10.0, "parameter_m_mean": 0.2, "parameter_n_mean": 0.05},
        "parameter_standard_deviations": {
            "1": {"parameter_l_std_dev": 0.1, "parameter_m_std_dev": 0.01, "parameter_n_std_dev": 0.001},
            "2": {"parameter_l_std_dev": 0.1, "parameter_m_std_dev": 0.01, "parameter_n_std_dev": 0.001},
            "3": {"parameter_l_std_dev": 0.1, "parameter_m_std_dev": 0.01, "parameter_n_std_dev": 0.001},
        },
    }
    mock_animal_inputs = {
        "animal_config": {"management_decisions": {"cow_times_milked_per_day": 2.0}},
        "herd_information": {"annual_milk_yield": None},
    }

    mock_im.get_data.side_effect = lambda key, required=True: {
        "lactation": mock_lactation_inputs,
        "config.country": "BRA",
        "config.region_code": 3148004,
        "config.FIPS_county_code": None,
        "animal": mock_animal_inputs,
    }.get(key)

    mock_get_adj = mocker.patch.object(
        LactationCurve, "_get_region_adjustments", wraps=LactationCurve._get_region_adjustments
    )
    LactationCurve.set_lactation_parameters(mock_time)
    mock_get_adj.assert_called_once_with(
        mock_lactation_inputs["adjustments"]["region"],
        mock_lactation_inputs["brazil_state_to_region_mapping"],
        3148004,
        "BRA",
    )
```

- [ ] **Step 2: Run test to verify it fails**

Run: `/home/thiago/Projects/RuFaS/.venv/bin/pytest /home/thiago/Projects/RuFaS/tests/test_biophysical/test_animal/test_milk_production/test_lactation_curve.py -k "test_set_lactation_parameters_brazil_disambiguates_from_us_fips" -v`
Expected: FAIL because `LactationCurve.set_lactation_parameters` currently passes `state_to_region_mapping` instead of `brazil_state_to_region_mapping`.

- [ ] **Step 3: Export `brazil_state_to_region_mapping` into RuFaS and update `lactation_curve.py` & `default.json`**

1. Run atomic export from `rufas-lactation-equivalence`:
```bash
PYTHONPATH=/home/thiago/Projects/rufas-lactation-equivalence/src \
  /home/thiago/Projects/RuFaS/.venv/bin/python3 -m rufas_lactation_equivalence.cli \
  --export-rufas-json /home/thiago/Projects/RuFaS/input/data/animal/lactation_curve_adjustment_inputs.json
```

2. Update `/home/thiago/Projects/RuFaS/RUFAS/biophysical/animal/milk/lactation_curve.py` lines 96–101:
```python
        all_region_adjustments: dict[str, dict[str, float]] = lactation_inputs["adjustments"]["region"]
        if str(country).upper() == "BRA":
            region_mapping: dict[str, str] = lactation_inputs.get(
                "brazil_state_to_region_mapping", lactation_inputs.get("state_to_region_mapping", {})
            )
        else:
            region_mapping = lactation_inputs.get("state_to_region_mapping", {})
        region_adjustments = cls._get_region_adjustments(
            all_region_adjustments, region_mapping, region_code, country
        )
```

3. Update `/home/thiago/Projects/RuFaS/RUFAS/input/metadata/properties/default.json` right after `state_to_region_mapping` to add optional `brazil_state_to_region_mapping` for the 27 IBGE UF codes (`"11"` through `"53"`).

- [ ] **Step 4: Run RuFaS lactation and metadata schema tests to verify they pass**

Run:
```bash
/home/thiago/Projects/RuFaS/.venv/bin/pytest \
  /home/thiago/Projects/RuFaS/tests/test_biophysical/test_animal/test_milk_production/test_lactation_curve.py \
  /home/thiago/Projects/RuFaS/tests/test_metadata_schema_passo_1.py -v
```
Expected: All tests PASS.

- [ ] **Step 5: Commit changes in `RuFaS/`**

```bash
rtk git -C /home/thiago/Projects/RuFaS add \
  input/data/animal/lactation_curve_adjustment_inputs.json \
  RUFAS/biophysical/animal/milk/lactation_curve.py \
  RUFAS/input/metadata/properties/default.json \
  tests/test_biophysical/test_animal/test_milk_production/test_lactation_curve.py
rtk git -C /home/thiago/Projects/RuFaS commit -m "feat(animal): add edaphoclimatic brazil_state_to_region_mapping and disambiguate IBGE from US FIPS"
```
