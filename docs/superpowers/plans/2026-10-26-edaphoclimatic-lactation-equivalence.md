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

### Task 3: Scientific API Data Fetcher (NASA POWER + ISRIC SoilGrids v2.0 WCS) and Reference Coordinates

**Files:**
- Create: `/home/thiago/Projects/rufas-lactation-equivalence/data/reference_coordinates.json`
- Create: `/home/thiago/Projects/rufas-lactation-equivalence/data/rufas_region_adjustments.json`
- Create: `/home/thiago/Projects/rufas-lactation-equivalence/src/rufas_lactation_equivalence/fetcher.py`
- Test: `/home/thiago/Projects/rufas-lactation-equivalence/tests/test_fetcher.py`
- Test: `/home/thiago/Projects/rufas-lactation-equivalence/tests/test_datasets.py`

**Interfaces:**
- Consumes: `data/reference_coordinates.json` (coordinates `(lat, lon)` for 48 US FIPS states and 27 Brazilian IBGE UFs).
- Produces:
  - `decode_soilgrids_tiff(tiff_bytes: bytes) -> float` (pure stdlib `zlib` + `struct` GeoTIFF decoder for ISRIC WCS 2.0.1)
  - `fetch_nasa_power_climate(lat: float, lon: float) -> dict[str, float]`
  - `fetch_isric_soilgrids(lat: float, lon: float) -> dict[str, float]`
  - `fetch_all_profiles(coords_file: Path, output_us_file: Path, output_br_file: Path) -> None`
  - Atomically generated `data/us_states_edaphoclimatic.json` and `data/brazil_ufs_edaphoclimatic.json` with full API provenance metadata.

- [ ] **Step 1: Write failing tests for GeoTIFF decoding, NASA POWER parsing, and dataset completeness**

Create `/home/thiago/Projects/rufas-lactation-equivalence/tests/test_fetcher.py`:

```python
import pytest

from rufas_lactation_equivalence.fetcher import (
    compute_thornthwaite_or_hargreaves_pet,
    parse_nasa_power_response,
    weighted_topsoil_average,
)


def test_parse_nasa_power_response_extracts_climate_metrics() -> None:
    """Verify extraction of T_mean, RH_mean, seasonality (max-min), precip_annual, and PET."""
    sample_payload = {
        "properties": {
            "parameter": {
                "T2M": {
                    "JAN": 23.0, "FEB": 23.2, "MAR": 22.8, "APR": 21.5,
                    "MAY": 19.0, "JUN": 17.5, "JUL": 17.2, "AUG": 19.5,
                    "SEP": 21.8, "OCT": 22.9, "NOV": 22.7, "DEC": 22.8,
                    "ANN": 21.16,
                },
                "RH2M": {"ANN": 70.5},
                "PRECTOTCORR": {"ANN": 4.0},
            }
        }
    }
    parsed = parse_nasa_power_response(sample_payload, lat=-18.5)
    assert pytest.approx(parsed["temp_mean_c"], rel=1e-3) == 21.16
    assert pytest.approx(parsed["rh_mean_pct"], rel=1e-3) == 70.5
    assert pytest.approx(parsed["temp_seasonality_c"], rel=1e-3) == (23.2 - 17.2)
    assert pytest.approx(parsed["precip_annual_mm"], rel=1e-3) == 4.0 * 365.0
    assert 800.0 <= parsed["pet_annual_mm"] <= 1800.0


def test_weighted_topsoil_average_0_to_30cm() -> None:
    """Verify 0-30 cm depth-weighted average (5cm * d1 + 10cm * d2 + 15cm * d3) / 30."""
    avg = weighted_topsoil_average(val_0_5=300.0, val_5_15=360.0, val_15_30=420.0)
    expected = (5.0 * 300.0 + 10.0 * 360.0 + 15.0 * 420.0) / 30.0
    assert pytest.approx(avg) == expected
```

- [ ] **Step 2: Run test to verify it fails**

Run: `PYTHONPATH=/home/thiago/Projects/rufas-lactation-equivalence/src /home/thiago/Projects/RuFaS/.venv/bin/pytest /home/thiago/Projects/rufas-lactation-equivalence/tests/test_fetcher.py -v`
Expected: FAIL with `ModuleNotFoundError: No module named 'rufas_lactation_equivalence.fetcher'`

- [ ] **Step 3: Implement `fetcher.py` and `data/reference_coordinates.json`, then execute live fetch to populate `data/`**

Create `/home/thiago/Projects/rufas-lactation-equivalence/src/rufas_lactation_equivalence/fetcher.py` querying:
- `https://power.larc.nasa.gov/api/temporal/climatology/point?parameters=T2M,RH2M,PRECTOTCORR&community=AG&longitude={lon}&latitude={lat}&format=JSON`
- `https://maps.isric.org/mapserv?map=/map/{prop}.map&SERVICE=WCS&VERSION=2.0.1&REQUEST=GetCoverage&COVERAGEID={prop}_{depth}_mean&FORMAT=image/tiff&SUBSET=X({lon_min},{lon_max})&SUBSET=Y({lat_min},{lat_max})&SUBSETTINGCRS=http://www.opengis.net/def/crs/EPSG/0/4326&OUTPUTCRS=http://www.opengis.net/def/crs/EPSG/0/4326`

- [ ] **Step 4: Run fetcher & dataset validation tests**

Run: `PYTHONPATH=/home/thiago/Projects/rufas-lactation-equivalence/src /home/thiago/Projects/RuFaS/.venv/bin/pytest /home/thiago/Projects/rufas-lactation-equivalence/tests/test_fetcher.py /home/thiago/Projects/rufas-lactation-equivalence/tests/test_datasets.py -v`
Expected: PASS

- [ ] **Step 5: Commit Task 3**

```bash
rtk git -C /home/thiago/Projects/rufas-lactation-equivalence add data/ src/rufas_lactation_equivalence/fetcher.py tests/test_fetcher.py tests/test_datasets.py
rtk git -C /home/thiago/Projects/rufas-lactation-equivalence commit -m "feat(fetcher): add NASA POWER and ISRIC SoilGrids v2.0 WCS scientific data fetcher and datasets"
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
