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
