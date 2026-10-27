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

    def __post_init__(self) -> None:
        if not (0.0 <= self.rh_mean_pct <= 100.0):
            raise ValueError(
                f"rh_mean_pct must be between 0.0 and 100.0, got {self.rh_mean_pct}"
            )
        if not (0.0 <= self.soil_clay_pct <= 100.0):
            raise ValueError(
                f"soil_clay_pct must be between 0.0 and 100.0, got {self.soil_clay_pct}"
            )
        if not (0.0 <= self.soil_sand_pct <= 100.0):
            raise ValueError(
                f"soil_sand_pct must be between 0.0 and 100.0, got {self.soil_sand_pct}"
            )
        if self.precip_annual_mm < 0.0:
            raise ValueError(
                f"precip_annual_mm must be non-negative, got {self.precip_annual_mm}"
            )
        if self.pet_annual_mm < 0.0:
            raise ValueError(
                f"pet_annual_mm must be non-negative, got {self.pet_annual_mm}"
            )
        if self.soil_ph <= 0.0:
            raise ValueError(f"soil_ph must be positive, got {self.soil_ph}")
        if self.soil_soc_g_kg < 0.0:
            raise ValueError(
                f"soil_soc_g_kg must be non-negative, got {self.soil_soc_g_kg}"
            )

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

    def __post_init__(self) -> None:
        fields = (
            ("thi_index", self.thi_index),
            ("temp_seasonality", self.temp_seasonality),
            ("moisture_index", self.moisture_index),
            ("soil_clay_fraction", self.soil_clay_fraction),
            ("soil_pawc_index", self.soil_pawc_index),
            ("soil_fertility_index", self.soil_fertility_index),
        )
        for name, val in fields:
            if not (0.0 <= val <= 1.0):
                raise ValueError(f"{name} must be within [0.0, 1.0], got {val}")

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

    def __post_init__(self) -> None:
        if not (0.0 <= self.gower_similarity <= 1.0):
            raise ValueError(
                f"gower_similarity must be within [0.0, 1.0], got {self.gower_similarity}"
            )
        if not (0.0 <= self.euclidean_similarity <= 1.0):
            raise ValueError(
                f"euclidean_similarity must be within [0.0, 1.0], got {self.euclidean_similarity}"
            )

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
