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
    build_dashboard: Path | None = None,
) -> dict[str, Any]:
    """
    Executes the complete equivalence calculation and atomically writes output artifacts.
    """
    data_dir = Path(data_dir)
    output_dir = Path(output_dir)
    regions_raw = json.loads((data_dir / "rufas_region_adjustments.json").read_text(encoding="utf-8"))
    us_raw = json.loads((data_dir / "us_states_edaphoclimatic.json").read_text(encoding="utf-8"))
    br_raw = json.loads((data_dir / "brazil_ufs_edaphoclimatic.json").read_text(encoding="utf-8"))

    region_adjustments: dict[str, dict[str, float]] = regions_raw["region_adjustments"]
    us_items = us_raw.get("profiles", us_raw.get("states", []))
    br_items = br_raw.get("profiles", br_raw.get("states", []))
    us_profiles = [RawEdaphoclimaticProfile.from_dict(item) for item in us_items]
    br_profiles = [RawEdaphoclimaticProfile.from_dict(item) for item in br_items]

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

    if build_dashboard is not None:
        from rufas_lactation_equivalence.dashboard_builder import build_and_save_dashboard

        build_and_save_dashboard(
            destination_path=build_dashboard,
            audit_report_path=output_dir / "equivalence_audit_report.json",
            edaphoclimatic_path=data_dir / "brazil_ufs_edaphoclimatic.json",
            geojson_path=data_dir / "brazil_ibge_ufs.geojson",
        )

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
    parser.add_argument(
        "--build-dashboard",
        type=Path,
        nargs="?",
        const=project_root / "docs" / "index.html",
        default=None,
        help="Export standalone interactive HTML dashboard.",
    )
    args = parser.parse_args()
    run_pipeline(
        data_dir=args.data_dir,
        output_dir=args.output_dir,
        export_rufas_json=args.export_rufas_json,
        build_dashboard=args.build_dashboard,
    )


if __name__ == "__main__":
    main()
