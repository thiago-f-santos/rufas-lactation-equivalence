"""
Dashboard Data Preprocessor and IBGE SVG Map Projection Generator.

Transforms biophysical edaphoclimatic profiles, Gower similarity matches,
and Wood lactation curve adjustment parameters into an integrated,
deterministic JSON payload with projected SVG vector geometries.
"""

from __future__ import annotations

import json
import math
from pathlib import Path
from typing import Any

BASE_WOOD_L: float = 19.90
BASE_WOOD_M: float = 0.247
BASE_WOOD_N: float = 0.003376

DEFAULT_LON_MIN: float = -74.0
DEFAULT_LON_MAX: float = -34.5
DEFAULT_LAT_MIN: float = -33.8
DEFAULT_LAT_MAX: float = 5.5

DEFAULT_SVG_WIDTH: float = 600.0
DEFAULT_SVG_HEIGHT: float = 600.0

REGION_COLORS: dict[str, str] = {
    "Delta": "#10b981",
    "Southeast": "#3b82f6",
    "Mountain": "#f59e0b",
    "Southern Plains": "#ef4444",
    "Appalachian": "#8b5cf6",
}


def project_coordinates_to_svg(
    lon: float,
    lat: float,
    width: float = DEFAULT_SVG_WIDTH,
    height: float = DEFAULT_SVG_HEIGHT,
    lon_min: float = DEFAULT_LON_MIN,
    lon_max: float = DEFAULT_LON_MAX,
    lat_min: float = DEFAULT_LAT_MIN,
    lat_max: float = DEFAULT_LAT_MAX,
) -> tuple[float, float]:
    """
    Projects geographic coordinates (lon, lat) to an SVG canvas (x, y)
    using Mercator projection bounded to the specified bounding box.

    Parameters
    ----------
    lon : float
        Longitude in decimal degrees.
    lat : float
        Latitude in decimal degrees.
    width : float
        Target SVG canvas width in user units.
    height : float
        Target SVG canvas height in user units.
    lon_min, lon_max : float
        Longitude bounds.
    lat_min, lat_max : float
        Latitude bounds.

    Returns
    -------
    tuple[float, float]
        Projected (x, y) coordinates rounded to 2 decimal places.
    """
    rad_lat = math.radians(lat)
    y_merc = math.log(math.tan(math.pi / 4.0 + rad_lat / 2.0))
    y_merc_min = math.log(math.tan(math.pi / 4.0 + math.radians(lat_min) / 2.0))
    y_merc_max = math.log(math.tan(math.pi / 4.0 + math.radians(lat_max) / 2.0))

    x = (lon - lon_min) / (lon_max - lon_min) * width
    y = (y_merc_max - y_merc) / (y_merc_max - y_merc_min) * height

    return round(x, 2), round(y, 2)


def convert_geometry_to_svg_path(
    geometry: dict[str, Any],
    width: float = DEFAULT_SVG_WIDTH,
    height: float = DEFAULT_SVG_HEIGHT,
    lon_min: float = DEFAULT_LON_MIN,
    lon_max: float = DEFAULT_LON_MAX,
    lat_min: float = DEFAULT_LAT_MIN,
    lat_max: float = DEFAULT_LAT_MAX,
) -> str:
    """
    Converts a GeoJSON Polygon or MultiPolygon geometry into an SVG path 'd' string.

    Parameters
    ----------
    geometry : dict[str, Any]
        GeoJSON geometry object with 'type' and 'coordinates'.
    width, height : float
        Target SVG canvas dimensions.

    Returns
    -------
    str
        SVG path data commands (M x y L x y ... Z).
    """
    geom_type = geometry.get("type")
    coords = geometry.get("coordinates", [])

    if geom_type == "Polygon":
        polygons = [coords]
    elif geom_type == "MultiPolygon":
        polygons = coords
    else:
        return ""

    path_segments: list[str] = []
    for poly in polygons:
        for ring in poly:
            if not ring:
                continue
            first_pt = ring[0]
            x0, y0 = project_coordinates_to_svg(
                first_pt[0], first_pt[1], width, height, lon_min, lon_max, lat_min, lat_max
            )
            ring_cmds = [f"M {x0} {y0}"]
            for pt in ring[1:]:
                x, y = project_coordinates_to_svg(
                    pt[0], pt[1], width, height, lon_min, lon_max, lat_min, lat_max
                )
                ring_cmds.append(f"L {x} {y}")
            ring_cmds.append("Z")
            path_segments.append(" ".join(ring_cmds))

    return " ".join(path_segments)


def compute_wood_curve(
    delta_l: float = 0.0,
    delta_m: float = 0.0,
    delta_n: float = 0.0,
) -> dict[str, Any]:
    """
    Calculates daily lactation yield Y(t) = l * t^m * exp(-n * t) for DIM 1..305,
    along with peak DIM, peak daily yield, and cumulative 305-day yield.

    Scaling factors:
        l_adj = l_base + delta_l
        m_adj = m_base + (delta_m * 1e-2)
        n_adj = n_base + (delta_n * 1e-4)

    Returns
    -------
    dict[str, Any]
        Dictionary containing 'dim', 'yield_kg', 'peak_day', 'peak_yield', 'total_305d', and 'params'.
    """
    l_adj = BASE_WOOD_L + delta_l
    m_adj = BASE_WOOD_M + (delta_m * 1e-2)
    n_adj = BASE_WOOD_N + (delta_n * 1e-4)

    t_peak = m_adj / n_adj if n_adj > 0 else 0.0
    y_peak = (
        l_adj * (t_peak**m_adj) * math.exp(-m_adj)
        if t_peak > 0
        else 0.0
    )

    dim = list(range(1, 306))
    yield_kg = [
        round(l_adj * (t**m_adj) * math.exp(-n_adj * t), 4)
        for t in dim
    ]
    total_305d = round(sum(yield_kg), 2)

    return {
        "dim": dim,
        "yield_kg": yield_kg,
        "peak_day": round(t_peak, 2),
        "peak_yield": round(y_peak, 2),
        "total_305d": total_305d,
        "params": {
            "l_adj": round(l_adj, 4),
            "m_adj": round(m_adj, 6),
            "n_adj": round(n_adj, 8),
            "delta_l": delta_l,
            "delta_m": delta_m,
            "delta_n": delta_n,
        },
    }


def _resolve_project_path(given: Path | str | None, relative_subpath: str) -> Path:
    if given is not None:
        p = Path(given)
        if p.exists():
            return p

    # Try relative to repo root
    repo_root = Path(__file__).resolve().parent.parent.parent
    candidate = repo_root / relative_subpath
    if candidate.exists():
        return candidate

    # Try current working directory
    cwd_candidate = Path.cwd() / relative_subpath
    if cwd_candidate.exists():
        return cwd_candidate

    raise FileNotFoundError(f"Could not locate required file: {relative_subpath}")


def generate_dashboard_payload(
    audit_report_path: Path | str | None = None,
    edaphoclimatic_path: Path | str | None = None,
    geojson_path: Path | str | None = None,
    svg_width: float = DEFAULT_SVG_WIDTH,
    svg_height: float = DEFAULT_SVG_HEIGHT,
) -> dict[str, Any]:
    """
    Constructs the comprehensive data payload powering the lactation equivalence dashboard.

    Merges:
    - 27 Brazilian UFs edaphoclimatic baseline profiles
    - IBGE official simplified SVG polygon geometries
    - Audit equivalence rankings, Gower similarity, and US state homologues
    - Calibrated Wood's lactation curve parameters and 305-day yields

    Returns
    -------
    dict[str, Any]
        Dictionary with 'summary' and 'ufs' mappings.
    """
    audit_p = _resolve_project_path(audit_report_path, "output/equivalence_audit_report.json")
    edap_p = _resolve_project_path(edaphoclimatic_path, "data/brazil_ufs_edaphoclimatic.json")
    geo_p = _resolve_project_path(geojson_path, "data/brazil_ibge_ufs.geojson")

    with open(audit_p, "r", encoding="utf-8") as f:
        audit_data = json.load(f)

    with open(edap_p, "r", encoding="utf-8") as f:
        edap_data = json.load(f)

    with open(geo_p, "r", encoding="utf-8") as f:
        geojson_data = json.load(f)

    geo_lookup: dict[str, dict[str, Any]] = {
        feat["properties"]["codarea"]: feat["geometry"]
        for feat in geojson_data.get("features", [])
        if "properties" in feat and "codarea" in feat["properties"]
    }

    raw_lookup: dict[str, dict[str, Any]] = {
        str(item["code"]): item
        for item in edap_data.get("profiles", [])
    }

    us_vector_lookup: dict[str, dict[str, float]] = {
        item["state_abbrev"]: item["normalized_vector"]
        for item in audit_data.get("us_reference_vectors", [])
        if "state_abbrev" in item and "normalized_vector" in item
    }

    state_to_region = audit_data.get("brazil_state_to_region_mapping", {})

    ufs: dict[str, dict[str, Any]] = {}
    for item in audit_data.get("brazil_equivalences", []):
        uf = str(item["uf"])
        code = str(item["ibge_code"])
        state_name = str(item["state_name"])
        region = state_to_region.get(code) or item.get("best_match", {}).get("rufas_region", "Southeast")
        region_color = REGION_COLORS.get(region, "#3b82f6")

        best_match = dict(item.get("best_match", {}))
        us_abbrev = best_match.get("us_abbrev", "")
        if "normalized_vector" not in best_match and us_abbrev in us_vector_lookup:
            best_match["normalized_vector"] = us_vector_lookup[us_abbrev]

        wood_adj = item.get("wood_adjustments", {"l": 0.0, "m": 0.0, "n": 0.0})
        wood_curve = compute_wood_curve(
            delta_l=float(wood_adj.get("l", 0.0)),
            delta_m=float(wood_adj.get("m", 0.0)),
            delta_n=float(wood_adj.get("n", 0.0)),
        )

        geom = geo_lookup.get(code, {})
        svg_path = convert_geometry_to_svg_path(
            geom,
            width=svg_width,
            height=svg_height,
        )

        raw_prof = dict(raw_lookup.get(code, {}))
        if "soil_sand_pct" in raw_prof and "soil_clay_pct" in raw_prof and "soil_soc_g_kg" in raw_prof:
            from rufas_lactation_equivalence.biophysics import calculate_pawc

            raw_prof["pawc"] = round(
                calculate_pawc(
                    float(raw_prof["soil_sand_pct"]),
                    float(raw_prof["soil_clay_pct"]),
                    float(raw_prof["soil_soc_g_kg"]),
                ),
                4,
            )

        ufs[uf] = {
            "uf": uf,
            "name": state_name,
            "ibge_code": code,
            "rufas_region": region,
            "region_color": region_color,
            "raw_profile": raw_prof,
            "normalized_vector": item.get("normalized_vector", {}),
            "best_match": best_match,
            "top_3_candidates": item.get("top_3_candidates", []),
            "wood_adjustments": wood_adj,
            "wood_curve": wood_curve,
            "svg_path": svg_path,
        }

    gower_scores = [
        u["best_match"]["gower_similarity"]
        for u in ufs.values()
        if "best_match" in u and "gower_similarity" in u["best_match"]
    ]

    distinct_regions = sorted(list(set(u["rufas_region"] for u in ufs.values())))

    summary: dict[str, Any] = {
        "total_ufs": len(ufs),
        "total_regions": len(distinct_regions),
        "regions": distinct_regions,
        "region_colors": REGION_COLORS,
        "avg_similarity": round(sum(gower_scores) / max(1, len(gower_scores)), 4),
        "min_similarity": round(min(gower_scores), 4) if gower_scores else 0.0,
        "max_similarity": round(max(gower_scores), 4) if gower_scores else 0.0,
        "base_wood_parameters": {
            "l": BASE_WOOD_L,
            "m": BASE_WOOD_M,
            "n": BASE_WOOD_N,
        },
        "neutral_wood_curve": compute_wood_curve(0.0, 0.0, 0.0),
    }

    return {
        "summary": summary,
        "ufs": ufs,
    }


def render_html_dashboard(payload: dict[str, Any] | None = None) -> str:
    """
    Renders the standalone interactive HTML5 dashboard document with modern CSS,
    SVG vector geometries, chart containers, and audit table.
    """
    from rufas_lactation_equivalence.template import render_html_dashboard as _render_html

    return _render_html(payload=payload)

