"""
Scientific API data fetcher for NASA POWER Agroclimatology and ISRIC SoilGrids v2.0 WCS.

Pure stdlib Python 3.12+ implementation without external dependencies.
"""

from __future__ import annotations

import concurrent.futures
import contextlib
import json
import math
from pathlib import Path
import socket
import struct
import subprocess
import time
import urllib.error
import urllib.request
import zlib

from rufas_lactation_equivalence.atomic_io import atomic_write_json
from rufas_lactation_equivalence.models import RawEdaphoclimaticProfile


@contextlib.contextmanager
def forced_ipv4():
    """Context manager to force IPv4 resolution to prevent IPv6 routing stalls on Linux hosts."""
    orig_getaddrinfo = socket.getaddrinfo

    def _ipv4_getaddrinfo(*args, **kwargs):
        responses = orig_getaddrinfo(*args, **kwargs)
        ipv4_responses = [r for r in responses if r[0] == socket.AF_INET]
        return ipv4_responses if ipv4_responses else responses

    socket.getaddrinfo = _ipv4_getaddrinfo
    try:
        yield
    finally:
        socket.getaddrinfo = orig_getaddrinfo


def compute_thornthwaite_or_hargreaves_pet(monthly_t2m: list[float], lat: float) -> float:
    """
    Computes annual Potential Evapotranspiration (PET, mm/year) using Thornthwaite (1948).

    Parameters
    ----------
    monthly_t2m : list[float]
        12 monthly mean temperatures (°C), from Jan to Dec.
    lat : float
        Latitude in decimal degrees (-90.0 to 90.0).

    Returns
    -------
    float
        Annual PET clamped to [300.0, 2600.0] mm/year.
    """
    if len(monthly_t2m) != 12:
        raise ValueError(f"Expected 12 monthly temperatures, got {len(monthly_t2m)}")

    days_in_month = (31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31)
    mid_days = (15, 45, 74, 105, 135, 165, 196, 227, 258, 288, 319, 349)

    # Annual heat index I = sum((T / 5) ^ 1.514) for T > 0
    heat_index = sum((max(0.0, t) / 5.0) ** 1.514 for t in monthly_t2m)
    i_safe = max(0.01, heat_index)

    # Thornthwaite cubic polynomial exponent a
    a = (
        6.75e-7 * (i_safe**3)
        - 7.71e-5 * (i_safe**2)
        + 1.792e-2 * i_safe
        + 0.49239
    )

    lat_rad = math.radians(lat)
    annual_pet = 0.0

    for m, temp_c in enumerate(monthly_t2m):
        if temp_c <= 0.0:
            continue

        j = mid_days[m]
        sol_dec = 0.409 * math.sin((2.0 * math.pi * j / 365.0) - 1.39)
        tan_prod = -math.tan(lat_rad) * math.tan(sol_dec)
        tan_prod = max(-1.0, min(1.0, tan_prod))
        omega_s = math.acos(tan_prod)
        daylight_hours = (24.0 / math.pi) * omega_s
        daylight_factor = (daylight_hours / 12.0) * (days_in_month[m] / 30.0)

        if temp_c >= 26.5:
            # Empirical high-temperature extension
            unadj_pet = -415.85 + 32.24 * temp_c - 0.43 * (temp_c**2)
        else:
            unadj_pet = 16.0 * ((10.0 * temp_c / i_safe) ** a)

        monthly_pet = daylight_factor * max(0.0, unadj_pet)
        annual_pet += monthly_pet

    return max(300.0, min(2600.0, annual_pet))


def parse_nasa_power_response(payload: dict, lat: float) -> dict[str, float]:
    """
    Extracts edaphoclimatic parameters from NASA POWER climatology JSON response.
    """
    params = payload.get("properties", {}).get("parameter", {})
    t2m = params.get("T2M", {})
    rh2m = params.get("RH2M", {})
    precip = params.get("PRECTOTCORR", {})

    temp_mean_c = float(t2m.get("ANN", 0.0))
    rh_mean_pct = float(rh2m.get("ANN", 0.0))

    months = ["JAN", "FEB", "MAR", "APR", "MAY", "JUN", "JUL", "AUG", "SEP", "OCT", "NOV", "DEC"]
    monthly_temps = [float(t2m.get(m, temp_mean_c)) for m in months]

    temp_seasonality_c = max(monthly_temps) - min(monthly_temps)
    precip_annual_mm = float(precip.get("ANN", 0.0)) * 365.0
    pet_annual_mm = compute_thornthwaite_or_hargreaves_pet(monthly_temps, lat)

    return {
        "temp_mean_c": round(temp_mean_c, 3),
        "rh_mean_pct": round(rh_mean_pct, 3),
        "temp_seasonality_c": round(temp_seasonality_c, 3),
        "precip_annual_mm": round(precip_annual_mm, 2),
        "pet_annual_mm": round(pet_annual_mm, 2),
    }


def weighted_topsoil_average(val_0_5: float, val_5_15: float, val_15_30: float) -> float:
    """Computes depth-weighted average across 0-30 cm: (5*d1 + 10*d2 + 15*d3) / 30."""
    return (5.0 * val_0_5 + 10.0 * val_5_15 + 15.0 * val_15_30) / 30.0


def decode_soilgrids_tiff(tiff_bytes: bytes) -> float:
    """
    Decodes single-band int16 deflate-compressed GeoTIFF from ISRIC SoilGrids WCS 2.0.1.

    Extracts valid soil pixels (> 0) and computes their mean.
    Raises ValueError if no valid > 0 pixels are present.
    """
    if len(tiff_bytes) < 8:
        raise ValueError(f"TIFF header too short ({len(tiff_bytes)} bytes)")

    endian_mark = tiff_bytes[:2]
    if endian_mark == b"II":
        endian = "<"
    elif endian_mark == b"MM":
        endian = ">"
    else:
        raise ValueError(f"Unknown TIFF endianness: {endian_mark}")

    magic, ifd_offset = struct.unpack(f"{endian}HI", tiff_bytes[2:8])
    if magic != 42:
        raise ValueError(f"Invalid TIFF magic number: {magic}")

    num_entries = struct.unpack(f"{endian}H", tiff_bytes[ifd_offset : ifd_offset + 2])[0]
    tags: dict[int, tuple[int, int, int]] = {}

    for i in range(num_entries):
        entry_offset = ifd_offset + 2 + i * 12
        tag, tag_type, count, val_or_offset = struct.unpack(
            f"{endian}HHII", tiff_bytes[entry_offset : entry_offset + 12]
        )
        tags[tag] = (tag_type, count, val_or_offset)

    width = tags[256][2] if 256 in tags else 0
    height = tags[257][2] if 257 in tags else 0
    predictor = tags.get(317, (0, 0, 1))[2]

    # Tile vs Strip
    if 322 in tags and 323 in tags:
        tile_w = tags[322][2]
        tile_h = tags[323][2]
        tile_offset = tags[324][2]
        tile_bytes = tags[325][2]
    elif 273 in tags and 279 in tags:
        tile_w = width
        tile_h = height
        tile_offset = tags[273][2]
        tile_bytes = tags[279][2]
    else:
        raise ValueError("TIFF contains neither Tile nor Strip offsets")

    raw_data = tiff_bytes[tile_offset : tile_offset + tile_bytes]
    decompressed = zlib.decompress(raw_data)

    num_shorts = len(decompressed) // 2
    unpacked = list(struct.unpack(f"{endian}{num_shorts}h", decompressed))

    if predictor == 2:
        for r in range(tile_h):
            row_start = r * tile_w
            for c in range(1, tile_w):
                idx = row_start + c
                if idx < len(unpacked):
                    diff = unpacked[idx]
                    prev = unpacked[idx - 1]
                    res = (diff + prev) & 0xFFFF
                    if res >= 0x8000:
                        res -= 0x10000
                    unpacked[idx] = res

    valid_values: list[int] = []
    for r in range(height):
        for c in range(width):
            idx = r * tile_w + c
            if idx < len(unpacked):
                val = unpacked[idx]
                if val > 0:
                    valid_values.append(val)

    if not valid_values:
        raise ValueError("No valid > 0 soil pixels found in TIFF")

    return float(sum(valid_values) / len(valid_values))


def _http_get(url: str, timeout: int = 15, max_retries: int = 3) -> bytes:
    """Reliable HTTP GET with IPv4 retry and subprocess curl fallback."""
    for attempt in range(max_retries):
        try:
            req = urllib.request.Request(
                url,
                headers={"User-Agent": "RuFaS-Research/1.0 (Lactation-Equivalence; Python-3.12)"},
            )
            with forced_ipv4():
                with urllib.request.urlopen(req, timeout=timeout) as response:
                    return response.read()
        except Exception:
            if attempt == max_retries - 1:
                # Fallback to system curl with forced IPv4
                cmd = ["curl", "-4", "-sS", "--retry", "2", "--max-time", str(timeout), url]
                try:
                    res = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, check=True)
                    return res.stdout
                except (FileNotFoundError, subprocess.CalledProcessError) as exc:
                    raise RuntimeError(f"Failed to fetch {url} via curl fallback: {exc}") from exc
            time.sleep(1.0 * (attempt + 1))
    raise RuntimeError(f"Failed to fetch {url} after {max_retries} attempts")


def fetch_nasa_power_climate(lat: float, lon: float, max_retries: int = 3) -> dict[str, float]:
    """Fetches and parses 30-year climatological normals from NASA POWER."""
    url = (
        f"https://power.larc.nasa.gov/api/temporal/climatology/point"
        f"?parameters=T2M,RH2M,PRECTOTCORR&community=AG"
        f"&longitude={lon:.4f}&latitude={lat:.4f}&format=JSON"
    )
    raw = _http_get(url, timeout=20, max_retries=max_retries)
    payload = json.loads(raw.decode("utf-8"))
    return parse_nasa_power_response(payload, lat=lat)


def fetch_isric_soilgrids(lat: float, lon: float, max_retries: int = 3) -> dict[str, float]:
    """
    Fetches 0-30 cm topsoil properties from ISRIC SoilGrids v2.0 WCS 2.0.1.

    Returns clay (%), sand (%), ph (pH), and soc (g/kg).
    """
    props = ("clay", "sand", "phh2o", "soc")
    depths = ("0-5cm", "5-15cm", "15-30cm")
    pairs = [(p, d) for p in props for d in depths]

    curr_lat, curr_lon = lat, lon

    for attempt in range(max_retries):
        try:
            results: dict[str, dict[str, float]] = {p: {} for p in props}

            def _fetch_layer(pd: tuple[str, str]) -> tuple[str, str, float]:
                p, d = pd
                url = (
                    f"https://maps.isric.org/mapserv?map=/map/{p}.map"
                    f"&SERVICE=WCS&VERSION=2.0.1&REQUEST=GetCoverage"
                    f"&COVERAGEID={p}_{d}_mean&FORMAT=image/tiff"
                    f"&SUBSET=X({curr_lon-0.02:.2f},{curr_lon+0.02:.2f})"
                    f"&SUBSET=Y({curr_lat-0.02:.2f},{curr_lat+0.02:.2f})"
                    f"&SUBSETTINGCRS=http://www.opengis.net/def/crs/EPSG/0/4326"
                    f"&OUTPUTCRS=http://www.opengis.net/def/crs/EPSG/0/4326"
                )
                tiff_bytes = _http_get(url, timeout=20, max_retries=3)
                val = decode_soilgrids_tiff(tiff_bytes)
                return p, d, val

            with concurrent.futures.ThreadPoolExecutor(max_workers=6) as pool:
                layer_results = list(pool.map(_fetch_layer, pairs))

            for p, d, val in layer_results:
                results[p][d] = val

            raw_clay = weighted_topsoil_average(
                results["clay"]["0-5cm"], results["clay"]["5-15cm"], results["clay"]["15-30cm"]
            )
            raw_sand = weighted_topsoil_average(
                results["sand"]["0-5cm"], results["sand"]["5-15cm"], results["sand"]["15-30cm"]
            )
            raw_phh2o = weighted_topsoil_average(
                results["phh2o"]["0-5cm"], results["phh2o"]["5-15cm"], results["phh2o"]["15-30cm"]
            )
            raw_soc = weighted_topsoil_average(
                results["soc"]["0-5cm"], results["soc"]["5-15cm"], results["soc"]["15-30cm"]
            )

            return {
                "soil_clay_pct": round(raw_clay / 10.0, 2),
                "soil_sand_pct": round(raw_sand / 10.0, 2),
                "soil_ph": round(raw_phh2o / 10.0, 2),
                "soil_soc_g_kg": round(raw_soc / 10.0, 2),
            }
        except (ValueError, RuntimeError):
            # Shift coordinate away from coastal/water pixel if needed
            curr_lat += 0.08
            curr_lon += 0.08
            if attempt == max_retries - 1:
                raise

    raise RuntimeError(f"Failed to fetch ISRIC SoilGrids for ({lat}, {lon})")


def fetch_profile(
    code: str,
    abbrev: str,
    name: str,
    country: str,
    rufas_region: str | None,
    lat: float,
    lon: float,
) -> RawEdaphoclimaticProfile:
    """Fetches climate and soil data to construct a RawEdaphoclimaticProfile."""
    climate = fetch_nasa_power_climate(lat=lat, lon=lon)
    soil = fetch_isric_soilgrids(lat=lat, lon=lon)

    return RawEdaphoclimaticProfile(
        code=code,
        abbrev=abbrev,
        name=name,
        country=country,
        rufas_region=rufas_region,
        temp_mean_c=climate["temp_mean_c"],
        rh_mean_pct=climate["rh_mean_pct"],
        temp_seasonality_c=climate["temp_seasonality_c"],
        precip_annual_mm=climate["precip_annual_mm"],
        pet_annual_mm=climate["pet_annual_mm"],
        soil_clay_pct=soil["soil_clay_pct"],
        soil_sand_pct=soil["soil_sand_pct"],
        soil_ph=soil["soil_ph"],
        soil_soc_g_kg=soil["soil_soc_g_kg"],
    )


def fetch_all_profiles(
    coords_file: Path,
    output_us_file: Path,
    output_br_file: Path,
    max_workers: int = 8,
) -> None:
    """
    Fetches edaphoclimatic profiles for all US states and Brazilian UFs in coords_file,
    and atomically writes the output JSON files with API provenance metadata.
    """
    with open(coords_file, "r", encoding="utf-8") as f:
        coords_data = json.load(f)

    us_states = coords_data["us_states"]
    brazil_ufs = coords_data["brazil_ufs"]

    metadata = {
        "sources": {
            "climate": "NASA POWER Agroclimatology 30-yr Climatology (T2M, RH2M, PRECTOTCORR)",
            "soil": "ISRIC SoilGrids v2.0 WCS 2.0.1 0-30cm Depth-Weighted Average (clay, sand, phh2o, soc)",
            "pet_model": "Thornthwaite (1948) Annual PET",
        },
        "generated_at": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
    }

    def _fetch_us(item: dict) -> dict:
        p = fetch_profile(
            code=item["code"],
            abbrev=item["abbrev"],
            name=item["name"],
            country="USA",
            rufas_region=item["rufas_region"],
            lat=item["lat"],
            lon=item["lon"],
        )
        return {
            "code": p.code,
            "abbrev": p.abbrev,
            "name": p.name,
            "country": p.country,
            "rufas_region": p.rufas_region,
            "temp_mean_c": p.temp_mean_c,
            "rh_mean_pct": p.rh_mean_pct,
            "temp_seasonality_c": p.temp_seasonality_c,
            "precip_annual_mm": p.precip_annual_mm,
            "pet_annual_mm": p.pet_annual_mm,
            "soil_clay_pct": p.soil_clay_pct,
            "soil_sand_pct": p.soil_sand_pct,
            "soil_ph": p.soil_ph,
            "soil_soc_g_kg": p.soil_soc_g_kg,
        }

    def _fetch_br(item: dict) -> dict:
        p = fetch_profile(
            code=item["code"],
            abbrev=item["abbrev"],
            name=item["name"],
            country="BRA",
            rufas_region=None,
            lat=item["lat"],
            lon=item["lon"],
        )
        return {
            "code": p.code,
            "abbrev": p.abbrev,
            "name": p.name,
            "country": p.country,
            "rufas_region": p.rufas_region,
            "temp_mean_c": p.temp_mean_c,
            "rh_mean_pct": p.rh_mean_pct,
            "temp_seasonality_c": p.temp_seasonality_c,
            "precip_annual_mm": p.precip_annual_mm,
            "pet_annual_mm": p.pet_annual_mm,
            "soil_clay_pct": p.soil_clay_pct,
            "soil_sand_pct": p.soil_sand_pct,
            "soil_ph": p.soil_ph,
            "soil_soc_g_kg": p.soil_soc_g_kg,
        }

    with concurrent.futures.ThreadPoolExecutor(max_workers=max_workers) as executor:
        us_profiles = list(executor.map(_fetch_us, sorted(us_states, key=lambda s: int(s["code"]))))
        br_profiles = list(executor.map(_fetch_br, sorted(brazil_ufs, key=lambda u: int(u["code"]))))

    us_payload = {"metadata": metadata, "profiles": us_profiles}
    br_payload = {"metadata": metadata, "profiles": br_profiles}

    # Atomic write
    atomic_write_json(output_us_file, us_payload)
    atomic_write_json(output_br_file, br_payload)
