import socket
import struct
import subprocess
from unittest.mock import patch
import zlib
import pytest

from rufas_lactation_equivalence.fetcher import (
    _http_get,
    compute_thornthwaite_or_hargreaves_pet,
    decode_soilgrids_tiff,
    forced_ipv4,
    parse_nasa_power_response,
    weighted_topsoil_average,
)


def test_parse_nasa_power_response_extracts_climate_metrics() -> None:
    """Verify extraction of T_mean, RH_mean, seasonality (max-min), precip_annual, and PET."""
    sample_payload = {
        "properties": {
            "parameter": {
                "T2M": {
                    "JAN": 23.0,
                    "FEB": 23.2,
                    "MAR": 22.8,
                    "APR": 21.5,
                    "MAY": 19.0,
                    "JUN": 17.5,
                    "JUL": 17.2,
                    "AUG": 19.5,
                    "SEP": 21.8,
                    "OCT": 22.9,
                    "NOV": 22.7,
                    "DEC": 22.8,
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


def test_compute_thornthwaite_or_hargreaves_pet_bounds_and_latitudes() -> None:
    """Verify PET calculation respects bounds [300.0, 2600.0] and behaves realistically."""
    # Cold temperate climate (e.g., Wisconsin / Minnesota)
    temperate_temps = [-8.0, -5.0, 1.0, 8.0, 15.0, 20.0, 23.0, 21.0, 16.0, 9.0, 2.0, -5.0]
    pet_temperate = compute_thornthwaite_or_hargreaves_pet(temperate_temps, lat=44.5)
    assert 400.0 <= pet_temperate <= 900.0

    # Tropical savanna climate (e.g., Minas Gerais / Goiás)
    tropical_temps = [24.0, 24.2, 24.0, 23.0, 21.5, 20.5, 20.5, 22.0, 24.5, 25.0, 24.5, 24.0]
    pet_tropical = compute_thornthwaite_or_hargreaves_pet(tropical_temps, lat=-18.5)
    assert 900.0 <= pet_tropical <= 1600.0

    # Extremely cold climate (clamped to lower bound 300.0)
    freezing_temps = [-20.0] * 12
    pet_cold = compute_thornthwaite_or_hargreaves_pet(freezing_temps, lat=65.0)
    assert pet_cold == 300.0

    # Extremely hot climate (clamped to upper bound 2600.0)
    hot_temps = [45.0] * 12
    pet_hot = compute_thornthwaite_or_hargreaves_pet(hot_temps, lat=25.0)
    assert pet_hot <= 2600.0


def _build_synthetic_geotiff(
    width: int,
    height: int,
    values: list[int],
    predictor: int = 2,
    tile_w: int = 16,
    tile_h: int = 16,
) -> bytes:
    """Builds a minimal valid little-endian GeoTIFF byte string for testing decode_soilgrids_tiff."""
    tile_pixels = [0] * (tile_w * tile_h)
    for r in range(height):
        for c in range(width):
            tile_pixels[r * tile_w + c] = values[r * width + c]

    raw_tile = list(tile_pixels)
    if predictor == 2:
        for r in range(tile_h):
            row_start = r * tile_w
            orig_row = list(raw_tile[row_start : row_start + tile_w])
            for c in range(1, tile_w):
                diff = (orig_row[c] - orig_row[c - 1]) & 0xFFFF
                if diff >= 0x8000:
                    diff -= 0x10000
                raw_tile[row_start + c] = diff

    packed_data = struct.pack(f"<{len(raw_tile)}h", *raw_tile)
    compressed_tile = zlib.compress(packed_data)

    ifd_offset = 8
    tags = [
        (256, 3, 1, width),  # ImageWidth
        (257, 3, 1, height),  # ImageLength
        (258, 3, 1, 16),  # BitsPerSample
        (259, 3, 1, 8),  # Compression (deflate)
        (317, 3, 1, predictor),  # Predictor
        (322, 3, 1, tile_w),  # TileWidth
        (323, 3, 1, tile_h),  # TileLength
        (324, 4, 1, 0),  # TileOffsets (placeholder)
        (325, 4, 1, len(compressed_tile)),  # TileByteCounts
    ]
    num_tags = len(tags)
    header_ifd_size = 8 + 2 + num_tags * 12 + 4
    tile_offset = header_ifd_size

    tags[7] = (324, 4, 1, tile_offset)

    buf = bytearray()
    buf.extend(struct.pack("<2sHI", b"II", 42, ifd_offset))
    buf.extend(struct.pack("<H", num_tags))
    for tag, tag_type, count, val_or_offset in tags:
        buf.extend(struct.pack("<HHII", tag, tag_type, count, val_or_offset))
    buf.extend(struct.pack("<I", 0))  # Next IFD offset
    buf.extend(compressed_tile)
    return bytes(buf)


def test_decode_soilgrids_tiff_synthetic() -> None:
    """Verify decode_soilgrids_tiff successfully decompresses, reverses differencing, and averages > 0."""
    w, h = 4, 4
    # Values with some nodata (-32768, which as unsigned 16-bit is 0x8000)
    values = [
        350, 360, 370, 380,
        340, -1,  360, 390,
        330, 340, 350, 360,
        320, 330, 340, 350,
    ]
    # Filter valid > 0
    valid_expected = [v for v in values if v > 0]
    expected_mean = sum(valid_expected) / len(valid_expected)

    tiff_bytes = _build_synthetic_geotiff(width=w, height=h, values=values, predictor=2, tile_w=16, tile_h=16)
    decoded_mean = decode_soilgrids_tiff(tiff_bytes)
    assert pytest.approx(decoded_mean) == expected_mean


def test_decode_soilgrids_tiff_all_nodata_raises() -> None:
    """Verify ValueError is raised if a tile contains no valid > 0 pixels."""
    values = [-1] * 16
    tiff_bytes = _build_synthetic_geotiff(width=4, height=4, values=values, predictor=2, tile_w=16, tile_h=16)
    with pytest.raises(ValueError, match="No valid > 0 soil pixels found"):
        decode_soilgrids_tiff(tiff_bytes)


def test_decode_soilgrids_tiff_big_endian() -> None:
    """Verify decode_soilgrids_tiff correctly decodes big-endian (MM) TIFFs."""
    w, h = 2, 2
    tile_w, tile_h = 2, 2
    values = [100, 200, 300, 400]
    packed = struct.pack(">4h", *values)
    compressed = zlib.compress(packed)
    ifd_offset = 8
    tags = [
        (256, 3, 1, w),
        (257, 3, 1, h),
        (258, 3, 1, 16),
        (259, 3, 1, 8),
        (317, 3, 1, 1),  # Predictor 1 (no differencing)
        (322, 3, 1, tile_w),
        (323, 3, 1, tile_h),
        (324, 4, 1, 0),
        (325, 4, 1, len(compressed)),
    ]
    header_ifd_size = 8 + 2 + len(tags) * 12 + 4
    tags[7] = (324, 4, 1, header_ifd_size)

    buf = bytearray()
    buf.extend(struct.pack(">2sHI", b"MM", 42, ifd_offset))
    buf.extend(struct.pack(">H", len(tags)))
    for tag, tag_type, count, val_or_offset in tags:
        buf.extend(struct.pack(">HHII", tag, tag_type, count, val_or_offset))
    buf.extend(struct.pack(">I", 0))
    buf.extend(compressed)

    decoded = decode_soilgrids_tiff(bytes(buf))
    assert pytest.approx(decoded) == 250.0


def test_decode_soilgrids_tiff_header_errors() -> None:
    """Verify TIFF header validation checks."""
    with pytest.raises(ValueError, match="TIFF header too short"):
        decode_soilgrids_tiff(b"short")

    with pytest.raises(ValueError, match="Unknown TIFF endianness"):
        decode_soilgrids_tiff(b"XX\x00\x2a\x00\x00\x00\x08")

    with pytest.raises(ValueError, match="Invalid TIFF magic number"):
        decode_soilgrids_tiff(b"II\x00\x99\x00\x00\x00\x08")


def test_forced_ipv4_context_manager() -> None:
    """Verify forced_ipv4 filters for AF_INET and restores original getaddrinfo afterwards."""
    orig_gai = socket.getaddrinfo

    dummy_responses = [
        (socket.AF_INET6, socket.SOCK_STREAM, 6, "", ("::1", 80, 0, 0)),
        (socket.AF_INET, socket.SOCK_STREAM, 6, "", ("127.0.0.1", 80)),
    ]

    with patch("socket.getaddrinfo", return_value=dummy_responses):
        before = socket.getaddrinfo
        with forced_ipv4():
            assert socket.getaddrinfo != before
            res = socket.getaddrinfo("localhost", 80)
            assert len(res) == 1
            assert res[0][0] == socket.AF_INET
        assert socket.getaddrinfo == before

    assert socket.getaddrinfo == orig_gai


def test_http_get_curl_fallback_error() -> None:
    """Verify _http_get raises RuntimeError cleanly when curl fallback fails."""
    with patch("urllib.request.urlopen", side_effect=Exception("network error")):
        with patch("subprocess.run", side_effect=FileNotFoundError("curl not found")):
            with pytest.raises(RuntimeError, match="Failed to fetch .* via curl fallback"):
                _http_get("http://example.com/test", timeout=1, max_retries=1)

        with patch(
            "subprocess.run",
            side_effect=subprocess.CalledProcessError(1, ["curl"], stderr=b"connection timed out"),
        ):
            with pytest.raises(RuntimeError, match="Failed to fetch .* via curl fallback"):
                _http_get("http://example.com/test", timeout=1, max_retries=1)

