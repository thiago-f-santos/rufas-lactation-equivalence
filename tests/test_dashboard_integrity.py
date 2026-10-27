"""Integrity tests for dashboard payload generation and map projection."""

from pathlib import Path
import pytest
from rufas_lactation_equivalence.dashboard_builder import (
    generate_dashboard_payload,
    project_coordinates_to_svg,
)


def test_project_coordinates_to_svg():
    # Lon: -74 to -34.8, Lat: -33.8 to 5.3 -> x: 0..600, y: 0..600
    x, y = project_coordinates_to_svg(-73.98, 5.27, width=600, height=600)
    assert 0 <= x <= 50
    assert 0 <= y <= 50

    x_se, y_se = project_coordinates_to_svg(-34.80, -33.74, width=600, height=600)
    assert 550 <= x_se <= 600
    assert 550 <= y_se <= 600


def test_generate_dashboard_payload():
    payload = generate_dashboard_payload()
    assert "ufs" in payload
    assert len(payload["ufs"]) == 27
    assert "summary" in payload

    sp = payload["ufs"]["SP"]
    assert sp["name"] == "São Paulo"
    assert sp["ibge_code"] == "35"
    assert sp["rufas_region"] == "Southeast"
    assert sp["best_match"]["us_abbrev"] == "GA"
    assert "svg_path" in sp
    assert len(sp["svg_path"]) > 20
    assert "wood_curve" in sp
    assert len(sp["wood_curve"]["dim"]) == 305
    assert len(sp["wood_curve"]["yield_kg"]) == 305
    assert sp["wood_curve"]["peak_day"] > 0
    assert sp["wood_curve"]["peak_yield"] > 0
