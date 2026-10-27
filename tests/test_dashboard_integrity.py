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


def test_render_html_dashboard():
    from rufas_lactation_equivalence.dashboard_builder import render_html_dashboard

    html = render_html_dashboard()
    assert "<!DOCTYPE html>" in html
    assert '<svg id="brazil-map"' in html
    assert '<canvas id="lactationChart"' in html
    assert '<canvas id="radarChart"' in html
    assert 'id="auditTable"' in html
    assert "--color-delta: #10b981" in html
    assert "--color-southeast: #3b82f6" in html
    assert "--color-mountain: #f59e0b" in html
    assert "--color-southern-plains: #ef4444" in html
    assert "--color-appalachian: #8b5cf6" in html


def test_svg_map_elements_and_interactivity():
    import re
    from rufas_lactation_equivalence.dashboard_builder import render_html_dashboard

    html = render_html_dashboard()
    ufs = [
        "AC",
        "AL",
        "AM",
        "AP",
        "BA",
        "CE",
        "DF",
        "ES",
        "GO",
        "MA",
        "MG",
        "MS",
        "MT",
        "PA",
        "PB",
        "PE",
        "PI",
        "PR",
        "RJ",
        "RN",
        "RO",
        "RR",
        "RS",
        "SC",
        "SE",
        "SP",
        "TO",
    ]
    for uf in ufs:
        assert f'id="uf-{uf}"' in html
        assert f'data-uf="{uf}"' in html
        match = re.search(rf'<path id="uf-{uf}"\s+data-uf="{uf}"\s+data-ibge="(\d+)"\s+data-region="([^"]+)"', html)
        assert match is not None, f"UF {uf} path does not have correct data-ibge and data-region attributes"

    # Tooltip and client-side map controller assertions
    assert 'id="map-tooltip"' in html
    assert "selectState(" in html
    assert "showTooltip(" in html
    assert "hideTooltip(" in html
    assert "filterByRegion(" in html

    # Tooltip content fields checking
    assert "Homólogo EUA:" in html
    assert "Similaridade Gower:" in html
    assert "Pico de Produção:" in html

    # Hovered state class handling
    assert ".uf-path.hovered" in html
    assert 'classList.add("hovered")' in html
    assert 'classList.remove("hovered")' in html

    # State select and region filter bindings
    assert 'id="stateSelect"' in html
    assert 'selectState(this.value)' in html
    assert 'classList.add("active")' in html
    assert 'classList.remove("active")' in html
    assert 'classList.add("faded")' in html
    assert 'classList.remove("faded")' in html


