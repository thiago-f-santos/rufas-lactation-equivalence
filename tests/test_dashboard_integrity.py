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


def test_chart_logic_and_wood_parameters():
    from rufas_lactation_equivalence.dashboard_builder import render_html_dashboard

    html = render_html_dashboard()
    # Script and core function bindings
    assert "https://cdn.jsdelivr.net/npm/chart.js" in html
    assert "updateLactationChart" in html
    assert "updateRadarChart" in html
    assert "renderSvgFallbackLactation" in html

    # Lactation curve KPIs
    assert "peak-day-val" in html
    assert "peak-yield-val" in html
    assert "total-305-val" in html
    assert "wood-deltas-val" in html

    # Radar chart 6 edaphoclimatic axes
    assert "THI" in html
    assert "Sazonalidade Térmica" in html
    assert "Índice Hídrico" in html
    assert "Argila Solo" in html
    assert "PAWC" in html
    assert "Fertilidade Química" in html

    # Neutral baseline and SVG fallback
    assert "neutral_wood_curve" in html
    assert "polyline" in html

    # Biophysical metric grid
    assert "bio-pawc-val" in html
    assert "bio-temp-val" in html
    assert "bio-rh-val" in html
    assert "bio-clay-val" in html


def test_audit_table_structure_and_sync():
    from rufas_lactation_equivalence.dashboard_builder import (
        generate_dashboard_payload,
        render_html_dashboard,
    )

    # 1. Verify payload contains region_br
    payload = generate_dashboard_payload()
    assert payload["ufs"]["SP"]["region_br"] == "Sudeste"
    assert payload["ufs"]["RS"]["region_br"] == "Sul"
    assert payload["ufs"]["AM"]["region_br"] == "Norte"
    assert payload["ufs"]["BA"]["region_br"] == "Nordeste"
    assert payload["ufs"]["GO"]["region_br"] == "Centro-Oeste"

    # 2. Verify HTML and table elements
    html = render_html_dashboard()
    assert '<table id="auditTable"' in html
    assert 'id="auditTableBody"' in html
    assert 'id="tableSearchInput"' in html
    assert "filterTable(" in html
    assert "sortTable(" in html
    assert "data-sort-key" in html
    assert "copyTableAsCSV" in html
    assert "exportJsonData" in html

    # 3. Check all 27 Brazilian UFs have row with click handler
    ufs = [
        "AC", "AL", "AM", "AP", "BA", "CE", "DF", "ES", "GO", "MA",
        "MG", "MS", "MT", "PA", "PB", "PE", "PI", "PR", "RJ", "RN",
        "RO", "RR", "RS", "SC", "SE", "SP", "TO",
    ]
    for uf in ufs:
        assert f'id="row-{uf}"' in html
        assert f'data-uf="{uf}"' in html
        assert f'onclick="selectState(\'{uf}\')"' in html

    # 4. Check table column sort keys
    expected_keys = [
        "ibge_code",
        "uf",
        "name",
        "region_br",
        "rufas_region",
        "us_homologue",
        "gower_similarity",
        "delta_l",
        "delta_m",
        "delta_n",
        "peak_yield",
        "total_305d",
    ]
    for key in expected_keys:
        assert f'data-sort-key="{key}"' in html

    # 5. Check active-row synchronization logic
    assert "active-row" in html
    assert 'id="row-SP" data-uf="SP" class="active-row"' in html
    assert 'document.querySelectorAll("#auditTableBody tr").forEach' in html
    assert 'tableRow.classList.add("active-row")' in html


def test_cli_dashboard_export(tmp_path):
    from rufas_lactation_equivalence.dashboard_builder import build_and_save_dashboard

    dest = tmp_path / "index.html"
    build_and_save_dashboard(dest)
    assert dest.exists()
    content = dest.read_text(encoding="utf-8")
    assert len(content) > 50000  # Comprehensive HTML with SVG and JSON payload
    assert "São Paulo" in content
    assert "Rio Grande do Sul" in content
    assert "Amazonas" in content


def test_cli_build_dashboard_flag(tmp_path, monkeypatch):
    from rufas_lactation_equivalence.cli import main

    out_file = tmp_path / "custom_docs" / "index.html"
    test_args = ["rufas-lactation-equivalence", "--build-dashboard", str(out_file)]
    monkeypatch.setattr("sys.argv", test_args)
    main()

    assert out_file.exists()
    content = out_file.read_text(encoding="utf-8")
    assert len(content) > 50000
    assert "São Paulo" in content


def test_us_states_reference_and_match_audit():
    from rufas_lactation_equivalence.dashboard_builder import (
        generate_dashboard_payload,
        render_html_dashboard,
    )

    payload = generate_dashboard_payload()
    assert "us_states" in payload
    assert len(payload["us_states"]) == 48

    ga = payload["us_states"]["GA"]
    assert ga["state_name"] == "Georgia"
    assert ga["rufas_region"] == "Southeast"
    assert "SP" in ga["mapped_ufs"]
    assert "normalized_vector" in ga
    assert "raw_profile" in ga
    assert "wood_adjustments" in ga

    # Check match audit candidates for SP
    sp = payload["ufs"]["SP"]
    assert "match_audit" in sp
    assert len(sp["match_audit"]["top_candidates"]) >= 5
    top1 = sp["match_audit"]["top_candidates"][0]
    assert top1["us_abbrev"] == "GA"
    assert "factor_differences" in top1
    assert "gower_delta_from_top1" in top1

    # Check HTML elements for US reference table, tabs and candidate comparison
    html = render_html_dashboard()
    assert 'id="tab-brazil"' in html
    assert 'id="tab-us"' in html
    assert 'id="usAuditTable"' in html
    assert 'id="candidatesRanking"' in html
    assert "compareCandidateRadar(" in html
    assert "switchAuditTab(" in html





