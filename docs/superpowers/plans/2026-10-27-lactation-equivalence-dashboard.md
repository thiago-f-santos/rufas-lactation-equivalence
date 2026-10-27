# Lactation Equivalence Interactive Dashboard Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Implement a self-contained, interactive HTML5 + CSS3 + JS dashboard with a vector SVG map of Brazil, detailed hover/click state exploration, Wood's lactation curves ($Y(t)$ DIM 1..305), and edaphoclimatic decomposition across all 27 Brazilian UFs mapped to US homologues.

**Architecture:** A standalone, zero-external-server HTML page (`docs/index.html`) generated deterministically by `scripts/build_dashboard.py` from `output/equivalence_audit_report.json` and `data/brazil_ufs_edaphoclimatic.json`. Official IBGE 27 UF geometries are projected into a crisp $600 \times 600$ SVG map. The client-side logic binds hover/click events to update Chart.js (with native SVG fallback) lactation curves, biophysical radar charts, US homologue cards, and a synchronized 27-UF audit table.

**Tech Stack:** Python 3.10+ (for data builder & validator), HTML5, Modern CSS3 (CSS Grid, Flexbox, CSS Variables), Vanilla JavaScript (ES6+), Chart.js (via CDN with inline SVG fallback), SVG 1.1, Pytest.

**Spec:** [`docs/superpowers/specs/2026-10-27-lactation-equivalence-dashboard-design.md`](file:///home/thiago/Projects/rufas-lactation-equivalence/docs/superpowers/specs/2026-10-27-lactation-equivalence-dashboard-design.md)

## Global Constraints

- Standalone portability: `docs/index.html` MUST work directly when opened via `file://` double-click in any modern browser without CORS errors.
- Data integrity: All 27 Brazilian UFs (IBGE codes `11` to `53`) and their mappings to RuFaS macro-regions (*Delta*, *Southeast*, *Mountain*, *Southern Plains*, *Appalachian*) MUST match `output/equivalence_audit_report.json` exactly.
- Wood curve precision: Regional adjustments must use exact scaling factors ($l_{\text{eff}} = 19.90 + \Delta l$, $m_{\text{eff}} = 0.247 + \Delta m \times 10^{-2}$, $n_{\text{eff}} = 0.003376 + \Delta n \times 10^{-4}$).
- Color consistency: Macro-region colors must be uniform across map, charts, and badges:
  - Delta: `#10b981` (Green)
  - Southeast: `#3b82f6` (Blue)
  - Mountain: `#f59e0b` (Amber)
  - Southern Plains: `#ef4444` (Red)
  - Appalachian: `#8b5cf6` (Purple)

---

### Task 1: Data Preprocessing and IBGE SVG Map Projection Generator

**Files:**
- Create: `src/rufas_lactation_equivalence/dashboard_builder.py`
- Create: `scripts/build_dashboard.py`
- Test: `tests/test_dashboard_integrity.py`

**Interfaces:**
- Consumes:
  - `output/equivalence_audit_report.json`
  - `data/brazil_ufs_edaphoclimatic.json`
- Produces:
  - `generate_dashboard_payload() -> dict`: Returns structured dictionary containing all 27 UFs, projected SVG paths, raw and normalized values, homologues, and Wood curve parameters.

- [ ] **Step 1: Write the failing test for dashboard data extraction and SVG projection**

```python
# tests/test_dashboard_integrity.py
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
```

- [ ] **Step 2: Run test to verify it fails**

Run: `rtk test tests/test_dashboard_integrity.py`
Expected: FAIL with `ModuleNotFoundError: No module named 'rufas_lactation_equivalence.dashboard_builder'`

- [ ] **Step 3: Implement `dashboard_builder.py` and `scripts/build_dashboard.py`**

Create `src/rufas_lactation_equivalence/dashboard_builder.py` implementing:
- `project_coordinates_to_svg(lon, lat, width, height)` using Mercator projection with bounds `lon: [-74.0, -34.5], lat: [-33.8, 5.5]`.
- Loading of official IBGE geometries (cached or local polygon representation).
- Extraction of raw edaphoclimatic variables from `data/brazil_ufs_edaphoclimatic.json`.
- Extraction of audit equivalences, Gower similarities, and Wood adjustments from `output/equivalence_audit_report.json`.
- Calculation of Wood curve daily yield $Y(t) = l \cdot t^m \cdot e^{-nt}$ for $t \in [1, 305]$, $t_{\text{peak}} = m/n$, $Y_{\text{peak}} = l (m/n)^m e^{-m}$, and 305-day cumulative sum.
- Compilation of the payload into `scripts/build_dashboard.py`.

- [ ] **Step 4: Run test to verify it passes**

Run: `rtk test tests/test_dashboard_integrity.py`
Expected: PASS

- [ ] **Step 5: Commit**

```bash
rtk git -C /home/thiago/Projects/rufas-lactation-equivalence add src/rufas_lactation_equivalence/dashboard_builder.py scripts/build_dashboard.py tests/test_dashboard_integrity.py
rtk git -C /home/thiago/Projects/rufas-lactation-equivalence commit -m "feat: add dashboard data preprocessor and SVG projection generator"
```

---

### Task 2: HTML5 Layout and Responsive Modern CSS Architecture

**Files:**
- Create: `src/rufas_lactation_equivalence/template.py` (or HTML template module)
- Modify: `src/rufas_lactation_equivalence/dashboard_builder.py`
- Modify: `tests/test_dashboard_integrity.py`

**Interfaces:**
- Consumes:
  - `generate_dashboard_payload()`
- Produces:
  - HTML structure with CSS tokens, cards, SVG container, Chart containers, and audit table skeleton.

- [ ] **Step 1: Write the failing test for HTML and CSS generation**

```python
# append to tests/test_dashboard_integrity.py
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
```

- [ ] **Step 2: Run test to verify it fails**

Run: `rtk test tests/test_dashboard_integrity.py::test_render_html_dashboard`
Expected: FAIL with `cannot import name 'render_html_dashboard'`

- [ ] **Step 3: Implement HTML layout and CSS stylesheet generator**

Implement `render_html_dashboard()` in `src/rufas_lactation_equivalence/dashboard_builder.py`:
- Root CSS tokens for light/scientific clean theme and dark contrast cards.
- Layout: Header with key metrics (27 UFs, 5 RuFaS regions, 81.3% avg similarity).
- Main view: 2-column responsive layout (`@media (min-width: 1024px)`).
  - Left column: Brazil SVG map container, region filter pill badges, zoom/reset controls.
  - Right column: Active UF header, US Homologue card, Wood Lactation Curve chart, 6 Edaphoclimatic Factors Radar chart, Biophysical metrics cards.
- Bottom view: Full-width responsive audit table with instant search bar and sorting headers.
- Embedded `<script id="dashboard-data" type="application/json">` containing complete payload.

- [ ] **Step 4: Run test to verify it passes**

Run: `rtk test tests/test_dashboard_integrity.py::test_render_html_dashboard`
Expected: PASS

- [ ] **Step 5: Commit**

```bash
rtk git -C /home/thiago/Projects/rufas-lactation-equivalence add src/rufas_lactation_equivalence/dashboard_builder.py tests/test_dashboard_integrity.py
rtk git -C /home/thiago/Projects/rufas-lactation-equivalence commit -m "feat: implement HTML5 template and modern responsive CSS layout"
```

---

### Task 3: Interactive SVG Brazil Map with Hover and Click Interactivity

**Files:**
- Modify: `src/rufas_lactation_equivalence/dashboard_builder.py`
- Modify: `tests/test_dashboard_integrity.py`

**Interfaces:**
- Consumes:
  - 27 UF SVG `<path>` elements with attributes `id="uf-{sigla}"`, `data-uf="{sigla}"`, `data-ibge="{code}"`, `data-region="{region}"`.
- Produces:
  - Client-side JavaScript handling hover tooltip, path highlighting, region filtering, and click selection state.

- [ ] **Step 1: Write test for SVG map paths and client script elements**

```python
# append to tests/test_dashboard_integrity.py
def test_svg_map_elements_and_interactivity():
    from rufas_lactation_equivalence.dashboard_builder import render_html_dashboard

    html = render_html_dashboard()
    for uf in [
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
    ]:
        assert f'id="uf-{uf}"' in html
        assert f'data-uf="{uf}"' in html

    assert 'id="map-tooltip"' in html
    assert "selectState(" in html
    assert "showTooltip(" in html
```

- [ ] **Step 2: Run test to verify it fails**

Run: `rtk test tests/test_dashboard_integrity.py::test_svg_map_elements_and_interactivity`
Expected: FAIL

- [ ] **Step 3: Implement SVG path rendering and client-side map controller**

In `src/rufas_lactation_equivalence/dashboard_builder.py`:
- Render all 27 `<path>` elements with appropriate fill colors based on their RuFaS region.
- Client JavaScript:
  - Event listener on mousemove/mouseenter for `#brazil-map path`:
    - Shows `#map-tooltip` positioned at cursor with UF, US Homologue, RuFaS region, Gower similarity, and peak milk yield.
    - Adds `.hovered` class to the path.
  - Event listener on click:
    - Calls `selectState(ufSigla)`.
    - Updates `.active` class on path, removing previous active state.
    - Sincroniza o dropdown seletor de UFs.
  - Region filter pills: clicking a pill (e.g. "Delta") highlights all Delta states and fades others.

- [ ] **Step 4: Run test to verify it passes**

Run: `rtk test tests/test_dashboard_integrity.py::test_svg_map_elements_and_interactivity`
Expected: PASS

- [ ] **Step 5: Commit**

```bash
rtk git -C /home/thiago/Projects/rufas-lactation-equivalence add src/rufas_lactation_equivalence/dashboard_builder.py tests/test_dashboard_integrity.py
rtk git -C /home/thiago/Projects/rufas-lactation-equivalence commit -m "feat: add interactive SVG map controller with hover tooltip and click selection"
```

---

### Task 4: Wood Lactation Curve and Edaphoclimatic Radar Chart Integration

**Files:**
- Modify: `src/rufas_lactation_equivalence/dashboard_builder.py`
- Modify: `tests/test_dashboard_integrity.py`

**Interfaces:**
- Consumes:
  - `payload["ufs"][uf]["wood_curve"]`
  - `payload["ufs"][uf]["normalized_vector"]`
  - `payload["ufs"][uf]["best_match"]["normalized_vector"]`
- Produces:
  - Responsive Chart.js line chart for lactation curve ($Y(t)$ DIM 1..305) with neutral baseline comparison.
  - Radar chart for 6 edaphoclimatic factors comparing Brazil UF vs US Homologue.
  - Inline SVG fallback renderer if Chart.js CDN is unavailable.

- [ ] **Step 1: Write test for curve and radar chart client logic**

```python
# append to tests/test_dashboard_integrity.py
def test_chart_logic_and_wood_parameters():
    from rufas_lactation_equivalence.dashboard_builder import render_html_dashboard

    html = render_html_dashboard()
    assert "https://cdn.jsdelivr.net/npm/chart.js" in html
    assert "updateLactationChart" in html
    assert "updateRadarChart" in html
    assert "renderSvgFallbackLactation" in html
    assert "peak-day-val" in html
    assert "peak-yield-val" in html
    assert "total-305-val" in html
```

- [ ] **Step 2: Run test to verify it fails**

Run: `rtk test tests/test_dashboard_integrity.py::test_chart_logic_and_wood_parameters`
Expected: FAIL

- [ ] **Step 3: Implement Chart.js bindings, KPI metric updates, and SVG fallback**

In `src/rufas_lactation_equivalence/dashboard_builder.py`:
- Client JavaScript:
  - `updateLactationChart(ufData)`:
    - Plots DIM 1 to 305 with selected UF adjusted curve ($Y_t$).
    - Plots neutral unadjusted RuFaS baseline ($l=19.9, m=0.247, n=0.003376$).
    - Highlights peak day marker.
    - Updates DOM KPI elements: Peak Day (DIM), Peak Yield (kg/d), Cumulative 305d Yield (kg), $\Delta l$, $\Delta m$, $\Delta n$.
  - `updateRadarChart(ufData)`:
    - Plots 6 axes: THI, Sazonalidade Térmica, Índice Hídrico, Argila Solo, PAWC, Fertilidade Química.
    - Dataset 1: UF selecionada (preenchimento suave na cor da região).
    - Dataset 2: Estado americano homólogo (linha tracejada).
  - `renderSvgFallbackLactation(container, ufData)`:
    - If Chart.js fails to load (offline), dynamically draws `<svg>` line chart using native polyline.
  - Updates biophysical raw metrics grid (Temperature °C, RH %, Annual Precip mm, Clay %, PAWC, pH, SOC g/kg).

- [ ] **Step 4: Run test to verify it passes**

Run: `rtk test tests/test_dashboard_integrity.py::test_chart_logic_and_wood_parameters`
Expected: PASS

- [ ] **Step 5: Commit**

```bash
rtk git -C /home/thiago/Projects/rufas-lactation-equivalence add src/rufas_lactation_equivalence/dashboard_builder.py tests/test_dashboard_integrity.py
rtk git -C /home/thiago/Projects/rufas-lactation-equivalence commit -m "feat: integrate interactive Wood lactation curve and edaphoclimatic radar charts"
```

---

### Task 5: Searchable, Sortable Audit Table & State Synchronization

**Files:**
- Modify: `src/rufas_lactation_equivalence/dashboard_builder.py`
- Modify: `tests/test_dashboard_integrity.py`

**Interfaces:**
- Consumes:
  - Full list of 27 UFs from embedded data payload.
- Produces:
  - Interactive table with instant filtering by UF/state/region and multi-column sorting.
  - Bidirectional click sync: clicking a table row triggers `selectState(uf)`.

- [ ] **Step 1: Write test for table generation and synchronization logic**

```python
# append to tests/test_dashboard_integrity.py
def test_audit_table_structure_and_sync():
    from rufas_lactation_equivalence.dashboard_builder import render_html_dashboard

    html = render_html_dashboard()
    assert 'id="auditTableBody"' in html
    assert 'id="tableSearchInput"' in html
    assert "filterTable(" in html
    assert "sortTable(" in html
    assert "data-sort-key" in html
```

- [ ] **Step 2: Run test to verify it fails**

Run: `rtk test tests/test_dashboard_integrity.py::test_audit_table_structure_and_sync`
Expected: FAIL

- [ ] **Step 3: Implement interactive table renderer and event bindings**

In `src/rufas_lactation_equivalence/dashboard_builder.py`:
- Generate table rows with: IBGE Code, UF badge, State Name, Region BR, US Homologue, RuFaS Region badge, Gower Similarity %, $\Delta l$, $\Delta m$, $\Delta n$, Peak Yield, 305d Total Yield.
- Client JavaScript:
  - `filterTable(query)`: filters rows based on search input.
  - `sortTable(columnKey, isNumeric)`: sorts table rows ascending/descending.
  - Row click handler: calls `selectState(uf)` and smoothly scrolls to top view if desired.
  - "Copy as CSV" and "Export JSON" utility button handlers.

- [ ] **Step 4: Run test to verify it passes**

Run: `rtk test tests/test_dashboard_integrity.py::test_audit_table_structure_and_sync`
Expected: PASS

- [ ] **Step 5: Commit**

```bash
rtk git -C /home/thiago/Projects/rufas-lactation-equivalence add src/rufas_lactation_equivalence/dashboard_builder.py tests/test_dashboard_integrity.py
rtk git -C /home/thiago/Projects/rufas-lactation-equivalence commit -m "feat: add interactive searchable and sortable audit table with map sync"
```

---

### Task 6: CLI Integration, Final Build of `docs/index.html` & Full Verification

**Files:**
- Modify: `src/rufas_lactation_equivalence/cli.py`
- Modify: `scripts/build_dashboard.py`
- Output: `docs/index.html`
- Test: `tests/test_dashboard_integrity.py`

**Interfaces:**
- Consumes:
  - CLI flag `--build-dashboard [PATH]`
- Produces:
  - Generated `docs/index.html` file tested for size, completeness, and offline execution.

- [ ] **Step 1: Write test for CLI dashboard export and generated file integrity**

```python
# append to tests/test_dashboard_integrity.py
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
```

- [ ] **Step 2: Run test to verify it fails**

Run: `rtk test tests/test_dashboard_integrity.py::test_cli_dashboard_export`
Expected: FAIL

- [ ] **Step 3: Wire CLI and build script, generate `docs/index.html`**

- Update `src/rufas_lactation_equivalence/cli.py` to add `--build-dashboard` argument.
- Execute `python3 scripts/build_dashboard.py` to generate `docs/index.html`.
- Verify the file is generated, valid UTF-8, and opens without errors.

- [ ] **Step 4: Run entire test suite to verify 100% pass**

Run: `rtk test`
Expected: All tests pass (existing tests + new integrity tests).

- [ ] **Step 5: Commit and Finalize**

```bash
rtk git -C /home/thiago/Projects/rufas-lactation-equivalence add docs/index.html src/rufas_lactation_equivalence/cli.py scripts/build_dashboard.py tests/test_dashboard_integrity.py
rtk git -C /home/thiago/Projects/rufas-lactation-equivalence commit -m "feat: build and export standalone interactive lactation equivalence dashboard"
```
