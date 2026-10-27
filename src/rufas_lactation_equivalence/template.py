"""
HTML5 Template and Modern Responsive CSS Architecture for RuFaS Lactation Equivalence.

Generates a standalone, zero-dependency, self-contained dashboard visualization
with vector SVG map of Brazil, Wood lactation curve charts, edaphoclimatic radar
breakdowns, and an interactive 27-UF audit table.
"""

from __future__ import annotations

import json
from typing import Any

CSS_STYLES = """
:root {
  --color-delta: #10b981;
  --color-southeast: #3b82f6;
  --color-mountain: #f59e0b;
  --color-southern-plains: #ef4444;
  --color-appalachian: #8b5cf6;

  --color-bg: #f8fafc;
  --color-surface: #ffffff;
  --color-surface-elevated: #ffffff;
  --color-surface-subtle: #f1f5f9;
  --color-surface-muted: #e2e8f0;
  --color-border: #e2e8f0;
  --color-border-subtle: #cbd5e1;
  --color-border-focus: #3b82f6;

  --color-text-main: #0f172a;
  --color-text-secondary: #334155;
  --color-text-muted: #64748b;
  --color-text-light: #94a3b8;
  --color-text-inverse: #ffffff;

  --color-primary: #1e293b;
  --color-primary-accent: #2563eb;
  --color-neutral-curve: #94a3b8;

  --shadow-xs: 0 1px 2px 0 rgba(0, 0, 0, 0.05);
  --shadow-sm: 0 1px 3px 0 rgba(0, 0, 0, 0.08), 0 1px 2px -1px rgba(0, 0, 0, 0.08);
  --shadow-md: 0 4px 6px -1px rgba(0, 0, 0, 0.08), 0 2px 4px -2px rgba(0, 0, 0, 0.08);
  --shadow-lg: 0 10px 15px -3px rgba(0, 0, 0, 0.08), 0 4px 6px -4px rgba(0, 0, 0, 0.04);
  --shadow-xl: 0 20px 25px -5px rgba(0, 0, 0, 0.1), 0 8px 10px -6px rgba(0, 0, 0, 0.1);

  --radius-xs: 4px;
  --radius-sm: 6px;
  --radius-md: 10px;
  --radius-lg: 14px;
  --radius-xl: 20px;
  --radius-full: 9999px;

  --font-sans: -apple-system, BlinkMacSystemFont, "Segoe UI", Roboto, "Helvetica Neue", Arial, sans-serif;
  --font-mono: ui-monospace, SFMono-Regular, Menlo, Monaco, Consolas, "Liberation Mono", "Courier New", monospace;
}

*, *::before, *::after {
  box-sizing: border-box;
  margin: 0;
  padding: 0;
}

body {
  font-family: var(--font-sans);
  background-color: var(--color-bg);
  color: var(--color-text-main);
  line-height: 1.5;
  -webkit-font-smoothing: antialiased;
  -moz-osx-font-smoothing: grayscale;
}

.dashboard-container {
  max-width: 1480px;
  margin: 0 auto;
  padding: 24px 20px 60px;
}

/* Header */
.dashboard-header {
  margin-bottom: 24px;
}

.header-top {
  display: flex;
  align-items: center;
  justify-content: space-between;
  flex-wrap: wrap;
  gap: 12px;
  margin-bottom: 8px;
}

.brand-badge {
  display: inline-flex;
  align-items: center;
  gap: 6px;
  background: #e0f2fe;
  color: #0369a1;
  font-size: 0.75rem;
  font-weight: 700;
  text-transform: uppercase;
  letter-spacing: 0.06em;
  padding: 4px 10px;
  border-radius: var(--radius-full);
}

.header-version {
  font-size: 0.8rem;
  color: var(--color-text-muted);
}

.dashboard-title {
  font-size: 1.85rem;
  font-weight: 800;
  color: var(--color-text-main);
  letter-spacing: -0.025em;
  line-height: 1.25;
}

.dashboard-subtitle {
  font-size: 0.95rem;
  color: var(--color-text-muted);
  max-width: 960px;
  margin-top: 6px;
  line-height: 1.5;
}

/* Summary indicators KPI strip */
.summary-kpi-grid {
  display: grid;
  grid-template-columns: repeat(auto-fit, minmax(210px, 1fr));
  gap: 16px;
  margin-top: 20px;
}

.kpi-card {
  background: var(--color-surface);
  border: 1px solid var(--color-border);
  border-radius: var(--radius-md);
  padding: 14px 18px;
  box-shadow: var(--shadow-xs);
  transition: transform 0.15s ease, box-shadow 0.15s ease;
}

.kpi-card:hover {
  transform: translateY(-2px);
  box-shadow: var(--shadow-sm);
}

.kpi-label {
  font-size: 0.72rem;
  font-weight: 700;
  text-transform: uppercase;
  letter-spacing: 0.06em;
  color: var(--color-text-muted);
}

.kpi-value {
  font-size: 1.65rem;
  font-weight: 800;
  color: var(--color-text-main);
  margin-top: 4px;
  letter-spacing: -0.02em;
}

.kpi-meta {
  font-size: 0.75rem;
  color: var(--color-text-muted);
  margin-top: 2px;
}

/* Main 2-column layout */
.dashboard-main {
  display: grid;
  grid-template-columns: 1fr;
  gap: 24px;
  margin-bottom: 28px;
}

@media (min-width: 1024px) {
  .dashboard-main {
    grid-template-columns: 480px 1fr;
    align-items: start;
  }
}

@media (min-width: 1280px) {
  .dashboard-main {
    grid-template-columns: 520px 1fr;
  }
}

/* Cards */
.dash-card {
  background: var(--color-surface);
  border: 1px solid var(--color-border);
  border-radius: var(--radius-lg);
  box-shadow: var(--shadow-xs);
  overflow: hidden;
  margin-bottom: 24px;
}

.dash-card:last-child {
  margin-bottom: 0;
}

.card-header {
  padding: 16px 20px;
  border-bottom: 1px solid var(--color-border);
  display: flex;
  align-items: center;
  justify-content: space-between;
  flex-wrap: wrap;
  gap: 8px;
}

.card-title {
  font-size: 1.05rem;
  font-weight: 700;
  color: var(--color-text-main);
  display: flex;
  align-items: center;
  gap: 8px;
}

.card-subtitle {
  font-size: 0.8rem;
  color: var(--color-text-muted);
  margin-top: 2px;
}

.card-body {
  padding: 20px;
}

/* Left Column: Map & Region Filters */
.region-filters {
  display: flex;
  flex-wrap: wrap;
  gap: 6px;
  margin-bottom: 16px;
}

.region-pill {
  appearance: none;
  border: 1px solid var(--color-border);
  background: var(--color-surface);
  color: var(--color-text-secondary);
  border-radius: var(--radius-full);
  padding: 5px 12px;
  font-size: 0.78rem;
  font-weight: 600;
  cursor: pointer;
  transition: all 0.15s ease;
  display: inline-flex;
  align-items: center;
  gap: 6px;
}

.region-pill:hover {
  background: var(--color-surface-subtle);
  border-color: var(--color-border-subtle);
}

.region-pill.active {
  background: var(--color-primary);
  color: #ffffff;
  border-color: var(--color-primary);
  box-shadow: var(--shadow-xs);
}

.region-color-dot {
  width: 9px;
  height: 9px;
  border-radius: 50%;
  display: inline-block;
}

.map-viewport {
  position: relative;
  width: 100%;
  aspect-ratio: 1 / 1;
  max-width: 500px;
  margin: 0 auto;
  background: #fdfefe;
  border-radius: var(--radius-md);
  border: 1px solid var(--color-surface-subtle);
  overflow: hidden;
  display: flex;
  align-items: center;
  justify-content: center;
}

#brazil-map {
  width: 100%;
  height: 100%;
  display: block;
}

.uf-path {
  cursor: pointer;
  stroke: #ffffff;
  stroke-width: 1.2px;
  stroke-linejoin: round;
  transition: fill 0.2s ease, stroke 0.2s ease, stroke-width 0.2s ease, opacity 0.2s ease, filter 0.2s ease;
}

.uf-path:hover,
.uf-path.hovered {
  stroke: #0f172a;
  stroke-width: 2.2px;
  filter: brightness(1.15);
}

.uf-path.active {
  stroke: #0f172a;
  stroke-width: 2.8px;
  filter: brightness(1.2) drop-shadow(0 2px 6px rgba(0,0,0,0.3));
}

.uf-path.faded {
  opacity: 0.22;
  fill: #cbd5e1 !important;
}

/* Floating Tooltip */
.map-tooltip {
  position: absolute;
  pointer-events: none;
  z-index: 100;
  background: rgba(15, 23, 42, 0.94);
  backdrop-filter: blur(8px);
  -webkit-backdrop-filter: blur(8px);
  color: #ffffff;
  padding: 10px 14px;
  border-radius: var(--radius-sm);
  font-size: 0.8rem;
  line-height: 1.45;
  box-shadow: var(--shadow-lg);
  border: 1px solid rgba(255, 255, 255, 0.15);
  transition: opacity 0.15s ease;
  max-width: 260px;
}

.map-tooltip .tooltip-header {
  font-weight: 700;
  font-size: 0.9rem;
  margin-bottom: 4px;
  display: flex;
  align-items: center;
  justify-content: space-between;
}

.map-tooltip .tooltip-row {
  display: flex;
  justify-content: space-between;
  gap: 12px;
  margin-top: 2px;
  color: #cbd5e1;
}

.map-tooltip .tooltip-row strong {
  color: #ffffff;
}

.map-controls-bar {
  display: flex;
  align-items: center;
  justify-content: space-between;
  margin-top: 12px;
  padding-top: 10px;
  border-top: 1px solid var(--color-surface-subtle);
  font-size: 0.78rem;
  color: var(--color-text-muted);
}

.btn-map-control {
  background: var(--color-surface-subtle);
  border: 1px solid var(--color-border);
  padding: 4px 10px;
  border-radius: var(--radius-sm);
  font-size: 0.75rem;
  cursor: pointer;
  transition: background 0.15s;
}

.btn-map-control:hover {
  background: var(--color-surface-muted);
}

/* Right Column: Active State, Charts & Biophysical Breakdown */
.active-state-bar {
  display: flex;
  align-items: center;
  justify-content: space-between;
  flex-wrap: wrap;
  gap: 12px;
  background: var(--color-surface-subtle);
  border-radius: var(--radius-md);
  padding: 12px 16px;
  margin-bottom: 20px;
}

.state-select-wrapper {
  display: flex;
  align-items: center;
  gap: 8px;
}

.state-select-label {
  font-size: 0.8rem;
  font-weight: 600;
  color: var(--color-text-muted);
}

.state-select {
  padding: 6px 12px;
  border: 1px solid var(--color-border);
  border-radius: var(--radius-sm);
  background: var(--color-surface);
  font-family: var(--font-sans);
  font-size: 0.85rem;
  font-weight: 600;
  color: var(--color-text-main);
  cursor: pointer;
  outline: none;
  transition: border-color 0.15s;
}

.state-select:focus {
  border-color: var(--color-border-focus);
}

.active-state-identity {
  display: flex;
  align-items: center;
  gap: 10px;
}

.state-pill-badge {
  background: var(--color-primary);
  color: #ffffff;
  padding: 3px 8px;
  border-radius: var(--radius-xs);
  font-weight: 700;
  font-size: 0.85rem;
}

.state-name-display {
  font-size: 1.15rem;
  font-weight: 700;
  color: var(--color-text-main);
}

.region-tag-badge {
  padding: 3px 10px;
  border-radius: var(--radius-full);
  font-size: 0.75rem;
  font-weight: 700;
  color: #ffffff;
}

/* US Homologue Card */
.homologue-banner {
  background: linear-gradient(135deg, #f8fafc 0%, #eff6ff 100%);
  border: 1px solid #dbeafe;
  border-radius: var(--radius-md);
  padding: 14px 18px;
  margin-bottom: 20px;
  display: flex;
  align-items: center;
  justify-content: space-between;
  flex-wrap: wrap;
  gap: 16px;
}

.homologue-info {
  display: flex;
  flex-direction: column;
}

.homologue-meta {
  font-size: 0.72rem;
  text-transform: uppercase;
  letter-spacing: 0.05em;
  font-weight: 700;
  color: #2563eb;
}

.homologue-title {
  font-size: 1.15rem;
  font-weight: 800;
  color: #1e3a8a;
  margin-top: 2px;
}

.homologue-scores {
  display: flex;
  align-items: center;
  gap: 20px;
}

.score-badge {
  display: flex;
  flex-direction: column;
  align-items: flex-end;
}

.score-badge .score-label {
  font-size: 0.7rem;
  color: #64748b;
  text-transform: uppercase;
  font-weight: 600;
}

.score-badge .score-val {
  font-size: 1.35rem;
  font-weight: 800;
  color: #0369a1;
}

/* Lactation Metrics Quick Strip */
.lactation-metrics-strip {
  display: grid;
  grid-template-columns: repeat(auto-fit, minmax(140px, 1fr));
  gap: 12px;
  margin-bottom: 16px;
}

.lactation-metric-box {
  background: var(--color-surface-subtle);
  border-radius: var(--radius-sm);
  padding: 10px 14px;
}

.lactation-metric-box .metric-label {
  font-size: 0.7rem;
  color: var(--color-text-muted);
  text-transform: uppercase;
  font-weight: 600;
}

.lactation-metric-box .metric-value-wrap {
  display: flex;
  align-items: baseline;
  gap: 4px;
  margin-top: 3px;
}

.lactation-metric-box .metric-num {
  font-size: 1.35rem;
  font-weight: 800;
  color: var(--color-text-main);
}

.lactation-metric-box .metric-unit {
  font-size: 0.75rem;
  color: var(--color-text-muted);
  font-weight: 500;
}

/* Chart Canvas Wrapper */
.chart-container-box {
  position: relative;
  width: 100%;
  height: 270px;
  margin-top: 8px;
}

.chart-container-box canvas {
  width: 100% !important;
  height: 100% !important;
}

.svg-chart-fallback {
  width: 100%;
  height: 100%;
  display: flex;
  align-items: center;
  justify-content: center;
  background: var(--color-surface-subtle);
  border-radius: var(--radius-sm);
}

/* Biophysical Raw Metrics Grid */
.biophysical-grid {
  display: grid;
  grid-template-columns: repeat(auto-fit, minmax(130px, 1fr));
  gap: 10px;
  margin-top: 18px;
}

.bio-metric-cell {
  background: var(--color-surface-subtle);
  border: 1px solid var(--color-border);
  border-radius: var(--radius-sm);
  padding: 8px 12px;
}

.bio-metric-label {
  font-size: 0.68rem;
  text-transform: uppercase;
  color: var(--color-text-muted);
  font-weight: 600;
  letter-spacing: 0.04em;
}

.bio-metric-val {
  font-size: 1.05rem;
  font-weight: 700;
  color: var(--color-text-main);
  margin-top: 2px;
}

/* Bottom Audit Table */
.table-card {
  background: var(--color-surface);
  border: 1px solid var(--color-border);
  border-radius: var(--radius-lg);
  box-shadow: var(--shadow-xs);
  overflow: hidden;
}

.table-toolbar {
  padding: 16px 20px;
  display: flex;
  align-items: center;
  justify-content: space-between;
  flex-wrap: wrap;
  gap: 16px;
  border-bottom: 1px solid var(--color-border);
}

.table-search-box {
  position: relative;
  flex: 1;
  max-width: 440px;
}

.table-search {
  width: 100%;
  padding: 8px 14px 8px 36px;
  font-family: var(--font-sans);
  font-size: 0.85rem;
  border: 1px solid var(--color-border);
  border-radius: var(--radius-md);
  background: var(--color-surface-subtle);
  color: var(--color-text-main);
  outline: none;
  transition: all 0.15s ease;
}

.table-search:focus {
  background: var(--color-surface);
  border-color: var(--color-border-focus);
  box-shadow: 0 0 0 3px rgba(59, 130, 246, 0.15);
}

.search-icon {
  position: absolute;
  left: 12px;
  top: 50%;
  transform: translateY(-50%);
  color: var(--color-text-light);
  font-size: 0.9rem;
  pointer-events: none;
}

.table-actions {
  display: flex;
  align-items: center;
  gap: 10px;
}

.btn-utility {
  background: var(--color-surface);
  border: 1px solid var(--color-border);
  color: var(--color-text-secondary);
  padding: 7px 14px;
  border-radius: var(--radius-sm);
  font-size: 0.8rem;
  font-weight: 600;
  cursor: pointer;
  transition: all 0.15s ease;
  display: inline-flex;
  align-items: center;
  gap: 6px;
}

.btn-utility:hover {
  background: var(--color-surface-subtle);
  border-color: var(--color-border-subtle);
}

.table-scroll-wrapper {
  overflow-x: auto;
  max-height: 540px;
}

.audit-table {
  width: 100%;
  border-collapse: collapse;
  font-size: 0.82rem;
  text-align: left;
}

.audit-table thead th {
  position: sticky;
  top: 0;
  background: #f8fafc;
  color: var(--color-text-muted);
  font-weight: 700;
  text-transform: uppercase;
  font-size: 0.7rem;
  letter-spacing: 0.05em;
  padding: 10px 14px;
  border-bottom: 1px solid var(--color-border);
  cursor: pointer;
  user-select: none;
  white-space: nowrap;
  transition: background 0.15s;
}

.audit-table thead th:hover {
  background: #f1f5f9;
  color: var(--color-text-main);
}

.audit-table tbody tr {
  border-bottom: 1px solid var(--color-surface-subtle);
  cursor: pointer;
  transition: background 0.15s ease;
}

.audit-table tbody tr:hover {
  background: #f8fafc;
}

.audit-table tbody tr.active-row {
  background: #eff6ff !important;
  font-weight: 600;
}

.audit-table td {
  padding: 10px 14px;
  white-space: nowrap;
  color: var(--color-text-main);
}

.table-badge {
  display: inline-block;
  padding: 2px 7px;
  border-radius: var(--radius-xs);
  font-weight: 700;
  font-size: 0.75rem;
}

.table-badge-outline {
  display: inline-block;
  padding: 1px 6px;
  border-radius: var(--radius-xs);
  border: 1px solid var(--color-border-subtle);
  background: var(--color-surface-subtle);
  color: var(--color-text-secondary);
  font-weight: 600;
  font-size: 0.72rem;
}

.audit-table thead th.sort-asc::after {
  content: " ▲";
  color: var(--color-primary-accent);
}

.audit-table thead th.sort-desc::after {
  content: " ▼";
  color: var(--color-primary-accent);
}

.table-num {
  font-family: var(--font-mono);
  text-align: right;
}

.audit-table th.table-num {
  text-align: right;
}

/* Footer */
.dashboard-footer {
  margin-top: 40px;
  padding-top: 20px;
  border-top: 1px solid var(--color-border);
  display: flex;
  align-items: center;
  justify-content: space-between;
  flex-wrap: wrap;
  gap: 12px;
  font-size: 0.8rem;
  color: var(--color-text-muted);
}
""".strip()


def render_html_dashboard(payload: dict[str, Any] | None = None) -> str:
    """
    Renders the standalone interactive HTML5 dashboard document.

    Embeds:
    - Inlined responsive modern CSS3 stylesheet with exact scientific tokens
    - IBGE SVG vector paths for 27 UFs
    - Chart.js containers and native fallback handlers
    - 27-UF sortable and searchable audit table
    - Complete deterministic JSON payload in `<script id="dashboard-data">`
    """
    if payload is None:
        from rufas_lactation_equivalence.dashboard_builder import generate_dashboard_payload

        payload = generate_dashboard_payload()

    ufs = payload.get("ufs", {})
    summary = payload.get("summary", {})
    total_ufs = summary.get("total_ufs", len(ufs))
    total_regions = summary.get("total_regions", 5)
    avg_similarity = summary.get("avg_similarity", 0.813)
    min_similarity = summary.get("min_similarity", 0.75)
    max_similarity = summary.get("max_similarity", 0.94)

    # Initial selected state (default to SP or first UF)
    default_uf = "SP" if "SP" in ufs else next(iter(ufs.keys()))
    default_state = ufs.get(default_uf, {})
    default_match = default_state.get("best_match", {})
    default_curve = default_state.get("wood_curve", {})
    default_wood_adj = default_state.get("wood_adjustments", {})
    default_raw = default_state.get("raw_profile", {})

    # Generate 27 SVG <path> elements
    svg_paths_html: list[str] = []
    for sigla, udata in sorted(ufs.items(), key=lambda x: x[0]):
        code = udata.get("ibge_code", "")
        region = udata.get("rufas_region", "Southeast")
        color = udata.get("region_color", "#3b82f6")
        path_d = udata.get("svg_path", "")
        active_cls = " active" if sigla == default_uf else ""
        svg_paths_html.append(
            f'<path id="uf-{sigla}" data-uf="{sigla}" data-ibge="{code}" '
            f'data-region="{region}" fill="{color}" d="{path_d}" '
            f'class="uf-path{active_cls}" role="graphics-symbol" aria-label="{udata.get("name", sigla)} ({sigla})"/>'
        )

    # Generate State select dropdown options
    select_options: list[str] = []
    for sigla, udata in sorted(ufs.items(), key=lambda x: x[1].get("name", x[0])):
        selected = ' selected="selected"' if sigla == default_uf else ""
        name = udata.get("name", sigla)
        select_options.append(f'<option value="{sigla}"{selected}>{name} ({sigla})</option>')

    # Generate Audit Table initial rows
    table_rows: list[str] = []
    for sigla, udata in sorted(ufs.items(), key=lambda x: x[1].get("ibge_code", "")):
        code = udata.get("ibge_code", "")
        name = udata.get("name", sigla)
        reg_br = udata.get("region_br", "")
        reg = udata.get("rufas_region", "Southeast")
        reg_col = udata.get("region_color", "#3b82f6")
        bmatch = udata.get("best_match", {})
        us_str = f"{bmatch.get('us_state', '')} ({bmatch.get('us_abbrev', '')})"
        gower = bmatch.get("gower_similarity", 0.0)
        wadj = udata.get("wood_adjustments", {})
        wcurve = udata.get("wood_curve", {})
        active_row_cls = ' class="active-row"' if sigla == default_uf else ""

        table_rows.append(
            f'<tr id="row-{sigla}" data-uf="{sigla}"{active_row_cls} onclick="selectState(\'{sigla}\')">\n'
            f'  <td class="table-num">{code}</td>\n'
            f'  <td><span class="table-badge" style="background:#1e293b; color:#ffffff;">{sigla}</span></td>\n'
            f'  <td><strong>{name}</strong></td>\n'
            f'  <td><span class="table-badge-outline">{reg_br}</span></td>\n'
            f'  <td><span class="table-badge" style="background:{reg_col}; color:#ffffff;">{reg}</span></td>\n'
            f'  <td>{us_str}</td>\n'
            f'  <td class="table-num"><strong>{gower * 100:.1f}%</strong></td>\n'
            f'  <td class="table-num">{wadj.get("l", 0.0):+.2f}</td>\n'
            f'  <td class="table-num">{wadj.get("m", 0.0):+.2f}</td>\n'
            f'  <td class="table-num">{wadj.get("n", 0.0):+.2f}</td>\n'
            f'  <td class="table-num">{wcurve.get("peak_day", 0.0):.1f}</td>\n'
            f'  <td class="table-num"><strong>{wcurve.get("peak_yield", 0.0):.1f}</strong></td>\n'
            f'  <td class="table-num"><strong>{wcurve.get("total_305d", 0.0):,.0f}</strong></td>\n'
            f"</tr>"
        )

    json_payload_str = json.dumps(payload, ensure_ascii=False)

    html = f"""<!DOCTYPE html>
<html lang="pt-BR">
<head>
  <meta charset="UTF-8">
  <meta name="viewport" content="width=device-width, initial-scale=1.0">
  <title>RuFaS Lactation Equivalence — Brasil ↔ EUA</title>
  <!-- Chart.js CDN with graceful inline SVG fallback -->
  <script src="https://cdn.jsdelivr.net/npm/chart.js"></script>
  <style>
{CSS_STYLES}
  </style>
</head>
<body>
  <div class="dashboard-container">
    <!-- Header -->
    <header class="dashboard-header">
      <div class="header-top">
        <span class="brand-badge">RuFaS Platform • Dairy Biophysics</span>
        <span class="header-version">Equivalence Engine v1.0 • Li et al. (2022)</span>
      </div>
      <h1 class="dashboard-title">Equivalência Edafoclimática & Curvas de Lactação (Brasil ↔ EUA)</h1>
      <p class="dashboard-subtitle">
        Mapeamento biofísico determinístico das 27 UFs brasileiras para macrorregiões do modelo RuFaS
        via Coeficiente de Similaridade de Gower com calibração das curvas de lactação de Wood (DIM 1..305).
      </p>

      <!-- Summary KPI Indicators -->
      <div class="summary-kpi-grid">
        <div class="kpi-card">
          <div class="kpi-label">Unidades Federativas</div>
          <div class="kpi-value">{total_ufs} UFs</div>
          <div class="kpi-meta">100% Malha IBGE coberta</div>
        </div>
        <div class="kpi-card">
          <div class="kpi-label">Macrorregiões RuFaS</div>
          <div class="kpi-value">{total_regions} Regiões</div>
          <div class="kpi-meta">Delta, SE, Mtn, S.Plains, Appal.</div>
        </div>
        <div class="kpi-card">
          <div class="kpi-label">Similaridade Média</div>
          <div class="kpi-value">{avg_similarity * 100:.1f}%</div>
          <div class="kpi-meta">Gower Multivariado (6 Fatores)</div>
        </div>
        <div class="kpi-card">
          <div class="kpi-label">Intervalo de Similaridade</div>
          <div class="kpi-value">{min_similarity * 100:.1f}% – {max_similarity * 100:.1f}%</div>
          <div class="kpi-meta">Piso calibrado > 75.0%</div>
        </div>
      </div>
    </header>

    <!-- Main 2-Column Responsive Layout -->
    <main class="dashboard-main">
      <!-- Left Column: Map & Filter Controls -->
      <section class="map-column">
        <div class="dash-card">
          <div class="card-header">
            <div>
              <h2 class="card-title">Mapa Biofísico do Brasil</h2>
              <div class="card-subtitle">Clique em uma UF para inspecionar os parâmetros</div>
            </div>
          </div>
          <div class="card-body">
            <!-- Region Filter Badges -->
            <div class="region-filters">
              <button type="button" class="region-pill active" data-region="ALL" onclick="filterByRegion('ALL')">
                Todas as Regiões
              </button>
              <button type="button" class="region-pill" data-region="Delta" onclick="filterByRegion('Delta')">
                <span class="region-color-dot" style="background: var(--color-delta);"></span> Delta
              </button>
              <button type="button" class="region-pill" data-region="Southeast" onclick="filterByRegion('Southeast')">
                <span class="region-color-dot" style="background: var(--color-southeast);"></span> Southeast
              </button>
              <button type="button" class="region-pill" data-region="Mountain" onclick="filterByRegion('Mountain')">
                <span class="region-color-dot" style="background: var(--color-mountain);"></span> Mountain
              </button>
              <button type="button" class="region-pill" data-region="Southern Plains" onclick="filterByRegion('Southern Plains')">
                <span class="region-color-dot" style="background: var(--color-southern-plains);"></span> Southern Plains
              </button>
              <button type="button" class="region-pill" data-region="Appalachian" onclick="filterByRegion('Appalachian')">
                <span class="region-color-dot" style="background: var(--color-appalachian);"></span> Appalachian
              </button>
            </div>

            <!-- Brazil SVG Map Viewport -->
            <div class="map-viewport">
              <svg id="brazil-map" viewBox="0 0 600 600" width="100%" height="auto" role="img" aria-label="Mapa Vetorial do Brasil com 27 UFs">
                <g id="ufs-group">
                  {' '.join(svg_paths_html)}
                </g>
              </svg>
              <!-- Hover Floating Tooltip -->
              <div id="map-tooltip" class="map-tooltip" style="display: none;"></div>
            </div>

            <div class="map-controls-bar">
              <span>Projeção: Mercator calibrada IBGE</span>
              <button type="button" class="btn-map-control" onclick="resetMapSelection()">Resetar Seleção</button>
            </div>
          </div>
        </div>
      </section>

      <!-- Right Column: Active State, Wood Lactation Curve & Radar -->
      <section class="detail-column">
        <!-- Active UF Identity & Homologue Banner -->
        <div class="dash-card">
          <div class="card-body">
            <div class="active-state-bar">
              <div class="state-select-wrapper">
                <label for="stateSelect" class="state-select-label">Selecionar UF:</label>
                <select id="stateSelect" class="state-select" onchange="selectState(this.value)">
                  {' '.join(select_options)}
                </select>
              </div>
              <div class="active-state-identity" id="active-state-header">
                <span class="state-pill-badge" id="active-uf-badge">{default_uf}</span>
                <span class="state-name-display" id="active-uf-name">{default_state.get('name', 'São Paulo')}</span>
                <span class="region-tag-badge" id="active-uf-region" style="background: {default_state.get('region_color', '#3b82f6')};">
                  {default_state.get('rufas_region', 'Southeast')}
                </span>
              </div>
            </div>

            <!-- US Homologue Card -->
            <div class="homologue-banner" id="homologue-card">
              <div class="homologue-info">
                <span class="homologue-meta">Top-1 Nearest Homologue nos EUA</span>
                <div class="homologue-title" id="homologue-title">
                  {default_match.get('us_state', 'Georgia')} ({default_match.get('us_abbrev', 'GA')})
                </div>
                <div style="font-size: 0.8rem; color: #475569; margin-top: 2px;" id="homologue-region-label">
                  Região Leiteira RuFaS: <strong>{default_match.get('rufas_region', 'Southeast')}</strong>
                </div>
              </div>
              <div class="homologue-scores">
                <div class="score-badge">
                  <span class="score-label">Similaridade Gower</span>
                  <span class="score-val" id="homologue-gower-val">{default_match.get('gower_similarity', 0.892) * 100:.1f}%</span>
                </div>
                <div class="score-badge">
                  <span class="score-label">Similaridade Euclidiana</span>
                  <span class="score-val" id="homologue-euclidean-val" style="color: #64748b; font-size: 1.15rem;">
                    {default_match.get('euclidean_similarity', 0.865) * 100:.1f}%
                  </span>
                </div>
              </div>
            </div>

            <!-- Wood Lactation Curve Card Section -->
            <div>
              <div style="display: flex; align-items: center; justify-content: space-between; margin-bottom: 8px;">
                <h3 class="card-title">Curva de Lactação de Wood (DIM 1..305)</h3>
                <span style="font-size: 0.75rem; color: var(--color-text-muted);">Y(t) = l · t^m · e^(-nt)</span>
              </div>

              <!-- Quick Lactation KPI Strip -->
              <div class="lactation-metrics-strip">
                <div class="lactation-metric-box">
                  <div class="metric-label">Dia do Pico</div>
                  <div class="metric-value-wrap">
                    <span id="peak-day-val" class="metric-num">{default_curve.get('peak_day', 73.16):.1f}</span>
                    <span class="metric-unit">DIM</span>
                  </div>
                </div>
                <div class="lactation-metric-box">
                  <div class="metric-label">Produção de Pico</div>
                  <div class="metric-value-wrap">
                    <span id="peak-yield-val" class="metric-num">{default_curve.get('peak_yield', 27.84):.1f}</span>
                    <span class="metric-unit">kg/dia</span>
                  </div>
                </div>
                <div class="lactation-metric-box">
                  <div class="metric-label">Acumulado 305d</div>
                  <div class="metric-value-wrap">
                    <span id="total-305-val" class="metric-num">{default_curve.get('total_305d', 6842.15):,.0f}</span>
                    <span class="metric-unit">kg</span>
                  </div>
                </div>
                <div class="lactation-metric-box">
                  <div class="metric-label">Ajustes Wood (Δl, Δm, Δn)</div>
                  <div class="metric-value-wrap">
                    <span id="wood-deltas-val" class="metric-num" style="font-size: 1.05rem;">
                      {default_wood_adj.get('l', 0.0):+.1f}, {default_wood_adj.get('m', 0.0):+.1f}, {default_wood_adj.get('n', 0.0):+.1f}
                    </span>
                  </div>
                </div>
              </div>

              <!-- Canvas / Fallback Container -->
              <div class="chart-container-box">
                <canvas id="lactationChart"></canvas>
                <div id="lactation-svg-fallback" class="svg-chart-fallback" style="display: none;"></div>
              </div>
            </div>

            <!-- Radar Chart & Biophysical Factors -->
            <div style="margin-top: 24px; padding-top: 20px; border-top: 1px solid var(--color-border);">
              <div style="display: flex; align-items: center; justify-content: space-between; margin-bottom: 8px;">
                <h3 class="card-title">Decomposição Edafoclimática dos 6 Fatores</h3>
                <span style="font-size: 0.75rem; color: var(--color-text-muted);">Comparativo Normalizado [0, 1]</span>
              </div>
              <div class="chart-container-box" style="height: 290px;">
                <canvas id="radarChart"></canvas>
              </div>

              <!-- Biophysical Raw Metrics Grid -->
              <div class="biophysical-grid" id="biophysical-metrics-grid">
                <div class="bio-metric-cell">
                  <div class="bio-metric-label">Temp. Média</div>
                  <div class="bio-metric-val" id="bio-temp-val">{default_raw.get('temp_mean_c', 21.2):.1f} °C</div>
                </div>
                <div class="bio-metric-cell">
                  <div class="bio-metric-label">Sazonalidade Térm.</div>
                  <div class="bio-metric-val" id="bio-season-val">{default_raw.get('temp_seasonality_c', 4.5):.1f} °C</div>
                </div>
                <div class="bio-metric-cell">
                  <div class="bio-metric-label">Precip. Anual</div>
                  <div class="bio-metric-val" id="bio-precip-val">{default_raw.get('precip_annual_mm', 1380.0):.0f} mm</div>
                </div>
                <div class="bio-metric-cell">
                  <div class="bio-metric-label">Umidade Rel.</div>
                  <div class="bio-metric-val" id="bio-rh-val">{default_raw.get('rh_mean_pct', 74.0):.1f} %</div>
                </div>
                <div class="bio-metric-cell">
                  <div class="bio-metric-label">Argila Solo</div>
                  <div class="bio-metric-val" id="bio-clay-val">{default_raw.get('soil_clay_pct', 32.5):.1f} %</div>
                </div>
                <div class="bio-metric-cell">
                  <div class="bio-metric-label">PAWC (Água Solo)</div>
                  <div class="bio-metric-val" id="bio-pawc-val">{default_raw.get('pawc', 0.12):.3f} cm³/cm³</div>
                </div>
                <div class="bio-metric-cell">
                  <div class="bio-metric-label">pH em Água</div>
                  <div class="bio-metric-val" id="bio-ph-val">{default_raw.get('soil_ph', 5.4):.2f}</div>
                </div>
                <div class="bio-metric-cell">
                  <div class="bio-metric-label">Carbono Orgânico</div>
                  <div class="bio-metric-val" id="bio-soc-val">{default_raw.get('soil_soc_g_kg', 18.2):.1f} g/kg</div>
                </div>
              </div>
            </div>

          </div>
        </div>
      </section>
    </main>

    <!-- Bottom View: 27-UF Responsive Audit Table -->
    <section class="table-section">
      <div class="table-card">
        <div class="table-toolbar">
          <div>
            <h2 class="card-title">Tabela de Auditoria das 27 UFs & Equivalências RuFaS</h2>
            <div class="card-subtitle">Clique em qualquer linha para sincronizar o mapa e os gráficos</div>
          </div>
          <div class="table-search-box">
            <span class="search-icon">🔍</span>
            <input type="text" id="tableSearchInput" class="table-search"
                   placeholder="Filtrar por UF, estado, região ou homólogo americano..."
                   onkeyup="filterTable(this.value)">
          </div>
          <div class="table-actions">
            <button type="button" class="btn-utility" id="btnCopyCsv" onclick="copyTableAsCSV()">Copiar CSV</button>
            <button type="button" class="btn-utility" id="btnExportJson" onclick="exportJsonData()">Exportar JSON</button>
          </div>
        </div>

        <div class="table-scroll-wrapper">
          <table id="auditTable" class="audit-table">
            <thead>
              <tr>
                <th data-sort-key="ibge_code" onclick="sortTable('ibge_code', true)" class="table-num">IBGE ⇕</th>
                <th data-sort-key="uf" onclick="sortTable('uf', false)">UF ⇕</th>
                <th data-sort-key="name" onclick="sortTable('name', false)">Estado Brasileiro ⇕</th>
                <th data-sort-key="region_br" onclick="sortTable('region_br', false)">Região BR ⇕</th>
                <th data-sort-key="rufas_region" onclick="sortTable('rufas_region', false)">Região RuFaS ⇕</th>
                <th data-sort-key="us_homologue" onclick="sortTable('us_homologue', false)">Top-1 Homólogo EUA ⇕</th>
                <th data-sort-key="gower_similarity" onclick="sortTable('gower_similarity', true)" class="table-num">Gower Sim. ⇕</th>
                <th data-sort-key="delta_l" onclick="sortTable('delta_l', true)" class="table-num">Δl ⇕</th>
                <th data-sort-key="delta_m" onclick="sortTable('delta_m', true)" class="table-num">Δm ⇕</th>
                <th data-sort-key="delta_n" onclick="sortTable('delta_n', true)" class="table-num">Δn ⇕</th>
                <th data-sort-key="peak_day" onclick="sortTable('peak_day', true)" class="table-num">Pico (DIM) ⇕</th>
                <th data-sort-key="peak_yield" onclick="sortTable('peak_yield', true)" class="table-num">Pico (kg/d) ⇕</th>
                <th data-sort-key="total_305d" onclick="sortTable('total_305d', true)" class="table-num">Total 305d (kg) ⇕</th>
              </tr>
            </thead>
            <tbody id="auditTableBody">
              {' '.join(table_rows)}
            </tbody>
          </table>
        </div>
      </div>
    </section>

    <!-- Footer -->
    <footer class="dashboard-footer">
      <div>RuFaS Brazilian Lactation Equivalence Engine • Documento estático autocontido</div>
      <div>Autores: Thiago & RuFaS Modeling Consortium • Referência: Li et al. (2022/2023)</div>
    </footer>
  </div>

  <!-- Embedded Deterministic Dashboard Data Payload -->
  <script id="dashboard-data" type="application/json">
{json_payload_str}
  </script>

  <!-- Client-side Interactive Application Controller Skeleton -->
  <script>
    /* RuFaS Lactation Equivalence Dashboard Client Controller */
    let dashboardPayload = null;
    let currentSelectedUf = "{default_uf}";
    let lactationChartInstance = null;
    let radarChartInstance = null;
    let sortState = {{ column: null, ascending: true }};

    function initDashboard() {{
      const dataElem = document.getElementById("dashboard-data");
      if (dataElem) {{
        try {{
          dashboardPayload = JSON.parse(dataElem.textContent);
        }} catch (e) {{
          console.error("Failed to parse dashboard data payload:", e);
        }}
      }}
      attachMapEvents();
      if (typeof updateLactationChart === "function") {{
        updateLactationChart(getSelectedUfData());
      }}
      if (typeof updateRadarChart === "function") {{
        updateRadarChart(getSelectedUfData());
      }}
    }}

    function getSelectedUfData() {{
      if (!dashboardPayload || !dashboardPayload.ufs) return null;
      return dashboardPayload.ufs[currentSelectedUf] || null;
    }}

    function selectState(uf) {{
      if (!uf || !dashboardPayload || !dashboardPayload.ufs[uf]) return;
      currentSelectedUf = uf;

      // Update Map Selection
      document.querySelectorAll(".uf-path").forEach(p => p.classList.remove("active"));
      const mapPath = document.getElementById(`uf-${{uf}}`);
      if (mapPath) mapPath.classList.add("active");

      // Update State Select Dropdown
      const selectElem = document.getElementById("stateSelect");
      if (selectElem && selectElem.value !== uf) {{
        selectElem.value = uf;
      }}

      // Update Table Active Row
      document.querySelectorAll("#auditTableBody tr").forEach(r => r.classList.remove("active-row"));
      const tableRow = document.getElementById(`row-${{uf}}`);
      if (tableRow) tableRow.classList.add("active-row");

      // Update Information Views & Charts
      const ufData = dashboardPayload.ufs[uf];
      updateActiveStateDetails(ufData);
      if (typeof updateLactationChart === "function") {{
        updateLactationChart(ufData);
      }}
      if (typeof updateRadarChart === "function") {{
        updateRadarChart(ufData);
      }}
    }}

    function updateActiveStateDetails(ufData) {{
      if (!ufData) return;
      const bmatch = ufData.best_match || {{}};
      const wcurve = ufData.wood_curve || {{}};
      const wadj = ufData.wood_adjustments || {{}};
      const raw = ufData.raw_profile || {{}};

      const badge = document.getElementById("active-uf-badge");
      if (badge) badge.textContent = ufData.uf;

      const nameElem = document.getElementById("active-uf-name");
      if (nameElem) nameElem.textContent = ufData.name;

      const regElem = document.getElementById("active-uf-region");
      if (regElem) {{
        regElem.textContent = ufData.rufas_region;
        regElem.style.background = ufData.region_color || "#3b82f6";
      }}

      const titleElem = document.getElementById("homologue-title");
      if (titleElem) {{
        titleElem.textContent = `${{bmatch.us_state || ''}} (${{bmatch.us_abbrev || ''}})`;
      }}

      const regLabel = document.getElementById("homologue-region-label");
      if (regLabel) {{
        regLabel.innerHTML = `Região Leiteira RuFaS: <strong>${{bmatch.rufas_region || ''}}</strong>`;
      }}

      const gowerVal = document.getElementById("homologue-gower-val");
      if (gowerVal && bmatch.gower_similarity !== undefined) {{
        gowerVal.textContent = `${{(bmatch.gower_similarity * 100).toFixed(1)}}%`;
      }}

      const euclidVal = document.getElementById("homologue-euclidean-val");
      if (euclidVal && bmatch.euclidean_similarity !== undefined) {{
        euclidVal.textContent = `${{(bmatch.euclidean_similarity * 100).toFixed(1)}}%`;
      }}

      const peakDay = document.getElementById("peak-day-val");
      if (peakDay && wcurve.peak_day !== undefined) peakDay.textContent = wcurve.peak_day.toFixed(1);

      const peakYield = document.getElementById("peak-yield-val");
      if (peakYield && wcurve.peak_yield !== undefined) peakYield.textContent = wcurve.peak_yield.toFixed(1);

      const total305 = document.getElementById("total-305-val");
      if (total305 && wcurve.total_305d !== undefined) total305.textContent = Math.round(wcurve.total_305d).toLocaleString('pt-BR');

      const woodDeltas = document.getElementById("wood-deltas-val");
      if (woodDeltas) {{
        const dl = (wadj.l || 0).toFixed(1);
        const dm = (wadj.m || 0).toFixed(1);
        const dn = (wadj.n || 0).toFixed(1);
        woodDeltas.textContent = `${{dl > 0 ? '+' : ''}}${{dl}}, ${{dm > 0 ? '+' : ''}}${{dm}}, ${{dn > 0 ? '+' : ''}}${{dn}}`;
      }}

      // Biophysical raw numbers
      const setVal = (id, val, suffix = '') => {{
        const el = document.getElementById(id);
        if (el) el.textContent = val !== undefined && val !== null ? `${{val}}${{suffix}}` : '--';
      }};
      setVal("bio-temp-val", raw.temp_mean_c?.toFixed(1), " °C");
      setVal("bio-season-val", raw.temp_seasonality_c?.toFixed(1), " °C");
      setVal("bio-precip-val", raw.precip_annual_mm ? Math.round(raw.precip_annual_mm) : null, " mm");
      setVal("bio-rh-val", raw.rh_mean_pct?.toFixed(1), " %");
      setVal("bio-clay-val", raw.soil_clay_pct?.toFixed(1), " %");
      const pawcNum = raw.pawc !== undefined && raw.pawc !== null
        ? raw.pawc.toFixed(3)
        : (ufData.normalized_vector && ufData.normalized_vector.soil_pawc_index !== undefined
            ? (ufData.normalized_vector.soil_pawc_index * 0.18 + 0.04).toFixed(3)
            : null);
      setVal("bio-pawc-val", pawcNum, " cm³/cm³");
      setVal("bio-ph-val", raw.soil_ph?.toFixed(2));
      setVal("bio-soc-val", raw.soil_soc_g_kg?.toFixed(1), " g/kg");
    }}

    function attachMapEvents() {{
      const map = document.getElementById("brazil-map");
      if (!map) return;

      map.querySelectorAll(".uf-path").forEach(path => {{
        path.addEventListener("mouseenter", (e) => {{
          const uf = path.getAttribute("data-uf");
          path.classList.add("hovered");
          showTooltip(e, uf);
        }});
        path.addEventListener("mousemove", (e) => {{
          const uf = path.getAttribute("data-uf");
          showTooltip(e, uf);
        }});
        path.addEventListener("mouseleave", () => {{
          path.classList.remove("hovered");
          hideTooltip();
        }});
        path.addEventListener("click", () => {{
          const uf = path.getAttribute("data-uf");
          if (uf) selectState(uf);
        }});
      }});
    }}

    function showTooltip(event, uf) {{
      const tooltip = document.getElementById("map-tooltip");
      if (!tooltip || !dashboardPayload || !dashboardPayload.ufs[uf]) return;
      const udata = dashboardPayload.ufs[uf];
      const bmatch = udata.best_match || {{}};
      const wcurve = udata.wood_curve || {{}};

      tooltip.innerHTML = `
        <div class="tooltip-header">
          <span>${{udata.name}} (${{uf}})</span>
          <span style="font-size:0.75rem; color:${{udata.region_color || '#38bdf8'}};">${{udata.rufas_region}}</span>
        </div>
        <div class="tooltip-row"><span>Homólogo EUA:</span> <strong>${{bmatch.us_state || ''}} (${{bmatch.us_abbrev || ''}})</strong></div>
        <div class="tooltip-row"><span>Similaridade Gower:</span> <strong>${{((bmatch.gower_similarity || 0) * 100).toFixed(1)}}%</strong></div>
        <div class="tooltip-row"><span>Pico de Produção:</span> <strong>${{wcurve.peak_yield?.toFixed(1)}} kg/d (DIM ${{wcurve.peak_day?.toFixed(0)}})</strong></div>
      `;

      tooltip.style.display = "block";
      const viewport = document.querySelector(".map-viewport");
      if (viewport) {{
        const rect = viewport.getBoundingClientRect();
        const x = event.clientX - rect.left + 14;
        const y = event.clientY - rect.top + 14;
        tooltip.style.left = `${{Math.min(x, rect.width - 240)}}px`;
        tooltip.style.top = `${{Math.min(y, rect.height - 120)}}px`;
      }}
    }}

    function hideTooltip() {{
      const tooltip = document.getElementById("map-tooltip");
      if (tooltip) tooltip.style.display = "none";
    }}

    function filterByRegion(region) {{
      document.querySelectorAll(".region-pill").forEach(btn => {{
        if (btn.getAttribute("data-region") === region) {{
          btn.classList.add("active");
        }} else {{
          btn.classList.remove("active");
        }}
      }});

      document.querySelectorAll(".uf-path").forEach(path => {{
        const pathReg = path.getAttribute("data-region");
        if (region === "ALL" || pathReg === region) {{
          path.classList.remove("faded");
        }} else {{
          path.classList.add("faded");
        }}
      }});
    }}

    function resetMapSelection() {{
      filterByRegion("ALL");
      selectState("{default_uf}");
    }}

    function updateLactationChart(ufData) {{
      if (!ufData) return;
      const canvas = document.getElementById("lactationChart");
      const fallback = document.getElementById("lactation-svg-fallback");
      if (!canvas) return;

      const wcurve = ufData.wood_curve || {{}};
      const wadj = ufData.wood_adjustments || {{}};

      // Update Quick Lactation KPI Strip
      const peakDayElem = document.getElementById("peak-day-val");
      if (peakDayElem && wcurve.peak_day !== undefined) {{
        peakDayElem.textContent = wcurve.peak_day.toFixed(1);
      }}

      const peakYieldElem = document.getElementById("peak-yield-val");
      if (peakYieldElem && wcurve.peak_yield !== undefined) {{
        peakYieldElem.textContent = wcurve.peak_yield.toFixed(1);
      }}

      const total305Elem = document.getElementById("total-305-val");
      if (total305Elem && wcurve.total_305d !== undefined) {{
        total305Elem.textContent = Math.round(wcurve.total_305d).toLocaleString('pt-BR');
      }}

      const woodDeltasElem = document.getElementById("wood-deltas-val");
      if (woodDeltasElem) {{
        const dl = (wadj.l || 0).toFixed(1);
        const dm = (wadj.m || 0).toFixed(1);
        const dn = (wadj.n || 0).toFixed(1);
        woodDeltasElem.textContent = `${{dl > 0 ? '+' : ''}}${{dl}}, ${{dm > 0 ? '+' : ''}}${{dm}}, ${{dn > 0 ? '+' : ''}}${{dn}}`;
      }}

      // Graceful fallback if Chart.js CDN is unavailable
      if (typeof Chart === "undefined") {{
        canvas.style.display = "none";
        if (fallback) {{
          fallback.style.display = "block";
          renderSvgFallbackLactation(fallback, ufData);
        }}
        return;
      }}

      canvas.style.display = "block";
      if (fallback) fallback.style.display = "none";

      const neutralYields = (dashboardPayload && dashboardPayload.summary && dashboardPayload.summary.neutral_wood_curve && dashboardPayload.summary.neutral_wood_curve.yield_kg)
        ? dashboardPayload.summary.neutral_wood_curve.yield_kg
        : Array.from({{ length: 305 }}, (_, i) => {{
            const t = i + 1;
            return +(19.90 * Math.pow(t, 0.247) * Math.exp(-0.003376 * t)).toFixed(4);
          }});

      const ufYields = wcurve.yield_kg || [];
      const dimLabels = wcurve.dim || Array.from({{ length: 305 }}, (_, i) => i + 1);
      const ufColor = ufData.region_color || "#2563eb";
      const peakDayNum = Math.round(wcurve.peak_day || 73);

      const datasets = [
        {{
          label: `${{ufData.name}} (${{ufData.uf}}) — Curva Ajustada`,
          data: ufYields,
          borderColor: ufColor,
          backgroundColor: ufColor + "18",
          borderWidth: 2.5,
          fill: true,
          tension: 0.3,
          pointRadius: (ctx) => (ctx.dataIndex === peakDayNum - 1 ? 6 : 0),
          pointHoverRadius: 7,
          pointBackgroundColor: (ctx) => (ctx.dataIndex === peakDayNum - 1 ? "#ef4444" : ufColor),
          pointBorderColor: "#ffffff",
          pointBorderWidth: 2,
        }},
        {{
          label: "Linha de Base RuFaS (l=19.9, m=0.247, n=0.003376)",
          data: neutralYields,
          borderColor: "#94a3b8",
          borderWidth: 1.75,
          borderDash: [5, 5],
          fill: false,
          tension: 0.3,
          pointRadius: 0,
          pointHoverRadius: 4,
        }}
      ];

      if (lactationChartInstance) {{
        lactationChartInstance.data.labels = dimLabels;
        lactationChartInstance.data.datasets = datasets;
        lactationChartInstance.update();
      }} else {{
        const ctx = canvas.getContext("2d");
        lactationChartInstance = new Chart(ctx, {{
          type: "line",
          data: {{
            labels: dimLabels,
            datasets: datasets,
          }},
          options: {{
            responsive: true,
            maintainAspectRatio: false,
            animation: {{ duration: 350 }},
            interaction: {{
              mode: "index",
              intersect: false,
            }},
            plugins: {{
              legend: {{
                position: "top",
                labels: {{
                  boxWidth: 14,
                  font: {{ size: 11, family: "var(--font-sans)" }},
                  color: "#334155",
                }},
              }},
              tooltip: {{
                callbacks: {{
                  title: (items) => `DIM ${{items[0].label}}`,
                  label: (item) => `${{item.dataset.label}}: ${{item.parsed.y.toFixed(2)}} kg/dia`,
                }},
              }},
            }},
            scales: {{
              x: {{
                title: {{
                  display: true,
                  text: "Dias em Lactação (DIM)",
                  font: {{ size: 11, weight: "bold" }},
                  color: "#64748b",
                }},
                grid: {{ color: "rgba(226, 232, 240, 0.6)" }},
                ticks: {{ maxTicksLimit: 12, color: "#64748b" }},
              }},
              y: {{
                title: {{
                  display: true,
                  text: "Produção (kg/vaca/dia)",
                  font: {{ size: 11, weight: "bold" }},
                  color: "#64748b",
                }},
                grid: {{ color: "rgba(226, 232, 240, 0.6)" }},
                suggestedMin: 12,
                suggestedMax: 34,
                ticks: {{ color: "#64748b" }},
              }},
            }},
          }},
        }});
      }}
    }}

    function updateRadarChart(ufData) {{
      if (!ufData) return;
      const canvas = document.getElementById("radarChart");
      if (!canvas || typeof Chart === "undefined") return;

      const ufColor = ufData.region_color || "#3b82f6";
      const bmatch = ufData.best_match || {{}};
      const ufVec = ufData.normalized_vector || {{}};
      const usVec = bmatch.normalized_vector || {{}};

      const labels = [
        "THI",
        "Sazonalidade Térmica",
        "Índice Hídrico",
        "Argila Solo",
        "PAWC",
        "Fertilidade Química"
      ];

      const ufScores = [
        ufVec.thi_index ?? 0.5,
        ufVec.temp_seasonality ?? 0.5,
        ufVec.moisture_index ?? 0.5,
        ufVec.soil_clay_fraction ?? 0.5,
        ufVec.soil_pawc_index ?? 0.5,
        ufVec.soil_fertility_index ?? 0.5,
      ];

      const usScores = [
        usVec.thi_index ?? 0.5,
        usVec.temp_seasonality ?? 0.5,
        usVec.moisture_index ?? 0.5,
        usVec.soil_clay_fraction ?? 0.5,
        usVec.soil_pawc_index ?? 0.5,
        usVec.soil_fertility_index ?? 0.5,
      ];

      const hexToRgba = (hex, alpha) => {{
        const clean = (hex || "#3b82f6").replace("#", "");
        const r = parseInt(clean.substring(0, 2), 16) || 59;
        const g = parseInt(clean.substring(2, 4), 16) || 130;
        const b = parseInt(clean.substring(4, 6), 16) || 246;
        return `rgba(${{r}}, ${{g}}, ${{b}}, ${{alpha}})`;
      }};

      const datasets = [
        {{
          label: `${{ufData.name}} (${{ufData.uf}})`,
          data: ufScores,
          backgroundColor: hexToRgba(ufColor, 0.25),
          borderColor: ufColor,
          pointBackgroundColor: ufColor,
          pointBorderColor: "#ffffff",
          pointHoverBackgroundColor: "#ffffff",
          pointHoverBorderColor: ufColor,
          pointRadius: 4,
          borderWidth: 2,
        }},
        {{
          label: `Homólogo: ${{bmatch.us_state || ''}} (${{bmatch.us_abbrev || ''}})`,
          data: usScores,
          backgroundColor: "rgba(100, 116, 139, 0.15)",
          borderColor: "#64748b",
          borderDash: [5, 5],
          pointBackgroundColor: "#64748b",
          pointBorderColor: "#ffffff",
          pointHoverBackgroundColor: "#ffffff",
          pointHoverBorderColor: "#64748b",
          pointRadius: 3,
          borderWidth: 2,
        }}
      ];

      if (radarChartInstance) {{
        radarChartInstance.data.labels = labels;
        radarChartInstance.data.datasets = datasets;
        radarChartInstance.update();
      }} else {{
        const ctx = canvas.getContext("2d");
        radarChartInstance = new Chart(ctx, {{
          type: "radar",
          data: {{
            labels: labels,
            datasets: datasets,
          }},
          options: {{
            responsive: true,
            maintainAspectRatio: false,
            animation: {{ duration: 350 }},
            scales: {{
              r: {{
                min: 0.0,
                max: 1.0,
                ticks: {{
                  stepSize: 0.2,
                  display: true,
                  backdropColor: "transparent",
                  color: "#94a3b8",
                  font: {{ size: 9 }},
                }},
                pointLabels: {{
                  font: {{ size: 10, weight: "bold", family: "var(--font-sans)" }},
                  color: "#334155",
                }},
                grid: {{
                  color: "rgba(203, 213, 225, 0.6)",
                }},
                angleLines: {{
                  color: "rgba(203, 213, 225, 0.6)",
                }},
              }},
            }},
            plugins: {{
              legend: {{
                position: "top",
                labels: {{
                  boxWidth: 12,
                  font: {{ size: 11, family: "var(--font-sans)" }},
                  color: "#334155",
                }},
              }},
              tooltip: {{
                callbacks: {{
                  label: (item) => `${{item.dataset.label}}: ${{(item.raw * 100).toFixed(1)}}% (${{item.raw.toFixed(3)}})`,
                }},
              }},
            }},
          }},
        }});
      }}
    }}

    function renderSvgFallbackLactation(container, ufData) {{
      const target = container || document.getElementById("lactation-svg-fallback");
      if (!target || !ufData) return;
      target.style.display = "block";

      const wcurve = ufData.wood_curve || {{}};
      const yields = wcurve.yield_kg || [];
      const neutralYields = (dashboardPayload && dashboardPayload.summary && dashboardPayload.summary.neutral_wood_curve)
        ? dashboardPayload.summary.neutral_wood_curve.yield_kg
        : [];
      if (!yields.length) return;

      const w = 560;
      const h = 250;
      const padLeft = 45;
      const padRight = 20;
      const padTop = 30;
      const padBottom = 35;
      const plotW = w - padLeft - padRight;
      const plotH = h - padTop - padBottom;

      const minY = 10;
      const maxY = 35;

      const getX = (dim) => padLeft + ((dim - 1) / 304) * plotW;
      const getY = (yVal) => padTop + plotH - ((yVal - minY) / (maxY - minY)) * plotH;

      const ptsUf = yields.map((val, idx) => `${{getX(idx + 1).toFixed(1)}},${{getY(val).toFixed(1)}}`).join(" ");
      const ptsBase = neutralYields.length
        ? neutralYields.map((val, idx) => `${{getX(idx + 1).toFixed(1)}},${{getY(val).toFixed(1)}}`).join(" ")
        : "";

      const peakDay = wcurve.peak_day || 73;
      const peakYield = wcurve.peak_yield || 28;
      const peakX = getX(peakDay).toFixed(1);
      const peakY = getY(peakYield).toFixed(1);
      const ufColor = ufData.region_color || "#2563eb";

      let gridLines = "";
      [50, 100, 150, 200, 250, 300].forEach(d => {{
        const gx = getX(d).toFixed(1);
        gridLines += `<line x1="${{gx}}" y1="${{padTop}}" x2="${{gx}}" y2="${{padTop + plotH}}" stroke="#e2e8f0" stroke-width="1"/>`;
        gridLines += `<text x="${{gx}}" y="${{padTop + plotH + 16}}" font-size="10" fill="#94a3b8" text-anchor="middle">${{d}}</text>`;
      }});

      [15, 20, 25, 30].forEach(yv => {{
        const gy = getY(yv).toFixed(1);
        gridLines += `<line x1="${{padLeft}}" y1="${{gy}}" x2="${{padLeft + plotW}}" y2="${{gy}}" stroke="#e2e8f0" stroke-width="1"/>`;
        gridLines += `<text x="${{padLeft - 6}}" y="${{parseFloat(gy) + 3}}" font-size="10" fill="#94a3b8" text-anchor="end">${{yv}}</text>`;
      }});

      target.innerHTML = `
        <svg viewBox="0 0 ${{w}} ${{h}}" style="width: 100%; height: 100%; display: block;" role="img" aria-label="Curva de lactação de Wood">
          <rect x="${{padLeft}}" y="${{padTop}}" width="${{plotW}}" height="${{plotH}}" fill="#ffffff" rx="4"/>
          ${{gridLines}}
          ${{ptsBase ? `<polyline points="${{ptsBase}}" fill="none" stroke="#94a3b8" stroke-width="1.8" stroke-dasharray="4,4"/>` : ''}}
          <polyline points="${{ptsUf}}" fill="none" stroke="${{ufColor}}" stroke-width="2.5"/>
          <circle cx="${{peakX}}" cy="${{peakY}}" r="5" fill="#ef4444" stroke="#ffffff" stroke-width="2"/>
          <text x="${{peakX}}" y="${{peakY - 8}}" text-anchor="middle" font-size="10" font-weight="bold" fill="#ef4444">
            Pico: ${{peakYield.toFixed(1)}} kg/d (${{peakDay.toFixed(0)}} DIM)
          </text>
          <line x1="${{padLeft}}" y1="${{padTop + plotH}}" x2="${{padLeft + plotW}}" y2="${{padTop + plotH}}" stroke="#cbd5e1" stroke-width="1.5"/>
          <line x1="${{padLeft}}" y1="${{padTop}}" x2="${{padLeft}}" y2="${{padTop + plotH}}" stroke="#cbd5e1" stroke-width="1.5"/>
          <text x="${{padLeft + plotW / 2}}" y="${{h - 4}}" text-anchor="middle" font-size="10" font-weight="bold" fill="#64748b">
            Dias em Lactação (DIM)
          </text>
          <text x="12" y="${{padTop + plotH / 2}}" text-anchor="middle" font-size="10" font-weight="bold" fill="#64748b" transform="rotate(-90 12 ${{padTop + plotH / 2}})">
            kg/dia
          </text>
        </svg>
      `;
    }}

    function filterTable(query) {{
      const q = (query || "").toLowerCase().trim();
      const rows = document.querySelectorAll("#auditTableBody tr");
      rows.forEach(row => {{
        const text = (row.textContent || row.innerText || "").toLowerCase();
        row.style.display = text.includes(q) ? "" : "none";
      }});
    }}

    function sortTable(columnKey, isNumeric) {{
      const tbody = document.getElementById("auditTableBody");
      if (!tbody) return;
      const rows = Array.from(tbody.querySelectorAll("tr"));

      const isAsc = sortState.column === columnKey ? !sortState.ascending : true;
      sortState = {{ column: columnKey, ascending: isAsc }};

      document.querySelectorAll("#auditTable thead th").forEach(th => {{
        th.classList.remove("sort-asc", "sort-desc");
        if (th.getAttribute("data-sort-key") === columnKey) {{
          th.classList.add(isAsc ? "sort-asc" : "sort-desc");
        }}
      }});

      rows.sort((a, b) => {{
        const ufA = a.getAttribute("data-uf");
        const ufB = b.getAttribute("data-uf");
        const dataA = dashboardPayload.ufs[ufA];
        const dataB = dashboardPayload.ufs[ufB];

        let valA, valB;
        if (columnKey === "ibge_code") {{
          valA = parseInt(dataA.ibge_code, 10);
          valB = parseInt(dataB.ibge_code, 10);
        }} else if (columnKey === "uf") {{
          valA = dataA.uf;
          valB = dataB.uf;
        }} else if (columnKey === "name") {{
          valA = dataA.name;
          valB = dataB.name;
        }} else if (columnKey === "region_br") {{
          valA = dataA.region_br || "";
          valB = dataB.region_br || "";
        }} else if (columnKey === "rufas_region") {{
          valA = dataA.rufas_region;
          valB = dataB.rufas_region;
        }} else if (columnKey === "us_homologue") {{
          valA = dataA.best_match?.us_state || "";
          valB = dataB.best_match?.us_state || "";
        }} else if (columnKey === "gower_similarity") {{
          valA = dataA.best_match?.gower_similarity || 0;
          valB = dataB.best_match?.gower_similarity || 0;
        }} else if (columnKey === "delta_l") {{
          valA = dataA.wood_adjustments?.l || 0;
          valB = dataB.wood_adjustments?.l || 0;
        }} else if (columnKey === "delta_m") {{
          valA = dataA.wood_adjustments?.m || 0;
          valB = dataB.wood_adjustments?.m || 0;
        }} else if (columnKey === "delta_n") {{
          valA = dataA.wood_adjustments?.n || 0;
          valB = dataB.wood_adjustments?.n || 0;
        }} else if (columnKey === "peak_day") {{
          valA = dataA.wood_curve?.peak_day || 0;
          valB = dataB.wood_curve?.peak_day || 0;
        }} else if (columnKey === "peak_yield") {{
          valA = dataA.wood_curve?.peak_yield || 0;
          valB = dataB.wood_curve?.peak_yield || 0;
        }} else if (columnKey === "total_305d") {{
          valA = dataA.wood_curve?.total_305d || 0;
          valB = dataB.wood_curve?.total_305d || 0;
        }}

        if (typeof valA === "string" && typeof valB === "string") {{
          const cmp = valA.localeCompare(valB, "pt-BR");
          return isAsc ? cmp : -cmp;
        }}

        if (valA < valB) return isAsc ? -1 : 1;
        if (valA > valB) return isAsc ? 1 : -1;
        return 0;
      }});

      rows.forEach(r => tbody.appendChild(r));
    }}

    function copyTableAsCSV() {{
      if (!dashboardPayload || !dashboardPayload.ufs) return;
      const headers = ["IBGE", "UF", "Estado", "Região_BR", "Região_RuFaS", "Homólogo_EUA", "Sigla_EUA", "Similaridade_Gower", "Delta_L", "Delta_M", "Delta_N", "Pico_DIM", "Pico_kg_d", "Total_305d_kg"];
      const lines = [headers.join(",")];

      Object.values(dashboardPayload.ufs).forEach(u => {{
        const b = u.best_match || {{}};
        const w = u.wood_curve || {{}};
        const a = u.wood_adjustments || {{}};
        lines.push([
          u.ibge_code,
          u.uf,
          `"${{u.name}}"`,
          `"${{u.region_br || ''}}"`,
          u.rufas_region,
          `"${{b.us_state || ''}}"`,
          b.us_abbrev || '',
          (b.gower_similarity || 0).toFixed(4),
          a.l || 0,
          a.m || 0,
          a.n || 0,
          w.peak_day || 0,
          w.peak_yield || 0,
          w.total_305d || 0
        ].join(","));
      }});

      if (navigator.clipboard && navigator.clipboard.writeText) {{
        navigator.clipboard.writeText(lines.join("\\n")).then(() => {{
          alert("Tabela copiada como CSV para a área de transferência!");
        }}).catch(() => {{
          alert("Não foi possível copiar automaticamente. Use a visualização em tela.");
        }});
      }} else {{
        const textarea = document.createElement("textarea");
        textarea.value = lines.join("\\n");
        document.body.appendChild(textarea);
        textarea.select();
        try {{
          document.execCommand("copy");
          alert("Tabela copiada como CSV para a área de transferência!");
        }} catch (e) {{
          alert("Não foi possível copiar automaticamente.");
        }}
        document.body.removeChild(textarea);
      }}
    }}
    const copyTableAsCsv = copyTableAsCSV;

    function exportJsonData() {{
      if (!dashboardPayload) return;
      const str = "data:text/json;charset=utf-8," + encodeURIComponent(JSON.stringify(dashboardPayload, null, 2));
      const downloadAnchor = document.createElement("a");
      downloadAnchor.setAttribute("href", str);
      downloadAnchor.setAttribute("download", "rufas_brazil_lactation_equivalence.json");
      document.body.appendChild(downloadAnchor);
      downloadAnchor.click();
      downloadAnchor.remove();
    }}
    const exportDataJson = exportJsonData;

    // Auto initialize on DOM ready
    if (document.readyState === "loading") {{
      document.addEventListener("DOMContentLoaded", initDashboard);
    }} else {{
      initDashboard();
    }}
  </script>
</body>
</html>
"""
    return html
