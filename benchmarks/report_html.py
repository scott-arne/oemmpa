"""Render the staged benchmark as a single, self-contained interactive HTML page.

``render_html(records, meta)`` returns one HTML document with **no external
dependencies** — inline CSS, inline JavaScript, and hand-drawn inline SVG charts,
with the benchmark data embedded as JSON. This matches the repository's offline /
corporate-proxy posture (no CDN fetches) and makes the report a single portable
artifact.

Design follows the ``dataviz`` skill: a validated categorical palette (blue / aqua
/ yellow for OEMMPA / mmpdb / RDKit), an ordinal blue ramp for the ordered pipeline
stages, status colors reserved for verdicts, direct labels plus a table view to
satisfy the light-mode relief rule, hover tooltips on every chart, and a selected
dark mode via ``prefers-color-scheme``.
"""

from __future__ import annotations

import json

# The categorical palette validated with the dataviz skill's validate_palette.js
# (light: CVD worst-adjacent dE 47.2, contrast relief handled by direct labels +
# table view; dark: all checks pass).
_CSS = """
:root {
  color-scheme: light dark;
  --page: #f9f9f7;
  --surface-1: #fcfcfb;
  --text-primary: #0b0b0b;
  --text-secondary: #52514e;
  --text-muted: #898781;
  --grid: #e1e0d9;
  --axis: #c3c2b7;
  --border: rgba(11,11,11,0.10);
  --series-oemmpa: #2a78d6;
  --series-mmpdb: #1baf7a;
  --series-rdkit: #eda100;
  --good: #0ca30c;
  --warning: #fab219;
  --serious: #ec835a;
  --critical: #d03b3b;
  --shadow: 0 1px 2px rgba(11,11,11,0.06), 0 2px 8px rgba(11,11,11,0.04);
}
@media (prefers-color-scheme: dark) {
  :root {
    --page: #0d0d0d;
    --surface-1: #1a1a19;
    --text-primary: #ffffff;
    --text-secondary: #c3c2b7;
    --text-muted: #898781;
    --grid: #2c2c2a;
    --axis: #383835;
    --border: rgba(255,255,255,0.10);
    --series-oemmpa: #3987e5;
    --series-mmpdb: #199e70;
    --series-rdkit: #c98500;
    --good: #0ca30c;
    --warning: #fab219;
    --serious: #ec835a;
    --critical: #d03b3b;
    --shadow: 0 1px 2px rgba(0,0,0,0.4), 0 2px 8px rgba(0,0,0,0.3);
  }
}
* { box-sizing: border-box; }
body {
  margin: 0;
  background: var(--page);
  color: var(--text-primary);
  font-family: system-ui, -apple-system, "Segoe UI", sans-serif;
  line-height: 1.5;
  -webkit-font-smoothing: antialiased;
}
.wrap { max-width: 1080px; margin: 0 auto; padding: 32px 24px 80px; }
header.report { margin-bottom: 8px; }
h1 { font-size: 30px; font-weight: 650; margin: 0 0 6px; letter-spacing: -0.01em; }
h2 { font-size: 20px; font-weight: 620; margin: 40px 0 4px; letter-spacing: -0.01em; }
p.lead { color: var(--text-secondary); margin: 0 0 4px; max-width: 70ch; }
p.explain { color: var(--text-secondary); font-size: 14px; margin: 6px 0 16px; max-width: 74ch; }
.card {
  background: var(--surface-1);
  border: 1px solid var(--border);
  border-radius: 12px;
  padding: 18px 18px 12px;
  box-shadow: var(--shadow);
  margin-bottom: 8px;
}
.insight {
  background: var(--surface-1);
  border: 1px solid var(--border);
  border-left: 3px solid var(--series-oemmpa);
  border-radius: 10px;
  padding: 12px 16px;
  margin: 4px 0 16px;
  font-size: 13.5px;
  color: var(--text-secondary);
  max-width: 84ch;
  box-shadow: var(--shadow);
}
.insight b { color: var(--text-primary); }
.chart-grid { display: grid; grid-template-columns: repeat(auto-fit, minmax(340px, 1fr)); gap: 14px; }
.chart-title { font-size: 13px; font-weight: 600; color: var(--text-secondary); margin: 2px 2px 10px; }
.meta-grid { display: flex; flex-wrap: wrap; gap: 6px 22px; font-size: 13px;
  color: var(--text-secondary); margin: 10px 0 4px; }
.meta-grid b { color: var(--text-primary); font-weight: 600; }
.kpis { display: grid; grid-template-columns: repeat(auto-fit, minmax(200px, 1fr));
  gap: 14px; margin: 12px 0 8px; }
.kpi { background: var(--surface-1); border: 1px solid var(--border);
  border-radius: 12px; padding: 16px 18px; box-shadow: var(--shadow); }
.kpi .label { font-size: 13px; color: var(--text-secondary); margin-bottom: 6px; }
.kpi .value { font-size: 34px; font-weight: 660; letter-spacing: -0.02em; }
.kpi .sub { font-size: 12px; color: var(--text-muted); margin-top: 4px; }
.controls { display: flex; align-items: center; gap: 10px; margin: 4px 0 12px;
  font-size: 14px; color: var(--text-secondary); }
select { font: inherit; padding: 5px 8px; border-radius: 8px;
  border: 1px solid var(--axis); background: var(--surface-1);
  color: var(--text-primary); }
.legend { display: flex; flex-wrap: wrap; gap: 6px 16px; margin: 4px 0 2px;
  font-size: 13px; }
.legend .item { display: inline-flex; align-items: center; gap: 7px; cursor: pointer;
  user-select: none; color: var(--text-secondary); }
.legend .item.off { opacity: 0.38; }
.legend .swatch { width: 14px; height: 3px; border-radius: 2px; display: inline-block; }
.legend .swatch.box { width: 12px; height: 12px; border-radius: 3px; }
svg { display: block; width: 100%; height: auto; overflow: visible; }
svg text { fill: var(--text-muted); font-size: 12px; }
svg text.value { fill: var(--text-secondary); font-weight: 600; }
svg text.axis-title { fill: var(--text-secondary); font-size: 12px; }
.viz-tooltip {
  position: fixed; pointer-events: none; z-index: 20;
  background: var(--surface-1); color: var(--text-primary);
  border: 1px solid var(--border); border-radius: 8px;
  box-shadow: var(--shadow); padding: 8px 10px; font-size: 12.5px;
  min-width: 120px; opacity: 0; transition: opacity 0.08s; }
.viz-tooltip .tt-title { font-weight: 650; margin-bottom: 5px; }
.viz-tooltip .tt-row { display: flex; align-items: center; gap: 7px;
  justify-content: space-between; }
.viz-tooltip .tt-key { display: inline-flex; align-items: center; gap: 6px;
  color: var(--text-secondary); }
.viz-tooltip .tt-line { width: 12px; height: 3px; border-radius: 2px; }
.viz-tooltip .tt-val { font-weight: 650; font-variant-numeric: tabular-nums; }
table.data { border-collapse: collapse; width: 100%; font-size: 13px;
  font-variant-numeric: tabular-nums; }
table.data th, table.data td { text-align: right; padding: 6px 10px;
  border-bottom: 1px solid var(--grid); }
table.data th:first-child, table.data td:first-child,
table.data th.l, table.data td.l { text-align: left; font-variant-numeric: normal; }
table.data th { color: var(--text-secondary); font-weight: 600;
  position: sticky; top: 0; background: var(--surface-1); }
.table-scroll { max-height: 460px; overflow: auto; border: 1px solid var(--border);
  border-radius: 10px; }
.caveats { font-size: 13.5px; color: var(--text-secondary); }
.caveats li { margin-bottom: 7px; max-width: 78ch; }
.dot { width: 9px; height: 9px; border-radius: 50%; display: inline-block; }
.foot { color: var(--text-muted); font-size: 12px; margin-top: 40px; }
"""

_JS = r"""
const RAW = JSON.parse(document.getElementById("benchmark-data").textContent);
const RECORDS = RAW.records || [];
const META = RAW.meta || {};
const NS = "http://www.w3.org/2000/svg";

const TOOL_ORDER = ["oemmpa", "mmpdb", "rdkit"];
const TOOL_LABEL = { oemmpa: "OEMMPA", mmpdb: "mmpdb", rdkit: "RDKit" };
const TOOL_VAR = { oemmpa: "--series-oemmpa", mmpdb: "--series-mmpdb", rdkit: "--series-rdkit" };
const STAGE_ORDER = ["load", "fragment", "enumerate", "transforms", "materialize", "persist"];
const STAGE_LABEL = {
  load: "Load", fragment: "Fragment", enumerate: "Enumerate",
  transforms: "Transforms", materialize: "Materialize", persist: "Persist",
};
const STAGE_RAMP_LIGHT = ["#86b6ef", "#5598e7", "#2a78d6", "#1c5cab", "#184f95", "#104281"];
const STAGE_RAMP_DARK = ["#b7d3f6", "#86b6ef", "#5598e7", "#3987e5", "#256abf", "#184f95"];

const isDark = () => matchMedia("(prefers-color-scheme: dark)").matches;
const cssVar = (name) => getComputedStyle(document.body).getPropertyValue(name).trim();
const toolColor = (tool) => cssVar(TOOL_VAR[tool]) || "#888";
const stageRamp = () => (isDark() ? STAGE_RAMP_DARK : STAGE_RAMP_LIGHT);
const stageColor = (stage) => {
  const i = STAGE_ORDER.indexOf(stage);
  return stageRamp()[i < 0 ? 0 : i % stageRamp().length];
};

function el(tag, attrs, parent) {
  const e = document.createElementNS(NS, tag);
  if (attrs) for (const k in attrs) e.setAttribute(k, attrs[k]);
  if (parent) parent.appendChild(e);
  return e;
}
function text(node, value) { node.textContent = value == null ? "" : String(value); return node; }

function fmtTime(s) {
  if (s == null) return "-";
  if (s < 1) return (s * 1000).toFixed(s < 0.1 ? 1 : 0) + " ms";
  if (s < 60) return s.toFixed(2) + " s";
  return (s / 60).toFixed(1) + " min";
}
function fmtCompact(n) {
  if (n == null) return "-";
  const a = Math.abs(n);
  if (a >= 1e6) return (n / 1e6).toFixed(1) + "M";
  if (a >= 1e3) return (n / 1e3).toFixed(a >= 1e4 ? 0 : 1) + "K";
  return String(Math.round(n));
}
function fmtNum(n) { return n == null ? "-" : n.toLocaleString("en-US"); }

const scaling = RECORDS.filter((r) => r.benchmark === "stage_scaling");
const parallel = RECORDS.filter((r) => r.benchmark === "stage_parallel");
const uniq = (a) => Array.from(new Set(a));
const sizes = uniq(scaling.filter((r) => r.variant === "filtered").map((r) => r.size)).sort((a, b) => a - b);
const toolsPresent = TOOL_ORDER.filter((t) => scaling.some((r) => r.tool === t && r.variant === "filtered"));

function totalSeconds(tool, variant, size) {
  const rows = scaling.filter((r) => r.tool === tool && r.variant === variant && r.size === size);
  if (!rows.length) return null;
  return rows.reduce((s, r) => s + (r.seconds || 0), 0);
}
function stageSeconds(tool, variant, size, stage) {
  const r = scaling.find((x) => x.tool === tool && x.variant === variant && x.size === size && x.stage === stage);
  return r ? r.seconds : null;
}
function moleculeCount(size) {
  const r = scaling.find((x) => x.size === size);
  return r ? r.molecule_count : size;
}
function pairCount(tool, variant, size) {
  const r = scaling.find((x) => x.tool === tool && x.variant === variant && x.size === size);
  return r ? r.pair_count : null;
}

// ---- shared tooltip ----
const tooltip = document.createElement("div");
tooltip.className = "viz-tooltip";
document.body.appendChild(tooltip);
function showTip(html, evt) {
  tooltip.replaceChildren(...html);
  tooltip.style.opacity = "1";
  const pad = 14;
  let x = evt.clientX + pad, y = evt.clientY + pad;
  const w = tooltip.offsetWidth, h = tooltip.offsetHeight;
  if (x + w > innerWidth) x = evt.clientX - w - pad;
  if (y + h > innerHeight) y = evt.clientY - h - pad;
  tooltip.style.left = x + "px";
  tooltip.style.top = y + "px";
}
function hideTip() { tooltip.style.opacity = "0"; }
function tipTitle(s) { const d = document.createElement("div"); d.className = "tt-title"; d.textContent = s; return d; }
function tipRow(color, key, val) {
  const row = document.createElement("div"); row.className = "tt-row";
  const k = document.createElement("span"); k.className = "tt-key";
  const line = document.createElement("span"); line.className = "tt-line";
  line.style.background = color; k.appendChild(line);
  k.appendChild(document.createTextNode(key));
  const v = document.createElement("span"); v.className = "tt-val"; v.textContent = val;
  row.appendChild(k); row.appendChild(v); return row;
}

// ---- scales ----
function scaleLinear(d0, d1, r0, r1) {
  const span = d1 - d0 || 1;
  return (v) => r0 + ((v - d0) / span) * (r1 - r0);
}
function scaleLog(d0, d1, r0, r1) {
  const l0 = Math.log10(d0 <= 0 ? 1 : d0), l1 = Math.log10(d1 <= 0 ? 1 : d1);
  const span = l1 - l0 || 1;
  return (v) => r0 + ((Math.log10(v <= 0 ? d0 : v) - l0) / span) * (r1 - r0);
}
function niceLogTicks(d0, d1) {
  const ticks = [];
  let p = Math.floor(Math.log10(d0));
  const end = Math.ceil(Math.log10(d1));
  for (; p <= end; p++) ticks.push(Math.pow(10, p));
  return ticks.filter((t) => t >= d0 * 0.9 && t <= d1 * 1.1);
}

const M = { top: 18, right: 76, bottom: 44, left: 60 };
const H = 320;

// ---- generic multi-series line chart (log or linear axes) ----
function lineChart(mount, cfg) {
  mount.replaceChildren();
  const W = Math.max(mount.clientWidth || 640, 320);
  const iw = W - M.left - M.right, ih = H - M.top - M.bottom;
  const svg = el("svg", { viewBox: `0 0 ${W} ${H}`, width: W, height: H }, mount);
  const active = cfg.active || {};

  const allX = uniq(cfg.series.flatMap((s) => s.points.map((p) => p.x))).sort((a, b) => a - b);
  const visSeries = cfg.series.filter((s) => active[s.key] !== false);
  const ys = visSeries.flatMap((s) => s.points.map((p) => p.y)).filter((v) => v != null && v > 0);
  if (!allX.length || !ys.length) { text(el("text", { x: M.left, y: M.top + 20 }, svg), "No data"); return; }
  const xd0 = allX[0], xd1 = allX[allX.length - 1];
  const idealMax = cfg.ideal ? cfg.ideal(xd1) : -Infinity;
  const yd0 = Math.min(...ys), yd1 = Math.max(Math.max(...ys), idealMax);
  const sx = cfg.xLog ? scaleLog(xd0, xd1, M.left, M.left + iw) : scaleLinear(xd0, xd1, M.left, M.left + iw);
  const yPad = cfg.yLog ? 1 : (yd1 - yd0) * 0.08 || 1;
  const sy = cfg.yLog
    ? scaleLog(yd0, yd1, M.top + ih, M.top)
    : scaleLinear(Math.min(0, yd0), yd1 + yPad, M.top + ih, M.top);

  // gridlines + y ticks
  const yTicks = cfg.yLog ? niceLogTicks(yd0, yd1) : linTicks(cfg.yMin != null ? cfg.yMin : Math.min(0, yd0), yd1 + yPad, 5);
  yTicks.forEach((t) => {
    const y = sy(t);
    el("line", { x1: M.left, y1: y, x2: M.left + iw, y2: y, stroke: "var(--grid)", "stroke-width": 1 }, svg);
    text(el("text", { x: M.left - 8, y: y + 4, "text-anchor": "end" }, svg), cfg.yFormat ? cfg.yFormat(t) : t);
  });
  // x ticks
  const xTicks = cfg.xLog ? niceLogTicks(xd0, xd1) : allX;
  xTicks.forEach((t) => {
    const x = sx(t);
    text(el("text", { x, y: M.top + ih + 20, "text-anchor": "middle" }, svg), cfg.xFormat ? cfg.xFormat(t) : t);
  });
  // axis titles
  text(el("text", { x: M.left + iw / 2, y: H - 6, "text-anchor": "middle", class: "axis-title" }, svg), cfg.xLabel || "");
  const yt = el("text", { x: 14, y: M.top + ih / 2, "text-anchor": "middle", class: "axis-title",
    transform: `rotate(-90 14 ${M.top + ih / 2})` }, svg);
  text(yt, cfg.yLabel || "");

  // ideal reference line (parallel chart)
  if (cfg.ideal) {
    const pts = allX.map((x) => `${sx(x)},${sy(cfg.ideal(x))}`).join(" ");
    el("polyline", { points: pts, fill: "none", stroke: "var(--axis)", "stroke-width": 1.5,
      "stroke-dasharray": "4 4" }, svg);
    const lx = allX[allX.length - 1];
    text(el("text", { x: sx(lx) + 6, y: sy(cfg.ideal(lx)) + 4 }, svg), "ideal");
  }

  // series lines
  visSeries.forEach((s) => {
    const color = s.color();
    const pts = s.points.filter((p) => p.y != null && (!cfg.yLog || p.y > 0));
    if (pts.length > 1) {
      const d = pts.map((p) => `${sx(p.x)},${sy(p.y)}`).join(" ");
      el("polyline", { points: d, fill: "none", stroke: color, "stroke-width": 2,
        "stroke-linejoin": "round", "stroke-linecap": "round",
        "stroke-dasharray": s.dashed ? "6 4" : "none" }, svg);
    }
    pts.forEach((p) => {
      el("circle", { cx: sx(p.x), cy: sy(p.y), r: 4, fill: color,
        stroke: "var(--surface-1)", "stroke-width": 2 }, svg);
    });
    // direct end label; a colored end-dot carries identity, the text stays in ink
    if (pts.length && cfg.directLabels !== false) {
      const last = pts[pts.length - 1];
      text(el("text", { x: sx(last.x) + 9, y: sy(last.y) + 4, class: "value" }, svg), s.label);
    }
  });

  // crosshair + tooltip (snap to nearest x)
  const hair = el("line", { x1: 0, y1: M.top, x2: 0, y2: M.top + ih, stroke: "var(--axis)",
    "stroke-width": 1, opacity: 0 }, svg);
  const hit = el("rect", { x: M.left, y: M.top, width: iw, height: ih, fill: "transparent" }, svg);
  hit.addEventListener("pointermove", (evt) => {
    const rect = svg.getBoundingClientRect();
    const px = (evt.clientX - rect.left) * (W / rect.width);
    let nx = allX[0], best = Infinity;
    allX.forEach((x) => { const dd = Math.abs(sx(x) - px); if (dd < best) { best = dd; nx = x; } });
    hair.setAttribute("x1", sx(nx)); hair.setAttribute("x2", sx(nx)); hair.setAttribute("opacity", 1);
    const rows = [tipTitle(cfg.xTipLabel ? cfg.xTipLabel(nx) : String(nx))];
    visSeries.forEach((s) => {
      const p = s.points.find((q) => q.x === nx);
      if (p && p.y != null) rows.push(tipRow(s.color(), s.label, cfg.yTip ? cfg.yTip(p.y) : cfg.yFormat(p.y)));
    });
    showTip(rows, evt);
  });
  hit.addEventListener("pointerleave", () => { hair.setAttribute("opacity", 0); hideTip(); });
}

function linTicks(d0, d1, n) {
  const step = niceStep((d1 - d0) / n);
  const ticks = [];
  for (let t = Math.ceil(d0 / step) * step; t <= d1 + 1e-9; t += step) ticks.push(Number(t.toFixed(6)));
  return ticks;
}
function niceStep(raw) {
  const p = Math.pow(10, Math.floor(Math.log10(raw)));
  const f = raw / p;
  return (f <= 1 ? 1 : f <= 2 ? 2 : f <= 5 ? 5 : 10) * p;
}

// ---- stacked bar chart (stage breakdown per tool) ----
function stackedBars(mount, size) {
  mount.replaceChildren();
  const W = Math.max(mount.clientWidth || 640, 320);
  const iw = W - M.left - M.right, ih = H - M.top - M.bottom;
  const svg = el("svg", { viewBox: `0 0 ${W} ${H}`, width: W, height: H }, mount);
  const bars = toolsPresent.map((tool) => {
    const segs = STAGE_ORDER
      .map((stage) => ({ stage, value: stageSeconds(tool, "filtered", size, stage) }))
      .filter((s) => s.value != null);
    return { tool, segs, total: segs.reduce((a, b) => a + b.value, 0) };
  }).filter((b) => b.segs.length);
  if (!bars.length) { text(el("text", { x: M.left, y: M.top + 20 }, svg), "No data"); return; }
  const maxTotal = Math.max(...bars.map((b) => b.total));
  const sy = scaleLinear(0, maxTotal * 1.12, M.top + ih, M.top);
  const band = iw / bars.length;
  const bw = Math.min(64, band * 0.5);

  linTicks(0, maxTotal * 1.12, 5).forEach((t) => {
    const y = sy(t);
    el("line", { x1: M.left, y1: y, x2: M.left + iw, y2: y, stroke: "var(--grid)", "stroke-width": 1 }, svg);
    text(el("text", { x: M.left - 8, y: y + 4, "text-anchor": "end" }, svg), fmtTime(t));
  });
  text(el("text", { x: 14, y: M.top + ih / 2, "text-anchor": "middle", class: "axis-title",
    transform: `rotate(-90 14 ${M.top + ih / 2})` }, svg), "seconds (end-to-end)");

  bars.forEach((b, i) => {
    const cx = M.left + band * i + band / 2;
    let cursor = 0;
    b.segs.forEach((seg) => {
      const y0 = sy(cursor), y1 = sy(cursor + seg.value);
      const gap = (y0 - y1) > 4 ? 2 : 0;
      const color = stageColor(seg.stage);
      const rect = el("rect", { x: cx - bw / 2, y: y1, width: bw, height: Math.max(0.5, y0 - y1 - gap),
        rx: 3, fill: color }, svg);
      rect.addEventListener("pointermove", (evt) => showTip(
        [tipTitle(`${TOOL_LABEL[b.tool]} · ${STAGE_LABEL[seg.stage]}`),
          tipRow(color, "time", fmtTime(seg.value)),
          tipRow(color, "share", ((seg.value / b.total) * 100).toFixed(0) + "%")], evt));
      rect.addEventListener("pointerleave", hideTip);
      cursor += seg.value;
    });
    text(el("text", { x: cx, y: sy(b.total) - 8, "text-anchor": "middle", class: "value" }, svg), fmtTime(b.total));
    text(el("text", { x: cx, y: M.top + ih + 20, "text-anchor": "middle" }, svg), TOOL_LABEL[b.tool]);
  });
}

// ---- KPI tiles ----
function renderKpis() {
  const mount = document.getElementById("kpis");
  const top = sizes[sizes.length - 1];
  const tiles = [];
  const oe = totalSeconds("oemmpa", "filtered", top);
  if (oe != null) tiles.push({ label: `OEMMPA end-to-end @ ${fmtCompact(top)} molecules`,
    value: fmtTime(oe), sub: "load → fragment → enumerate → transforms → materialize → persist" });
  const mm = totalSeconds("mmpdb", "filtered", top);
  if (oe != null && mm != null) tiles.push({ label: `Faster than mmpdb @ ${fmtCompact(top)}`,
    value: (mm / oe).toFixed(1) + "×", sub: `${fmtTime(mm)} vs ${fmtTime(oe)} (fragment+index vs full pipeline)` });
  const rk = totalSeconds("rdkit", "filtered", top);
  if (oe != null && rk != null) tiles.push({ label: `Faster than RDKit @ ${fmtCompact(top)}`,
    value: (rk / oe).toFixed(1) + "×", sub: `equal-work (max-variable-heavies=10)` });
  const totals = parallel.filter((r) => r.stage === "total" && r.speedup != null);
  if (totals.length) {
    const best = totals.reduce((a, b) => (b.speedup > a.speedup ? b : a));
    tiles.push({ label: "Best OEMMPA parallel speedup",
      value: best.speedup.toFixed(1) + "×", sub: `${best.threads} threads, end-to-end @ ${fmtCompact(best.size)} molecules` });
  }
  mount.replaceChildren();
  tiles.forEach((t) => {
    const card = document.createElement("div"); card.className = "kpi";
    const l = document.createElement("div"); l.className = "label"; l.textContent = t.label;
    const v = document.createElement("div"); v.className = "value"; v.textContent = t.value;
    const s = document.createElement("div"); s.className = "sub"; s.textContent = t.sub;
    card.append(l, v, s); mount.appendChild(card);
  });
}

// ---- legends with toggle ----
function buildLegend(mount, items, state, onToggle) {
  mount.replaceChildren();
  items.forEach((it) => {
    const span = document.createElement("span");
    span.className = "item" + (state[it.key] === false ? " off" : "");
    const sw = document.createElement("span");
    sw.className = "swatch" + (it.box ? " box" : "");
    sw.style.background = it.color();
    span.append(sw, document.createTextNode(it.label));
    span.addEventListener("click", () => { state[it.key] = state[it.key] === false; onToggle(); });
    mount.appendChild(span);
  });
}

// ---- size scaling chart ----
const scalingActive = {};
function renderScaling() {
  const mount = document.getElementById("chart-scaling");
  const series = [];
  toolsPresent.forEach((tool) => {
    series.push({ key: tool, label: TOOL_LABEL[tool], color: () => toolColor(tool),
      points: sizes.map((s) => ({ x: s, y: totalSeconds(tool, "filtered", s) })).filter((p) => p.y != null) });
  });
  const unSizes = uniq(scaling.filter((r) => r.tool === "rdkit" && r.variant === "unfiltered").map((r) => r.size)).sort((a, b) => a - b);
  if (unSizes.length) {
    series.push({ key: "rdkit_unfiltered", label: "RDKit (native)", color: () => toolColor("rdkit"), dashed: true,
      points: unSizes.map((s) => ({ x: s, y: totalSeconds("rdkit", "unfiltered", s) })).filter((p) => p.y != null) });
  }
  buildLegend(document.getElementById("legend-scaling"),
    series.map((s) => ({ key: s.key, label: s.label, color: s.color })), scalingActive, renderScaling);
  lineChart(mount, {
    series, active: scalingActive, xLog: true, yLog: true,
    xLabel: "molecules (log scale)", yLabel: "seconds (log scale)",
    xFormat: fmtCompact, yFormat: fmtTime,
    xTipLabel: (x) => fmtCompact(x) + " molecules", yTip: fmtTime,
  });
}

// ---- parallel chart ----
const parallelStages = ["fragment", "enumerate", "materialize", "persist", "total"];
const parallelActive = {};
function renderParallel() {
  const mount = document.getElementById("chart-parallel");
  const psize = uniq(parallel.map((r) => r.size)).sort((a, b) => b - a)[0];
  const threads = uniq(parallel.map((r) => r.threads)).sort((a, b) => a - b);
  const ramp = stageRamp();
  const series = parallelStages.map((stage, i) => ({
    key: stage,
    label: stage === "total" ? "Total" : STAGE_LABEL[stage],
    color: () => (stage === "total" ? toolColor("oemmpa") : ramp[STAGE_ORDER.indexOf(stage) % ramp.length]),
    points: threads.map((t) => {
      const r = parallel.find((x) => x.threads === t && x.stage === stage && x.size === psize);
      return { x: t, y: r ? r.speedup : null };
    }).filter((p) => p.y != null),
  })).filter((s) => s.points.length);
  buildLegend(document.getElementById("legend-parallel"),
    series.map((s) => ({ key: s.key, label: s.label, color: s.color })), parallelActive, renderParallel);
  lineChart(mount, {
    series, active: parallelActive, xLog: false, yLog: false, yMin: 0,
    ideal: (x) => x, xLabel: "worker threads", yLabel: "speedup vs 1 thread",
    xFormat: (t) => t + "T", yFormat: (v) => v.toFixed(1) + "×",
    xTipLabel: (x) => x + " threads", yTip: (v) => v.toFixed(2) + "×",
  });
  const cap = document.getElementById("parallel-caption");
  if (cap) cap.textContent = `OEMMPA at ${fmtCompact(psize)} molecules. Fragment and enumerate scale with threads; materialize and persist stay flat (single-threaded by design).`;
}

// ---- throughput chart (molecules/sec, comparable across tools) ----
function renderThroughput() {
  const mount = document.getElementById("chart-throughput");
  const series = toolsPresent.map((tool) => ({
    key: tool, label: TOOL_LABEL[tool], color: () => toolColor(tool),
    points: sizes.map((s) => {
      const t = totalSeconds(tool, "filtered", s);
      return { x: s, y: t ? moleculeCount(s) / t : null };
    }).filter((p) => p.y != null),
  }));
  lineChart(mount, {
    series, xLog: true, yLog: false, yMin: 0,
    xLabel: "molecules (log scale)", yLabel: "molecules / second",
    xFormat: fmtCompact, yFormat: fmtCompact,
    xTipLabel: (x) => fmtCompact(x) + " molecules", yTip: (v) => fmtCompact(v) + " mol/s",
  });
}

// ---- throughput per unit of MMP work (pairs/sec, rises with scale) ----
function renderThroughputPairs() {
  const mount = document.getElementById("chart-throughput-pairs");
  const series = toolsPresent.map((tool) => ({
    key: tool, label: TOOL_LABEL[tool], color: () => toolColor(tool),
    points: sizes.map((s) => {
      const t = totalSeconds(tool, "filtered", s);
      const p = pairCount(tool, "filtered", s);
      return { x: s, y: t && p ? p / t : null };
    }).filter((p) => p.y != null),
  }));
  lineChart(mount, {
    series, xLog: true, yLog: false, yMin: 0,
    xLabel: "molecules (log scale)", yLabel: "pairs / second",
    xFormat: fmtCompact, yFormat: fmtCompact,
    xTipLabel: (x) => fmtCompact(x) + " molecules", yTip: (v) => fmtCompact(v) + " pairs/s",
  });
}

// ---- dynamic throughput note (computed from the loaded data, not hard-coded) ----
function metricTrend(tool, kind) {
  const pts = sizes.map((s) => {
    const t = totalSeconds(tool, "filtered", s);
    if (!t) return null;
    const num = kind === "mol" ? moleculeCount(s) : pairCount(tool, "filtered", s);
    return num ? num / t : null;
  }).filter((v) => v != null && v > 0);
  return pts.length >= 2 ? pts[pts.length - 1] / pts[0] : null;
}
function fmtTrend(ratio) {
  if (ratio == null) return "— (single size)";
  const arrow = ratio > 1.05 ? "▲" : ratio < 0.95 ? "▼" : "→";
  const factor = ratio >= 1 ? ratio.toFixed(1) : ratio.toFixed(2);
  return `${arrow} ${factor}× over the sweep`;
}
function renderThroughputNote() {
  const mount = document.getElementById("throughput-note");
  if (!mount) return;
  mount.replaceChildren();
  const lead = document.createElement("div");
  lead.style.marginBottom = "8px";
  lead.textContent =
    "Trends across this run's size range (per-molecule vs per-pair). Diverging " +
    "arrows for a tool — molecules/s down while pairs/s up — mean its cost " +
    "is pair-dominated and pairs are growing faster than molecules, not that it is " +
    "getting slower.";
  mount.appendChild(lead);
  toolsPresent.forEach((tool) => {
    const row = document.createElement("div");
    row.style.display = "flex";
    row.style.alignItems = "center";
    row.style.gap = "8px";
    row.style.margin = "3px 0";
    const dot = document.createElement("span");
    dot.className = "dot";
    dot.style.background = toolColor(tool);
    const label = document.createElement("span");
    label.textContent =
      `${TOOL_LABEL[tool]} — molecules/s ${fmtTrend(metricTrend(tool, "mol"))}` +
      `; pairs/s ${fmtTrend(metricTrend(tool, "pair"))}`;
    row.append(dot, label);
    mount.appendChild(row);
  });
}

// ---- stage-breakdown size selector ----
function initStageBreakdown() {
  const sel = document.getElementById("size-select");
  sizes.forEach((s) => {
    const opt = document.createElement("option");
    opt.value = String(s); opt.textContent = fmtCompact(s) + " molecules";
    sel.appendChild(opt);
  });
  sel.value = String(sizes[sizes.length - 1]);
  const draw = () => stackedBars(document.getElementById("chart-stages"), Number(sel.value));
  sel.addEventListener("change", draw);
  // stage legend (static, ordinal ramp)
  const leg = document.getElementById("legend-stages");
  STAGE_ORDER.forEach((stage) => {
    const span = document.createElement("span"); span.className = "item";
    const sw = document.createElement("span"); sw.className = "swatch box";
    sw.style.background = stageColor(stage);
    span.append(sw, document.createTextNode(STAGE_LABEL[stage]));
    leg.appendChild(span);
  });
  draw();
}

// ---- data table ----
function renderTable() {
  const body = document.getElementById("table-body");
  body.replaceChildren();
  scaling.slice().sort((a, b) =>
    a.size - b.size || TOOL_ORDER.indexOf(a.tool) - TOOL_ORDER.indexOf(b.tool) ||
    STAGE_ORDER.indexOf(a.stage) - STAGE_ORDER.indexOf(b.stage)
  ).forEach((r) => {
    const tr = document.createElement("tr");
    const cells = [
      [fmtNum(r.molecule_count), "l"], [TOOL_LABEL[r.tool], "l"], [r.variant, "l"],
      [STAGE_LABEL[r.stage] || r.stage, "l"], [fmtTime(r.seconds), ""],
      [fmtNum(r.pair_count), ""],
    ];
    cells.forEach(([v, cls]) => {
      const td = document.createElement("td"); if (cls) td.className = cls; td.textContent = v; tr.appendChild(td);
    });
    body.appendChild(tr);
  });
}

// ---- metadata line ----
function renderMeta() {
  const mount = document.getElementById("meta");
  const pairs = [
    ["Host", `${META.platform || "?"} · ${META.cpu_count || "?"} cores`],
    ["OEMMPA", META.oemmpa_version || "?"],
    ["RDKit", META.rdkit_version || "n/a"],
    ["mmpdb", META.mmpdb_available ? "available" : "n/a"],
    ["Filters", `max-variable-heavies=${(META.filters || {}).max_variable_heavies}, non-symmetric`],
    ["Generated", META.generated_at || "n/a"],
  ];
  mount.replaceChildren();
  pairs.forEach(([k, v]) => {
    const span = document.createElement("span");
    const b = document.createElement("b"); b.textContent = k + ": ";
    span.append(b, document.createTextNode(v));
    mount.appendChild(span);
  });
}

function renderAll() {
  renderMeta();
  renderKpis();
  initStageBreakdownOnce();
  renderScaling();
  renderParallel();
  renderThroughput();
  renderThroughputPairs();
  renderThroughputNote();
  renderTable();
}
let stageInit = false;
function initStageBreakdownOnce() {
  if (!stageInit) { initStageBreakdown(); stageInit = true; }
  else stackedBars(document.getElementById("chart-stages"), Number(document.getElementById("size-select").value));
}

renderAll();
let resizeTimer = null;
addEventListener("resize", () => { clearTimeout(resizeTimer); resizeTimer = setTimeout(renderAll, 150); });
matchMedia("(prefers-color-scheme: dark)").addEventListener("change", renderAll);
"""


def _page_body():
    """Return the static body HTML (charts are hydrated by the embedded JS)."""
    return """
  <div class="wrap">
    <header class="report">
      <h1>OEMMPA performance &mdash; staged benchmark</h1>
      <p class="lead">Per-stage timing of OEMMPA against RDKit and mmpdb across a
        molecule-count sweep and a thread sweep, on a public SureChEMBL corpus.</p>
      <div class="meta-grid" id="meta"></div>
    </header>

    <p class="explain">Every timing is the <b>minimum over repeated runs</b> (warm
      timing: the fastest observed run, which most cleanly reflects compute cost by
      excluding one-off scheduling noise). For a like-for-like comparison all three
      tools cap the variable fragment at <b>10 heavy atoms</b> &mdash; mmpdb's index
      default &mdash; so the pair surface is equal work. Absolute pair <i>counts</i>
      still differ by tool definition (see caveats); read the <i>speed</i>, not the
      raw counts.</p>

    <div class="kpis" id="kpis"></div>

    <h2>Stage breakdown</h2>
    <p class="explain">Where the end-to-end time goes, per tool, as a stack of
      pipeline stages. <b>Load</b> parses and desalts input; <b>fragment</b> cuts
      bonds; <b>enumerate</b> forms matched pairs (mmpdb's <code>index</code> also
      builds rules and writes SQLite here); <b>transforms</b> groups pairs into
      rules; <b>materialize</b> exports rows; <b>persist</b> writes DuckDB. mmpdb
      and RDKit expose only the shared core (fragment + enumerate).
      <i>How to read it &mdash;</i> a taller bar is more end-to-end time; each
      segment is one stage, so compare where each tool spends its time. Change the
      corpus size to see how the mix shifts as the dataset grows.</p>
    <div class="controls">
      <label for="size-select">Corpus size</label>
      <select id="size-select"></select>
    </div>
    <div class="card"><div id="chart-stages"></div></div>
    <div class="legend" id="legend-stages"></div>

    <h2>Scaling with corpus size</h2>
    <p class="explain">End-to-end wall time vs molecule count, both axes log-scaled.
      <i>How to read it &mdash;</i> lower is faster; on log-log axes a straight line
      means power-law growth and a shallower slope means better scaling, so compare
      slopes between tools. The dashed line is native (unfiltered) RDKit, shown only
      at small sizes because its pair count explodes past a few thousand molecules.
      Click a legend entry to toggle a series.</p>
    <div class="legend" id="legend-scaling"></div>
    <div class="card"><div id="chart-scaling"></div></div>

    <h2>Parallel scaling</h2>
    <p class="explain">OEMMPA speedup as worker threads increase, relative to a
      single thread. The dashed <i>ideal</i> line is perfect linear scaling; the gap
      to it is Amdahl's law &mdash; the serial fraction (materialize, persist) caps
      the achievable speedup. Click a legend entry to toggle a stage.</p>
    <div class="legend" id="legend-parallel"></div>
    <div class="card"><div id="chart-parallel"></div></div>
    <p class="explain" id="parallel-caption"></p>

    <h2>Throughput</h2>
    <p class="explain">The same runs, normalized two ways: by <b>input molecules</b>
      (left) and by <b>matched pairs produced</b> (right).
      <i>How to read it &mdash;</i> follow the <b>shape</b> of each tool's own curve,
      not the absolute heights between tools (pair counts are defined differently per
      tool, so per-pair heights are not comparable across tools). A per-molecule rate
      that falls while the per-pair rate rises is the signature of pair-dominated
      work &mdash; matched pairs grow faster than molecules as the corpus grows &mdash;
      rather than a slowdown. The note below is computed live from the loaded data.</p>
    <div class="insight" id="throughput-note"></div>
    <div class="chart-grid">
      <div class="card">
        <div class="chart-title">Molecules / second &mdash; normalized by input size</div>
        <div id="chart-throughput"></div>
      </div>
      <div class="card">
        <div class="chart-title">Pairs / second &mdash; normalized by MMP work
          (per-tool trend; pair definitions differ across tools)</div>
        <div id="chart-throughput-pairs"></div>
      </div>
    </div>

    <h2>All measurements</h2>
    <p class="explain">Every stage timing behind the charts. The table is the
      accessible fallback for the light-mode categorical colors (relief rule) and
      the source of truth for exact numbers.</p>
    <div class="table-scroll">
      <table class="data">
        <thead><tr>
          <th class="l">Molecules</th><th class="l">Tool</th><th class="l">Variant</th>
          <th class="l">Stage</th><th>Time</th><th>Pairs</th>
        </tr></thead>
        <tbody id="table-body"></tbody>
      </table>
    </div>

    <h2>Method &amp; caveats</h2>
    <ul class="caveats">
      <li><b>Pair counts are not comparable across tools.</b> mmpdb stores one row
        per (pair, environment radius), inflating its count several-fold; OEMMPA
        counts non-symmetric pairs once; RDKit counts differ again. Only the
        <i>timings</i> are compared like-for-like.</li>
      <li><b>mmpdb <code>index</code> bundles stages.</b> Its single index step
        covers enumeration, rule building, and SQLite persistence &mdash; so it is
        charted under <i>enumerate</i> and legitimately does more than OEMMPA's
        enumerate alone.</li>
      <li><b>mmpdb data stops where its RDKit backend aborts.</b> On raw
        SureChEMBL patent data, mmpdb's indexer (via RDKit 2026.03.3) hits a
        fatal canonical-ranking assertion on certain extended-polyene structures
        (retinoid / secosteroid analogs), so mmpdb produces no result at the
        larger sizes here. OEMMPA (OpenEye backend) and RDKit's own rdMMPA process
        the full corpus &mdash; a robustness difference, not a harness limitation.</li>
      <li><b>Native RDKit has no variable-size gate.</b> Unfiltered, its pair count
        grows explosively (tens of thousands of pairs at only a few hundred
        molecules), so it is captured only at small sizes for reference; the primary
        comparison filters it to match mmpdb/OEMMPA.</li>
      <li><b>Warm timing.</b> Each point is the fastest of repeated runs; repeats are
        reduced automatically at the largest sizes to bound total runtime.</li>
      <li><b>Corpus.</b> Deterministically sampled from the public SureChEMBL
        parquet (150&ndash;450 Da, connected), nested across sizes for
        comparability.</li>
    </ul>

    <p class="foot">Self-contained report &mdash; no network access required. Charts
      rendered as inline SVG; palette validated for color-vision deficiency.</p>
  </div>
"""


def render_html(records, meta):
    """Render the staged-benchmark records and metadata as one HTML document.

    :param records: List of ``stage_scaling`` / ``stage_parallel`` record dicts.
    :param meta: Run-metadata dict (host, versions, filters, timestamp).
    :returns: A complete, self-contained HTML document string.
    """
    payload = json.dumps({"records": records, "meta": meta})
    # Neutralize any "</script>"-like sequence inside the embedded JSON.
    payload = payload.replace("</", "<\\/")
    return (
        "<!doctype html>\n"
        '<html lang="en">\n<head>\n'
        '<meta charset="utf-8">\n'
        '<meta name="viewport" content="width=device-width, initial-scale=1">\n'
        "<title>OEMMPA staged benchmark</title>\n"
        "<style>" + _CSS + "</style>\n"
        "</head>\n"
        '<body data-palette="#2a78d6,#1baf7a,#eda100">\n'
        + _page_body()
        + '\n<script type="application/json" id="benchmark-data">'
        + payload
        + "</script>\n"
        + "<script>" + _JS + "</script>\n"
        + "</body>\n</html>\n"
    )
