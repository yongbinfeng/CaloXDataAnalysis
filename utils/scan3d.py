"""Interactive 3D stacked heatmaps: the same per-run 2D calorimeter-face maps
as make_adc_position_maps.py / make_toa_position_maps.py, but instead of one
page per scan position, every run's flat map is drawn as a solid tile plane
at its table position along a length axis, all stacked and visible at once
(no slider, no scrolling - just the same heatmaps stacked in X).

Each channel becomes one flat rectangular tile (two triangles, flat-shaded so
the colour doesn't interpolate across the tile) positioned at its physical
(x, y) tower location and at z = that run's table position. All runs of one
channel family are merged into a single mesh3d trace; a button group toggles
which family (CerQuartz / CerPlastic / Sci) is visible, since their ADC/TOA
ranges differ by more than an order of magnitude.

Self-contained HTML page, Plotly.js from CDN, no plotly Python package
required.
"""

from __future__ import annotations

import json
import math
from pathlib import Path
from string import Template
from typing import Any, Dict, List, Sequence, Tuple

import numpy as np

from channels.channel_map import build_drs_boards
from utils.visualization import _X_SCALE, _Y_SCALE

CATEGORY_NAMES = ("CerQuartz", "CerPlastic", "Sci")

_PLOTLY_CDN = "https://cdn.plot.ly/plotly-2.35.2.min.js"


def build_channel_geometry(map_run: int = 1897) -> Dict[str, List[Dict[str, Any]]]:
    """Per-category list of {branch, x, y, is6mm} in the tb2026 channel map."""
    geometry: Dict[str, List[Dict[str, Any]]] = {name: [] for name in CATEGORY_NAMES}
    drsboards = build_drs_boards(run_number=map_run)
    for _, board in drsboards.items():
        for channel in board:
            if channel.is_reference:
                continue
            if not channel.isCer:
                category = "Sci"
            elif channel.isQuartz:
                category = "CerQuartz"
            else:
                category = "CerPlastic"
            geometry[category].append({
                "branch": channel.get_channel_name(blsub=False),
                "x": channel.i_tower_x * _X_SCALE,
                "y": channel.i_tower_y * _Y_SCALE,
                "is6mm": bool(channel.is6mm),
            })
    return geometry


def percentile_range(
    values: Sequence[float], percentile: float = 2.0,
) -> Tuple[float, float]:
    """Same trimmed range as make_{adc,toa}_position_maps.py's z_range_for_category."""
    pool = sorted(v for v in values if math.isfinite(v))
    if not pool:
        return 0.0, 1.0
    frac = max(0.0, min(0.45, percentile / 100.0))
    lo = pool[int(frac * (len(pool) - 1))]
    hi = pool[int((1.0 - frac) * (len(pool) - 1))]
    if hi <= lo:
        lo, hi = pool[0], max(pool[-1], pool[0] + 1.0)
    return lo, hi


def compute_toa_references(
    data: Dict[str, Any], branches: List[str], min_adc: float, min_ref_hits: int,
) -> Dict[str, float]:
    """Per-branch median _Time_cfd_mcp over this run's clean, real hits.

    _Time_cfd_mcp's absolute scale is dominated by per-channel cable-length
    offsets, not real timing. This is the single-run analogue of the
    --relative flag in make_toa_position_maps.py (which averages over
    repeated scan runs at the same position): each channel's own typical
    arrival time, estimated from every other event in this run where it had
    a real pulse, so a single event's deviation from that is what's left.
    Channels with fewer than min_ref_hits such events are dropped (no
    stable baseline to subtract).

    `data` is the dict returned by RDataFrame.AsNumpy for the whole run,
    keyed by `"{branch}_energy"` / `"{branch}_Time_cfd_mcp"`.
    """
    references: Dict[str, float] = {}
    for branch in branches:
        energy_arr = np.asarray(data["{}_energy".format(branch)], dtype=float)
        toa_arr = np.asarray(data["{}_Time_cfd_mcp".format(branch)], dtype=float)
        mask = np.isfinite(energy_arr) & np.isfinite(toa_arr) & (energy_arr >= min_adc)
        if int(mask.sum()) >= min_ref_hits:
            references[branch] = float(np.median(toa_arr[mask]))
    return references


def _tile_half_widths(is6mm: bool) -> Tuple[float, float]:
    """Half-width/height of one channel's tile, in cm.

    Matches the bin sizes visualizeDRSBoards uses for the two DRS histogram
    grids: the x pitch is the same for both, the 3 mm channels get a y grid
    4x finer than the 6 mm one.
    """
    hx = _X_SCALE / 2.0
    hy = (_Y_SCALE / 4.0 if not is6mm else _Y_SCALE) / 2.0
    return hx, hy


def _category_mesh(
    channels: List[Dict[str, Any]],
    runs: List[int],
    positions: Dict[int, float],
    values: Dict[int, Dict[str, float]],
    value_label: str,
    value_fmt: str,
) -> Dict[str, Any]:
    vx: List[float] = []
    vy: List[float] = []
    vz: List[float] = []
    fi: List[int] = []
    fj: List[int] = []
    fk: List[int] = []
    intensity: List[float] = []
    text: List[str] = []
    all_values: List[float] = []

    for run in runs:
        per_run = values.get(run, {})
        pos = positions[run]
        for channel in channels:
            branch = channel["branch"]
            val = per_run.get(branch)
            if val is None or not math.isfinite(val):
                continue
            hx, hy = _tile_half_widths(channel["is6mm"])
            cx, cy = channel["x"], channel["y"]
            size_label = "6mm" if channel["is6mm"] else "3mm"
            hover = "run {}<br>{} ({})<br>X = {:.1f} mm<br>{} = {}".format(
                run, branch, size_label, pos, value_label, value_fmt.format(val),
            )
            base = len(vx)
            for corner_x, corner_y in (
                (cx - hx, cy - hy), (cx + hx, cy - hy),
                (cx + hx, cy + hy), (cx - hx, cy + hy),
            ):
                vx.append(corner_x)
                vy.append(corner_y)
                vz.append(pos)
                intensity.append(val)
                text.append(hover)
            fi.extend((base, base))
            fj.extend((base + 1, base + 2))
            fk.extend((base + 2, base + 3))
            all_values.append(val)

    lo, hi = percentile_range(all_values)
    return {
        "x": vx, "y": vy, "z": vz,
        "i": fi, "j": fj, "k": fk,
        "intensity": intensity, "text": text,
        "range": [lo, hi],
    }


def write_scan3d_html(
    output_path: Path,
    page_title: str,
    subtitle: str,
    geometry: Dict[str, List[Dict[str, Any]]],
    runs: List[int],
    positions: Dict[int, float],
    values_by_category: Dict[str, Dict[int, Dict[str, float]]],
    value_label: str,
    value_fmt: str = "{:.1f}",
    z_axis_title: str = "table scan position X [mm]",
    colorscale: str = "Turbo",
) -> None:
    """Write one interactive HTML page covering every category in `geometry`.

    values_by_category[category][run][branch] = observable value for that
    run/channel. runs/positions give the scan axis (the "length" dimension);
    each run's map is one flat tile plane at z = positions[run], and every
    plane is drawn at once - stacked heatmaps, not a slice-by-slice slider.
    """
    ordered_runs = sorted(runs, key=lambda r: (positions[r], r))
    categories = [name for name in CATEGORY_NAMES if geometry.get(name)]

    meshes = {
        name: _category_mesh(
            geometry[name], ordered_runs, positions,
            values_by_category.get(name, {}), value_label, value_fmt,
        )
        for name in categories
    }

    payload = {
        "categories": categories,
        "data": meshes,
        "value_label": value_label,
        "z_axis_title": z_axis_title,
        "colorscale": colorscale,
        "page_title": page_title,
        "subtitle": subtitle,
    }

    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    output_path.write_text(_HTML_TEMPLATE.substitute(
        cdn=_PLOTLY_CDN,
        title=page_title,
        payload_json=json.dumps(payload),
    ))


_HTML_TEMPLATE = Template("""<!DOCTYPE html>
<html>
<head>
<meta charset="utf-8">
<title>$title</title>
<script src="$cdn"></script>
<style>
  body { font-family: -apple-system, Helvetica, Arial, sans-serif; margin: 0;
         background: #111; color: #eee; }
  #header { padding: 10px 16px; }
  #header h1 { font-size: 18px; margin: 0 0 4px 0; }
  #header p { font-size: 13px; margin: 0; color: #aaa; }
  #plot { width: 100vw; height: calc(100vh - 60px); }
</style>
</head>
<body>
<div id="header">
  <h1>$title</h1>
  <p id="subtitle"></p>
</div>
<div id="plot"></div>
<script>
const PAYLOAD = $payload_json;

function makeTrace(cat, active) {
  const d = PAYLOAD.data[cat];
  return {
    type: "mesh3d", name: cat,
    x: d.x, y: d.y, z: d.z,
    i: d.i, j: d.j, k: d.k,
    intensity: d.intensity, text: d.text, hoverinfo: "text",
    flatshading: true,
    colorscale: PAYLOAD.colorscale, cmin: d.range[0], cmax: d.range[1],
    showscale: true,
    colorbar: { title: PAYLOAD.value_label,
      titlefont: { color: "#eee" }, tickfont: { color: "#eee" } },
    visible: active,
    lighting: { ambient: 1, diffuse: 0, specular: 0 },
  };
}

const traces = PAYLOAD.categories.map((cat, i) => makeTrace(cat, i === 0));

const catButtons = PAYLOAD.categories.map((cat, i) => {
  const visible = PAYLOAD.categories.map((_, j) => j === i);
  return { label: cat, method: "restyle", args: [{ visible: visible }] };
});

const layout = {
  paper_bgcolor: "#111", plot_bgcolor: "#111",
  font: { color: "#eee" },
  margin: { l: 0, r: 0, t: 10, b: 0 },
  scene: {
    xaxis: { title: "channel X [cm]", color: "#ccc", gridcolor: "#333" },
    yaxis: { title: "channel Y [cm]", color: "#ccc", gridcolor: "#333" },
    zaxis: { title: PAYLOAD.z_axis_title, color: "#ccc", gridcolor: "#333" },
    aspectmode: "manual", aspectratio: { x: 1.2, y: 1.0, z: 2.2 },
  },
  updatemenus: [
    {
      type: "buttons", direction: "left", x: 0.02, y: 1.05, xanchor: "left",
      buttons: catButtons,
      bgcolor: "#222", bordercolor: "#555", font: { color: "#eee" },
    },
  ],
};

document.getElementById("subtitle").textContent = PAYLOAD.subtitle;

Plotly.newPlot("plot", traces, layout, { responsive: true });
</script>
</body>
</html>
""")


def _tile_marker_px(is6mm: bool) -> float:
    """Rough on-screen square size for one channel's tile, in pixels.

    Not data-space accurate (a 2D scatter marker's size is fixed in pixels,
    so it drifts from the true tile footprint on zoom) but keeps the 6 mm /
    3 mm size difference visible without the complexity of drawing real
    per-tile rectangles.
    """
    return 16.0 if is6mm else 8.0


def _category_gallery(
    channels: List[Dict[str, Any]],
    events: List[int],
    values: Dict[int, Dict[str, Tuple[float, float]]],
    value_label: str,
    value_fmt: str,
) -> Dict[str, Any]:
    """values[event][branch] = (colour_value, adc_value) for hover context."""
    by_branch = {c["branch"]: c for c in channels}
    frames: Dict[int, Dict[str, Any]] = {}
    all_values: List[float] = []

    for event in events:
        per_event = values.get(event, {})
        fx: List[float] = []
        fy: List[float] = []
        fc: List[float] = []
        ft: List[str] = []
        fs: List[float] = []
        for branch, (color_value, adc_value) in per_event.items():
            geo = by_branch.get(branch)
            if geo is None or not math.isfinite(color_value):
                continue
            size_label = "6mm" if geo["is6mm"] else "3mm"
            fx.append(geo["x"])
            fy.append(geo["y"])
            fc.append(color_value)
            fs.append(_tile_marker_px(geo["is6mm"]))
            ft.append(
                "evt {}<br>{} ({})<br>{} = {}<br>ADC = {:.1f}".format(
                    event, branch, size_label, value_label,
                    value_fmt.format(color_value), adc_value,
                )
            )
            all_values.append(color_value)
        frames[event] = {"x": fx, "y": fy, "c": fc, "t": ft, "s": fs}

    lo, hi = percentile_range(all_values)
    return {"frames": frames, "range": [lo, hi]}


def write_event_heatmap_gallery_html(
    output_path: Path,
    page_title: str,
    subtitle: str,
    geometry: Dict[str, List[Dict[str, Any]]],
    events: List[int],
    values_by_category: Dict[str, Dict[int, Dict[str, Tuple[float, float]]]],
    value_label: str,
    value_fmt: str = "{:.2f}",
    colorscale: str = "Turbo",
) -> None:
    """Write a 2D calorimeter-face heatmap with a slider that pages through
    events, one flat frame at a time (no stacking - each event replaces the
    last). values_by_category[category][event][branch] = (colour, adc).
    """
    categories = [name for name in CATEGORY_NAMES if geometry.get(name)]
    galleries = {
        name: _category_gallery(
            geometry[name], events, values_by_category.get(name, {}),
            value_label, value_fmt,
        )
        for name in categories
    }

    payload = {
        "categories": categories,
        "events": events,
        "data": galleries,
        "value_label": value_label,
        "colorscale": colorscale,
        "page_title": page_title,
        "subtitle": subtitle,
    }

    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    output_path.write_text(_GALLERY_HTML_TEMPLATE.substitute(
        cdn=_PLOTLY_CDN,
        title=page_title,
        payload_json=json.dumps(payload),
    ))


_GALLERY_HTML_TEMPLATE = Template("""<!DOCTYPE html>
<html>
<head>
<meta charset="utf-8">
<title>$title</title>
<script src="$cdn"></script>
<style>
  body { font-family: -apple-system, Helvetica, Arial, sans-serif; margin: 0;
         background: #111; color: #eee; }
  #header { padding: 10px 16px; }
  #header h1 { font-size: 18px; margin: 0 0 4px 0; }
  #header p { font-size: 13px; margin: 0; color: #aaa; }
  #plot { width: 100vw; height: calc(100vh - 60px); }
</style>
</head>
<body>
<div id="header">
  <h1>$title</h1>
  <p id="subtitle"></p>
</div>
<div id="plot"></div>
<script>
const PAYLOAD = $payload_json;

function makeMarker(color, size, range) {
  return {
    symbol: "square", size: size, color: color,
    colorscale: PAYLOAD.colorscale, cmin: range[0], cmax: range[1],
    showscale: true, line: { width: 0 },
    colorbar: { title: PAYLOAD.value_label,
      titlefont: { color: "#eee" }, tickfont: { color: "#eee" } },
  };
}

function makeTrace(cat, active) {
  const d = PAYLOAD.data[cat];
  const first_event = PAYLOAD.events[0];
  const fr = d.frames[first_event];
  return {
    type: "scattergl", mode: "markers", name: cat,
    x: fr.x, y: fr.y, text: fr.t, hoverinfo: "text",
    visible: active,
    marker: makeMarker(fr.c, fr.s, d.range),
  };
}

const traces = PAYLOAD.categories.map((cat, i) => makeTrace(cat, i === 0));

const frames = PAYLOAD.events.map((ev) => {
  const data = PAYLOAD.categories.map((cat) => {
    const fr = PAYLOAD.data[cat].frames[ev];
    const range = PAYLOAD.data[cat].range;
    return {
      x: fr.x, y: fr.y, text: fr.t,
      marker: makeMarker(fr.c, fr.s, range),
    };
  });
  const trace_indices = PAYLOAD.categories.map((_, i) => i);
  return { name: String(ev), data: data, traces: trace_indices };
});

const sliderSteps = PAYLOAD.events.map((ev) => ({
  method: "animate",
  label: "evt " + ev,
  args: [[String(ev)], { mode: "immediate",
    frame: { duration: 0, redraw: true }, transition: { duration: 0 } }],
}));

const catButtons = PAYLOAD.categories.map((cat, i) => {
  const visible = PAYLOAD.categories.map((_, j) => j === i);
  return { label: cat, method: "restyle", args: [{ visible: visible }] };
});

const layout = {
  paper_bgcolor: "#111", plot_bgcolor: "#111",
  font: { color: "#eee" },
  margin: { l: 60, r: 20, t: 10, b: 60 },
  xaxis: { title: "channel X [cm]", color: "#ccc", gridcolor: "#333",
    zeroline: false },
  yaxis: { title: "channel Y [cm]", color: "#ccc", gridcolor: "#333",
    zeroline: false, scaleanchor: "x", scaleratio: 1 },
  updatemenus: [
    {
      type: "buttons", direction: "left", x: 0.02, y: 1.08, xanchor: "left",
      buttons: catButtons,
      bgcolor: "#222", bordercolor: "#555", font: { color: "#eee" },
    },
    {
      type: "buttons", direction: "left", x: 0.5, y: 1.08, xanchor: "left",
      buttons: [
        { label: "Play", method: "animate",
          args: [null, { fromcurrent: true,
            frame: { duration: 400, redraw: true },
            transition: { duration: 0 } }] },
        { label: "Pause", method: "animate",
          args: [[null], { mode: "immediate",
            frame: { duration: 0, redraw: false } }] },
      ],
      bgcolor: "#222", bordercolor: "#555", font: { color: "#eee" },
    },
  ],
  sliders: [{
    active: 0, x: 0.05, y: 0, len: 0.9,
    currentvalue: { prefix: "event: ", font: { color: "#eee" } },
    pad: { t: 40 },
    font: { color: "#eee" },
    steps: sliderSteps,
  }],
};

document.getElementById("subtitle").textContent = PAYLOAD.subtitle;

Plotly.newPlot("plot", traces, layout, { responsive: true }).then(() => {
  Plotly.addFrames("plot", frames);
});
</script>
</body>
</html>
""")
