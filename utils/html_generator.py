import os
import re
from datetime import datetime
from zoneinfo import ZoneInfo

JSROOT_CDN = "https://root.cern/js/latest/modules/main.mjs"


def generate_html(png_files, png_dir, plots_per_row=4, output_html="view_plots.html", title="", intro_text=""):
    html_dir = os.path.dirname(os.path.abspath(output_html))
    if not title:
        title = os.path.splitext(os.path.basename(output_html))[0].replace('_', ' ')
    real_png_files = [f for f in png_files if f != "NEWLINE"]
    png_paths = [os.path.join(png_dir, f) for f in real_png_files]
    rel_paths = [os.path.relpath(p, start=html_dir) for p in png_paths]
    path_map = dict(zip(real_png_files, rel_paths))

    timestamp = datetime.now(ZoneInfo("Europe/Zurich")
                             ).strftime("%B %d, %Y, %I:%M %p %Z")

    # --- Intro Text Processing (Bold and Bullets) ---
    # 1. Bold text: **word** -> <b>word</b>
    processed_intro = re.sub(r'\*\*(.*?)\*\*', r'<b>\1</b>', intro_text)

    # 2. Bullet points: Lines starting with * or -
    lines = processed_intro.split('\n')
    formatted_lines = []
    in_list = False

    for line in lines:
        stripped = line.strip()
        if stripped.startswith(('* ', '- ')):
            if not in_list:
                formatted_lines.append(
                    '<ul style="margin: 5px 0; padding-left: 20px;">')
                in_list = True
            # Remove the marker (* or -) and wrap in <li>
            content = stripped[2:].strip()
            formatted_lines.append(f'<li>{content}</li>')
        else:
            if in_list:
                formatted_lines.append('</ul>')
                in_list = False
            if stripped.startswith('<'):
                formatted_lines.append(stripped)  # raw HTML block, no <br>
            else:
                formatted_lines.append(line + "<br>")

    if in_list:
        formatted_lines.append('</ul>')

    formatted_intro = "".join(formatted_lines)

    html_header = f"""<!DOCTYPE html>
<html>
<head>
  <meta charset="UTF-8">
  <title>{title}</title>
  <style>
    body {{ font-family: -apple-system, BlinkMacSystemFont, "Segoe UI", Roboto, sans-serif; padding: 10px 20px; margin: 0; line-height: 1.4; color: #333; }}
    .top-bar {{ display: flex; align-items: center; gap: 30px; padding: 10px 0; border-bottom: 1px solid #eee; margin-bottom: 15px; position: sticky; top: 0; background: white; z-index: 100; }}
    .title-group {{ display: flex; flex-direction: column; min-width: fit-content; }}
    h1 {{ font-size: 20px; margin: 0; color: #1a1a1a; white-space: nowrap; }}
    .timestamp {{ font-size: 11px; color: #888; white-space: nowrap; }}
    .controls {{ display: flex; align-items: center; gap: 12px; }}
    input[type="text"] {{ font-size: 14px; padding: 4px 10px; width: 250px; border: 1px solid #ccc; border-radius: 4px; }}
    .intro-text {{ background: #f8f9fa; border-left: 4px solid #007bff; padding: 10px 16px; margin-bottom: 15px; font-size: 14px; border-radius: 0 4px 4px 0; }}
    .intro-text li {{ margin-bottom: 4px; }}
    .grid {{ display: grid; grid-template-columns: repeat({plots_per_row}, 1fr); gap: 15px; }}
    .plot {{ border: 1px solid #e1e4e8; padding: 8px; text-align: center; border-radius: 6px; background: #fff; }}
    .filename {{ font-weight: 600; font-size: 12px; margin-bottom: 8px; word-break: break-all; color: #444; }}
    img {{ max-width: 100%; height: auto; border-radius: 2px; }}
    button {{ padding: 4px 12px; cursor: pointer; font-size: 13px; }}
    @media (max-width: 1100px) {{ .grid {{ grid-template-columns: repeat(2, 1fr); }} }}
    @media (max-width: 800px) {{ .top-bar {{ flex-direction: column; align-items: flex-start; gap: 10px; }} }}
  </style>
</head>
<body>
  <div class="top-bar">
    <div class="title-group">
      <h1>{title}</h1>
      <div class="timestamp">{timestamp}</div>
    </div>
    <div class="controls">
      <input type="text" id="filterInput" placeholder="Filter by filename..." onkeyup="filterPlots()">
      <label><input type="checkbox" id="regexToggle" onchange="filterPlots()"> Regex</label>
      <label><input type="checkbox" id="caseToggle" onchange="filterPlots()"> Case</label>
      <button onclick="clearFilter()">Clear</button>
    </div>
  </div>
  {f'<div class="intro-text">{formatted_intro}</div>' if intro_text else ''}
  <div id="plotContainer">
"""

    html_body = ""
    current_row_plots = []

    def create_grid_html(plots):
        if not plots:
            return ""
        grid_html = '    <div class="grid">\n'
        for filename in plots:
            rel_path = path_map[filename]
            grid_html += f"""      <div class="plot" data-filename="{filename}">
        <div class="filename">{filename}</div>
        <a href="{rel_path}" target="_blank">
          <img src="{rel_path}" alt="{filename}">
        </a>
      </div>\n"""
        grid_html += '    </div>\n'
        return grid_html

    for filename in png_files:
        if filename == "NEWLINE":
            html_body += create_grid_html(current_row_plots)
            current_row_plots = []
        else:
            current_row_plots.append(filename)
    html_body += create_grid_html(current_row_plots)

    html_footer = """  </div>
  <script>
    function getQueryParams() {
      const params = new URLSearchParams(window.location.search);
      return { filter: params.get("filter") || "", regex: params.get("regex") === "1", caseSensitive: params.get("case") === "1" };
    }
    function updateURLParams(filter, regex, caseSensitive) {
      const params = new URLSearchParams();
      if (filter) params.set("filter", filter);
      if (regex) params.set("regex", "1");
      if (caseSensitive) params.set("case", "1");
      history.replaceState(null, "", `${window.location.pathname}?${params.toString()}`);
    }
    function filterPlots() {
      const input = document.getElementById("filterInput").value.trim();
      const useRegex = document.getElementById("regexToggle").checked;
      const caseSensitive = document.getElementById("caseToggle").checked;
      const plots = document.getElementsByClassName("plot");
      updateURLParams(input, useRegex, caseSensitive);
      for (let plot of plots) {
        const name = plot.getAttribute("data-filename");
        let match = false;
        if (useRegex) {
          try {
            const regex = new RegExp(input, caseSensitive ? "" : "i");
            match = regex.test(name);
          } catch (e) { match = false; }
        } else {
          match = caseSensitive ? name.includes(input) : name.toLowerCase().includes(input.toLowerCase());
        }
        plot.style.display = match ? "" : "none";
      }
    }
    function clearFilter() {
      document.getElementById("filterInput").value = "";
      document.getElementById("regexToggle").checked = false;
      document.getElementById("caseToggle").checked = false;
      filterPlots();
    }
    window.addEventListener("DOMContentLoaded", () => {
      const { filter, regex, caseSensitive } = getQueryParams();
      document.getElementById("filterInput").value = filter;
      document.getElementById("regexToggle").checked = regex;
      document.getElementById("caseToggle").checked = caseSensitive;
      filterPlots();
    });
  </script>
</body>
</html>
"""
    os.makedirs(os.path.dirname(output_html) or '.', exist_ok=True)
    with open(output_html, "w") as f:
        f.write(html_header + html_body + html_footer)
    return os.path.abspath(output_html)


def generate_jsroot_html(canvas_keys, canvas_jsons, plots_per_row=4, output_html="view_plots.html", title="", intro_text=""):
    """
    Generate a self-contained HTML gallery that renders ROOT TCanvas objects via JSROOT.
    Canvas data is embedded as JSON directly in the HTML — no external file loading required,
    works with file:// and any HTTP server including GitHub Pages.

    Args:
        canvas_keys: Ordered list of canvas key names, or "NEWLINE" markers
        canvas_jsons: Dict {key: json_str} from ROOT.TBufferJSON.ToJSON(canvas).Data()
        plots_per_row: Number of plots per row
        output_html: Output HTML file path
        title: Page title
        intro_text: Optional intro text (supports **bold** and bullet lists)
    """
    if not title:
        title = os.path.splitext(os.path.basename(output_html))[0].replace('_', ' ')

    timestamp = datetime.now(ZoneInfo("Europe/Zurich")).strftime("%B %d, %Y, %I:%M %p %Z")

    processed_intro = re.sub(r'\*\*(.*?)\*\*', r'<b>\1</b>', intro_text)
    lines = processed_intro.split('\n')
    formatted_lines = []
    in_list = False
    for line in lines:
        stripped = line.strip()
        if stripped.startswith(('* ', '- ')):
            if not in_list:
                formatted_lines.append('<ul style="margin: 5px 0; padding-left: 20px;">')
                in_list = True
            formatted_lines.append(f'<li>{stripped[2:].strip()}</li>')
        else:
            if in_list:
                formatted_lines.append('</ul>')
                in_list = False
            if stripped.startswith('<'):
                formatted_lines.append(stripped)
            else:
                formatted_lines.append(line + "<br>")
    if in_list:
        formatted_lines.append('</ul>')
    formatted_intro = "".join(formatted_lines)

    real_keys = [k for k in canvas_keys if k != "NEWLINE"]
    key_to_idx = {k: i for i, k in enumerate(real_keys)}

    # Embed each canvas JSON in a <script type="application/json"> block.
    # These are inert (not executed), safe to embed raw JSON, and trivially
    # readable from JS via element.textContent.  The only escape needed is
    # "</script>" which cannot legally appear in JSON but we guard anyway.
    json_blocks = ""
    for idx, key in enumerate(real_keys):
        json_str = canvas_jsons.get(key, "{}")
        safe_json = json_str.replace("</script>", "<\\/script>")
        json_blocks += f'  <script type="application/json" id="jsdata_{idx}">{safe_json}</script>\n'

    def make_grid(keys):
        if not keys:
            return ""
        out = '    <div class="grid">\n'
        for idx, key in keys:
            out += f"""      <div class="plot" data-filename="{key}">
        <div class="plot-header">
          <div class="filename">{key}</div>
          <button class="pdf-btn" onclick="savePDF({idx}, '{key}')">PDF</button>
          <button class="pdf-btn" onclick="savePNG({idx}, '{key}')">PNG</button>
        </div>
        <div id="jsroot_{idx}" class="jsroot-canvas" data-idx="{idx}">
          <div class="jsroot-loading">Loading...</div>
        </div>
      </div>\n"""
        out += '    </div>\n'
        return out

    current_group = []
    body = ""
    for entry in canvas_keys:
        if entry == "NEWLINE":
            body += make_grid(current_group)
            current_group = []
        else:
            current_group.append((key_to_idx[entry], entry))
    body += make_grid(current_group)

    html = f"""<!DOCTYPE html>
<html>
<head>
  <meta charset="UTF-8">
  <title>{title}</title>
  <script src="https://cdn.jsdelivr.net/npm/jspdf@2.5.1/dist/jspdf.umd.min.js"></script>
  <script src="https://cdn.jsdelivr.net/npm/svg2pdf.js@2.2.3/dist/svg2pdf.umd.min.js"></script>
  <style>
    body {{ font-family: -apple-system, BlinkMacSystemFont, "Segoe UI", Roboto, sans-serif; padding: 10px 20px; margin: 0; line-height: 1.4; color: #333; }}
    .top-bar {{ display: flex; align-items: center; gap: 30px; padding: 10px 0 0; border-bottom: 1px solid #eee; margin-bottom: 15px; position: sticky; top: 0; background: white; z-index: 100; flex-wrap: wrap; }}
    .title-group {{ display: flex; flex-direction: column; min-width: fit-content; }}
    h1 {{ font-size: 20px; margin: 0; color: #1a1a1a; white-space: nowrap; }}
    .timestamp {{ font-size: 11px; color: #888; white-space: nowrap; }}
    .controls {{ display: flex; align-items: center; gap: 12px; padding-bottom: 8px; }}
    input[type="text"] {{ font-size: 14px; padding: 4px 10px; width: 250px; border: 1px solid #ccc; border-radius: 4px; }}
    .intro-text {{ background: #f8f9fa; border-left: 4px solid #007bff; padding: 10px 16px; margin-bottom: 15px; font-size: 14px; border-radius: 0 4px 4px 0; }}
    .intro-text li {{ margin-bottom: 4px; }}
    .grid {{ display: grid; grid-template-columns: repeat({plots_per_row}, 1fr); gap: 15px; }}
    .plot {{ border: 1px solid #e1e4e8; padding: 8px; text-align: center; border-radius: 6px; background: #fff; }}
    .plot-header {{ display: flex; justify-content: space-between; align-items: flex-start; margin-bottom: 6px; gap: 4px; }}
    .filename {{ font-weight: 600; font-size: 12px; word-break: break-all; color: #444; text-align: left; }}
    .pdf-btn {{ flex-shrink: 0; padding: 1px 7px; font-size: 11px; cursor: pointer; border: 1px solid #bbb; border-radius: 3px; background: #f5f5f5; color: #555; white-space: nowrap; }}
    .pdf-btn:hover {{ background: #e0e0e0; }}
    .jsroot-canvas {{ width: 100%; min-height: 200px; cursor: zoom-in; }}
    .jsroot-loading {{ display: flex; align-items: center; justify-content: center; height: 100%; color: #aaa; font-size: 13px; }}
    button {{ padding: 4px 12px; cursor: pointer; font-size: 13px; }}
    #enlargeModal {{ display: none; position: fixed; inset: 0; background: rgba(0,0,0,0.75); z-index: 9999; align-items: center; justify-content: center; }}
    #enlargeModal.open {{ display: flex; }}
    #enlargeOuter {{ display: flex; flex-direction: column; align-items: flex-end; gap: 6px; }}
    #enlargeBtnBar {{ display: flex; gap: 6px; }}
    #enlargeBtnBar button {{ padding: 3px 12px; font-size: 12px; cursor: pointer; border: 1px solid #bbb; border-radius: 3px; background: rgba(255,255,255,0.9); color: #333; }}
    #enlargeBtnBar button:hover {{ background: #fff; }}
    #enlargeInner {{ background: #fff; box-shadow: 0 8px 40px rgba(0,0,0,0.5); }}
    .axis-bar {{ display: flex; align-items: center; gap: 6px; flex-basis: 100%; padding: 5px 0 8px; font-size: 12px; color: #555; flex-wrap: wrap; }}
    .axis-bar input[type="number"] {{ width: 72px; padding: 2px 4px; font-size: 12px; border: 1px solid #ccc; border-radius: 3px; }}
    .axis-bar label {{ display: flex; align-items: center; gap: 3px; cursor: pointer; }}
    .axis-ctrl-btn {{ padding: 2px 8px; font-size: 12px; cursor: pointer; border: 1px solid #bbb; border-radius: 3px; background: #f5f5f5; color: #444; }}
    .axis-ctrl-btn:hover {{ background: #e0e0e0; }}
    .col-control {{ display: flex; align-items: center; gap: 4px; font-size: 13px; color: #555; white-space: nowrap; }}
    .col-btn {{ width: 24px; height: 24px; border: 1px solid #ccc; border-radius: 3px; cursor: pointer; font-size: 15px; line-height: 1; background: #f5f5f5; color: #333; display: flex; align-items: center; justify-content: center; padding: 0; }}
    .col-btn:hover {{ background: #e0e0e0; }}
    .col-label {{ min-width: 18px; text-align: center; font-weight: 600; }}
    @media (max-width: 800px) {{ .top-bar {{ gap: 10px; }} }}
  </style>
{json_blocks}</head>
<body>
  <div class="top-bar">
    <div class="title-group">
      <h1>{title}</h1>
      <div class="timestamp">{timestamp}</div>
    </div>
    <div class="controls">
      <input type="text" id="filterInput" placeholder="Filter by canvas name..." onkeyup="filterPlots()">
      <label><input type="checkbox" id="regexToggle" onchange="filterPlots()"> Regex</label>
      <label><input type="checkbox" id="caseToggle" onchange="filterPlots()"> Case</label>
      <button onclick="clearFilter()">Clear</button>
      <div class="col-control">
        Cols:
        <button class="col-btn" onclick="changeCols(1)" title="Fewer columns (larger plots)">+</button>
        <span class="col-label" id="colLabel">{plots_per_row}</span>
        <button class="col-btn" onclick="changeCols(-1)" title="More columns (smaller plots)">−</button>
      </div>
    </div>
    <div class="axis-bar">
      X&nbsp;<input type="number" id="g_xmin" placeholder="xmin" step="any" onkeydown="if(event.key==='Enter')applyGlobalRange()">
      &ndash;<input type="number" id="g_xmax" placeholder="xmax" step="any" onkeydown="if(event.key==='Enter')applyGlobalRange()">
      <label><input type="checkbox" id="g_logx">&nbsp;log</label>
      &emsp;Y&nbsp;<input type="number" id="g_ymin" placeholder="ymin" step="any" onkeydown="if(event.key==='Enter')applyGlobalRange()">
      &ndash;<input type="number" id="g_ymax" placeholder="ymax" step="any" onkeydown="if(event.key==='Enter')applyGlobalRange()">
      <label><input type="checkbox" id="g_logy">&nbsp;log</label>
      &emsp;Z&nbsp;<input type="number" id="g_zmin" placeholder="zmin" step="any" onkeydown="if(event.key==='Enter')applyGlobalRange()">
      &ndash;<input type="number" id="g_zmax" placeholder="zmax" step="any" onkeydown="if(event.key==='Enter')applyGlobalRange()">
      <label><input type="checkbox" id="g_logz">&nbsp;log</label>
      <button class="axis-ctrl-btn" onclick="applyGlobalRange()">Apply all</button>
      <button class="axis-ctrl-btn" onclick="resetGlobalRange()">Reset all</button>
    </div>
  </div>
  {f'<div class="intro-text">{formatted_intro}</div>' if intro_text else ''}
  <div id="plotContainer">
{body}  </div>

  <!-- Enlarge modal: intercepts JSROOT dblclick, draws at correct aspect ratio -->
  <div id="enlargeModal">
    <div id="enlargeOuter">
      <div id="enlargeBtnBar">
        <button onclick="saveEnlargePDF()">PDF</button>
        <button onclick="saveEnlargePNG()">PNG</button>
        <button onclick="closeEnlarge()" title="Close (Esc)">&#x2715;</button>
      </div>
      <div id="enlargeInner"></div>
    </div>
  </div>

  <script type="module">
    import {{ parse, draw, resize, makeImage }} from '{JSROOT_CDN}';

    const painters = new Map();

    async function loadCanvas(div) {{
      if (div.dataset.loaded) return;
      div.dataset.loaded = '1';
      const jsonEl = document.getElementById('jsdata_' + div.dataset.idx);
      if (!jsonEl) {{
        div.innerHTML = '<span style="color:red;font-size:12px">No data</span>';
        return;
      }}
      try {{
        const obj = parse(jsonEl.textContent);
        // Set per-canvas aspect ratio from the ROOT canvas dimensions so
        // the y/x ratio is preserved when columns (and thus div width) change.
        const cw = obj.fWindowWidth || obj.fCw || 4;
        const ch = obj.fWindowHeight || obj.fCh || 3;
        div.style.aspectRatio = cw + '/' + ch;
        div.innerHTML = '';
        const painter = await draw(div, obj, '');
        painters.set(parseInt(div.dataset.idx), painter);
        // Auto-apply any global range/log that was set before this canvas loaded
        if (_globalRangeActive()) await _applyRangeToPainter(painter);
      }} catch(e) {{
        div.innerHTML = '<span style="color:red;font-size:12px">Error: ' + e + '</span>';
      }}
    }}

    // Return a cleaned clone of the live JSROOT SVG plus its natural dimensions.
    // Strips interactive UI nodes and paints the bottom-left margin white to
    // cover the JSROOT colour-select indicator (lives in a <g transform=...> so
    // element-removal heuristics miss it).
    function _liveSvgClone(div, origW, origH) {{
      const liveSvg = div.querySelector('svg');
      if (!liveSvg) return null;
      const clone = liveSvg.cloneNode(true);
      clone.querySelectorAll('[class*="jsroot_interactive"]').forEach(el => el.remove());
      clone.querySelectorAll('foreignObject').forEach(el => el.remove());
      const curW = parseFloat(liveSvg.getAttribute('width')) || origW;
      const curH = parseFloat(liveSvg.getAttribute('height')) || origH;
      // Use viewBox coordinate space if present so the cover rect lands correctly.
      let svgW = curW, svgH = curH;
      const vb = liveSvg.getAttribute('viewBox');
      if (vb) {{
        const parts = vb.trim().split(/[\s,]+/);
        if (parts.length >= 4) {{ svgW = parseFloat(parts[2]); svgH = parseFloat(parts[3]); }}
      }}
      const ns = 'http://www.w3.org/2000/svg';
      const cover = document.createElementNS(ns, 'rect');
      cover.setAttribute('x', '0');
      cover.setAttribute('y', String(svgH * 0.93));
      cover.setAttribute('width', String(svgW * 0.10));
      cover.setAttribute('height', String(svgH * 0.07));
      cover.setAttribute('fill', '#ffffff');
      cover.setAttribute('fill-opacity', '1');
      cover.setAttribute('stroke', 'none');
      clone.appendChild(cover);
      if (!clone.getAttribute('viewBox')) clone.setAttribute('viewBox', `0 0 ${{curW}} ${{curH}}`);
      return {{ clone, curW, curH }};
    }}

    // Build a scaled SVG Blob from the cleaned live DOM rendering.
    function _liveSvgBlob(div, origW, origH, renderW, renderH) {{
      const result = _liveSvgClone(div, origW, origH);
      if (!result) return null;
      const {{ clone }} = result;
      clone.setAttribute('width', renderW);
      clone.setAttribute('height', renderH);
      return new Blob([new XMLSerializer().serializeToString(clone)], {{ type: 'image/svg+xml' }});
    }}

    // PDF export: true vector PDF via svg2pdf.js applied to the live DOM SVG.
    // The browser's rendering is already correct; we avoid JSROOT's SVG backend
    // which caused text-positioning bugs in earlier approaches.
    window.savePDF = async function(idx, key) {{
      const div = document.getElementById('jsroot_' + idx);
      if (div && !div.dataset.loaded) await loadCanvas(div);
      const jsonEl = document.getElementById('jsdata_' + idx);
      if (!jsonEl) {{ alert('No canvas data.'); return; }}
      try {{
        const obj = parse(jsonEl.textContent);
        const origW = obj.fWindowWidth || obj.fCw || 800;
        const origH = obj.fWindowHeight || obj.fCh || 600;
        const result = _liveSvgClone(div, origW, origH);
        if (!result) {{ alert('Canvas SVG not found — scroll to load it first.'); return; }}
        const {{ clone, curW, curH }} = result;
        clone.setAttribute('width', origW);
        clone.setAttribute('height', origH);
        // Attach off-screen so svg2pdf.js can resolve inherited styles/fonts.
        const host = document.createElement('div');
        host.style.cssText = 'position:absolute;left:-9999px;top:-9999px;width:0;height:0;overflow:hidden';
        host.appendChild(clone);
        document.body.appendChild(host);
        try {{
          const {{ jsPDF }} = window.jspdf;
          const pdf = new jsPDF({{
            orientation: origW > origH ? 'landscape' : 'portrait',
            unit: 'pt',
            format: [origW, origH]
          }});
          const {{ svg2pdf }} = window.svg2pdf;
          await svg2pdf(clone, pdf, {{ x: 0, y: 0, width: origW, height: origH }});
          pdf.save(key + '.pdf');
        }} finally {{
          document.body.removeChild(host);
        }}
      }} catch(e) {{ alert('PDF export failed: ' + e); }}
    }};

    window.savePNG = async function(idx, key) {{
      const div = document.getElementById('jsroot_' + idx);
      if (div && !div.dataset.loaded) await loadCanvas(div);
      const jsonEl = document.getElementById('jsdata_' + idx);
      if (!jsonEl) {{ alert('No canvas data.'); return; }}
      try {{
        const obj = parse(jsonEl.textContent);
        const origW = obj.fWindowWidth || obj.fCw || 800;
        const origH = obj.fWindowHeight || obj.fCh || 600;
        const scale = 2;
        const blob = _liveSvgBlob(div, origW, origH, origW * scale, origH * scale);
        if (!blob) {{ alert('Canvas SVG not found — scroll to load it first.'); return; }}
        const url = URL.createObjectURL(blob);
        await new Promise((resolve, reject) => {{
          const img = new Image();
          img.onload = () => {{
            const canvas = document.createElement('canvas');
            canvas.width = origW * scale; canvas.height = origH * scale;
            const ctx = canvas.getContext('2d');
            ctx.fillStyle = 'white'; ctx.fillRect(0, 0, canvas.width, canvas.height);
            ctx.drawImage(img, 0, 0);
            // Erase any JSROOT corner indicator that survived the SVG cleanup.
            ctx.fillStyle = 'white';
            ctx.fillRect(0, Math.floor(origH * 0.93) * scale, Math.ceil(origW * 0.10) * scale, Math.ceil(origH * 0.07) * scale);
            URL.revokeObjectURL(url);
            const a = document.createElement('a');
            a.href = canvas.toDataURL('image/png'); a.download = key + '.png';
            document.body.appendChild(a); a.click(); document.body.removeChild(a);
            resolve();
          }};
          img.onerror = reject;
          img.src = url;
        }});
      }} catch(e) {{ alert('PNG export failed: ' + e); }}
    }};

    // --- Enlarge modal: draw at correct aspect ratio ---
    let _enlargePainter = null;
    let _enlargeIdx = null;
    let _enlargeKey = null;

    window.showEnlarge = async function(idx, key) {{
      _enlargeIdx = idx;
      _enlargeKey = key;
      const jsonEl = document.getElementById('jsdata_' + idx);
      if (!jsonEl) return;
      const obj = parse(jsonEl.textContent);
      const cw = obj.fWindowWidth || obj.fCw || 800;
      const ch = obj.fWindowHeight || obj.fCh || 600;

      // Fit within 90vw × 88vh while preserving ratio
      const maxW = window.innerWidth  * 0.90;
      const maxH = window.innerHeight * 0.88;
      let dispW, dispH;
      if (maxW / maxH > cw / ch) {{
        dispH = maxH; dispW = maxH * cw / ch;
      }} else {{
        dispW = maxW; dispH = maxW * ch / cw;
      }}

      const inner = document.getElementById('enlargeInner');
      inner.innerHTML = '';
      const canvasDiv = document.createElement('div');
      canvasDiv.className = 'enlarge-jsroot';
      canvasDiv.style.cssText = `width:${{dispW}}px;height:${{dispH}}px;`;
      inner.style.width  = dispW + 'px';
      inner.style.height = dispH + 'px';
      inner.appendChild(canvasDiv);

      document.getElementById('enlargeModal').classList.add('open');
      _enlargePainter = await draw(canvasDiv, obj, '');
    }};

    window.closeEnlarge = function() {{
      document.getElementById('enlargeModal').classList.remove('open');
      document.getElementById('enlargeInner').innerHTML = '';
      _enlargePainter = null;
      _enlargeIdx = null;
      _enlargeKey = null;
    }};

    // Save from the enlarged modal canvas (full-size rendering = titles always visible)
    window.saveEnlargePDF = async function() {{
      if (_enlargeIdx === null) return;
      const div = document.getElementById('enlargeInner').querySelector('.enlarge-jsroot');
      if (!div) {{ alert('Canvas not ready.'); return; }}
      const jsonEl = document.getElementById('jsdata_' + _enlargeIdx);
      if (!jsonEl) return;
      try {{
        const obj = parse(jsonEl.textContent);
        const origW = obj.fWindowWidth || obj.fCw || 800;
        const origH = obj.fWindowHeight || obj.fCh || 600;
        const result = _liveSvgClone(div, origW, origH);
        if (!result) {{ alert('SVG not found.'); return; }}
        const {{ clone }} = result;
        clone.setAttribute('width', origW);
        clone.setAttribute('height', origH);
        const host = document.createElement('div');
        host.style.cssText = 'position:absolute;left:-9999px;top:-9999px;width:0;height:0;overflow:hidden';
        host.appendChild(clone);
        document.body.appendChild(host);
        try {{
          const {{ jsPDF }} = window.jspdf;
          const pdf = new jsPDF({{
            orientation: origW > origH ? 'landscape' : 'portrait',
            unit: 'pt', format: [origW, origH]
          }});
          const {{ svg2pdf }} = window.svg2pdf;
          await svg2pdf(clone, pdf, {{ x: 0, y: 0, width: origW, height: origH }});
          pdf.save((_enlargeKey || 'plot') + '.pdf');
        }} finally {{ document.body.removeChild(host); }}
      }} catch(e) {{ alert('PDF export failed: ' + e); }}
    }};

    window.saveEnlargePNG = async function() {{
      if (_enlargeIdx === null) return;
      const div = document.getElementById('enlargeInner').querySelector('.enlarge-jsroot');
      if (!div) {{ alert('Canvas not ready.'); return; }}
      const jsonEl = document.getElementById('jsdata_' + _enlargeIdx);
      if (!jsonEl) return;
      try {{
        const obj = parse(jsonEl.textContent);
        const origW = obj.fWindowWidth || obj.fCw || 800;
        const origH = obj.fWindowHeight || obj.fCh || 600;
        const scale = 2;
        const blob = _liveSvgBlob(div, origW, origH, origW * scale, origH * scale);
        if (!blob) {{ alert('SVG not found.'); return; }}
        const url = URL.createObjectURL(blob);
        await new Promise((resolve, reject) => {{
          const img = new Image();
          img.onload = () => {{
            const canvas = document.createElement('canvas');
            canvas.width = origW * scale; canvas.height = origH * scale;
            const ctx = canvas.getContext('2d');
            ctx.fillStyle = 'white'; ctx.fillRect(0, 0, canvas.width, canvas.height);
            ctx.drawImage(img, 0, 0);
            ctx.fillStyle = 'white';
            ctx.fillRect(0, Math.floor(origH * 0.93) * scale, Math.ceil(origW * 0.10) * scale, Math.ceil(origH * 0.07) * scale);
            URL.revokeObjectURL(url);
            const a = document.createElement('a');
            a.href = canvas.toDataURL('image/png');
            a.download = (_enlargeKey || 'plot') + '.png';
            document.body.appendChild(a); a.click(); document.body.removeChild(a);
            resolve();
          }};
          img.onerror = reject;
          img.src = url;
        }});
      }} catch(e) {{ alert('PNG export failed: ' + e); }}
    }};

    // Close on backdrop click or Escape
    document.getElementById('enlargeModal').addEventListener('click', function(e) {{
      if (e.target === this) closeEnlarge();
    }});
    document.addEventListener('keydown', function(e) {{
      if (e.key === 'Escape') closeEnlarge();
    }});

    const observer = new IntersectionObserver((entries) => {{
      for (const entry of entries) {{
        if (entry.isIntersecting) loadCanvas(entry.target);
      }}
    }}, {{ rootMargin: '300px' }});

    document.querySelectorAll('.jsroot-canvas').forEach(div => {{
      observer.observe(div);
      // Intercept dblclick in capture phase so JSROOT's own handler never fires
      div.addEventListener('dblclick', function(e) {{
        e.stopPropagation();
        e.preventDefault();
        const plot = div.closest('.plot');
        const key = plot ? plot.getAttribute('data-filename') : '';
        showEnlarge(parseInt(div.dataset.idx), key);
      }}, true);
    }});

    function getQueryParams() {{
      const p = new URLSearchParams(window.location.search);
      return {{ filter: p.get('filter') || '', regex: p.get('regex') === '1', caseSensitive: p.get('case') === '1' }};
    }}
    function updateURLParams(filter, regex, cs) {{
      const p = new URLSearchParams();
      if (filter) p.set('filter', filter);
      if (regex) p.set('regex', '1');
      if (cs) p.set('case', '1');
      history.replaceState(null, '', window.location.pathname + '?' + p.toString());
    }}

    window.filterPlots = function() {{
      const input = document.getElementById('filterInput').value.trim();
      const useRegex = document.getElementById('regexToggle').checked;
      const cs = document.getElementById('caseToggle').checked;
      updateURLParams(input, useRegex, cs);
      for (const plot of document.getElementsByClassName('plot')) {{
        const name = plot.getAttribute('data-filename');
        let match = false;
        if (useRegex) {{
          try {{ match = new RegExp(input, cs ? '' : 'i').test(name); }} catch(e) {{ match = false; }}
        }} else {{
          match = cs ? name.includes(input) : name.toLowerCase().includes(input.toLowerCase());
        }}
        plot.style.display = match ? '' : 'none';
        if (match) {{
          const jsDiv = plot.querySelector('.jsroot-canvas');
          if (jsDiv && jsDiv.dataset.loaded) resize(jsDiv);
        }}
      }}
    }};

    window.clearFilter = function() {{
      document.getElementById('filterInput').value = '';
      document.getElementById('regexToggle').checked = false;
      document.getElementById('caseToggle').checked = false;
      window.filterPlots();
    }};

    function _globalRangeActive() {{
      const ids = ['g_xmin','g_xmax','g_ymin','g_ymax','g_zmin','g_zmax'];
      return ids.some(id => document.getElementById(id)?.value.trim()) ||
             ['g_logx','g_logy','g_logz'].some(id => document.getElementById(id)?.checked);
    }}

    async function _applyRangeToPainter(painter) {{
      const fp = painter.getFramePainter?.();
      if (!fp) return;
      const get = id => {{ const v = parseFloat(document.getElementById(id)?.value); return isNaN(v) ? undefined : v; }};
      const logx = !!document.getElementById('g_logx')?.checked;
      const logy = !!document.getElementById('g_logy')?.checked;
      const logz = !!document.getElementById('g_logz')?.checked;
      fp.logx = logx; fp.logy = logy; fp.logz = logz;
      const pp = painter.getPadPainter?.() ?? painter;
      const rootPad = pp.getRootPad?.();
      if (rootPad) {{ rootPad.fLogx = logx ? 1 : 0; rootPad.fLogy = logy ? 1 : 0; rootPad.fLogz = logz ? 1 : 0; }}
      if (fp.zoom) {{
        await fp.zoom(
          get('g_xmin') ?? fp.scale_xmin, get('g_xmax') ?? fp.scale_xmax,
          get('g_ymin') ?? fp.scale_ymin, get('g_ymax') ?? fp.scale_ymax,
          get('g_zmin') ?? fp.scale_zmin, get('g_zmax') ?? fp.scale_zmax
        );
      }}
    }}

    window.applyGlobalRange = async function() {{
      for (const painter of painters.values()) await _applyRangeToPainter(painter);
    }};

    window.resetGlobalRange = async function() {{
      ['g_xmin','g_xmax','g_ymin','g_ymax','g_zmin','g_zmax'].forEach(id => {{ const el = document.getElementById(id); if (el) el.value = ''; }});
      document.getElementById('g_logx').checked = false;
      document.getElementById('g_logy').checked = false;
      document.getElementById('g_logz').checked = false;
      for (const painter of painters.values()) {{
        const fp = painter.getFramePainter?.();
        if (!fp) continue;
        const pp = painter.getPadPainter?.() ?? painter;
        const rootPad = pp.getRootPad?.();
        if (rootPad) {{ fp.logx = !!rootPad.fLogx; fp.logy = !!rootPad.fLogy; fp.logz = !!rootPad.fLogz; }}
        if (fp.unzoom) await fp.unzoom('xyz');
      }}
    }};

    let _currentCols = {plots_per_row};
    window.changeCols = function(delta) {{
      _currentCols = Math.max(1, Math.min(8, _currentCols - delta));
      document.getElementById('colLabel').textContent = _currentCols;
      document.querySelectorAll('.grid').forEach(g => {{
        g.style.gridTemplateColumns = `repeat(${{_currentCols}}, 1fr)`;
      }});
      // Re-render loaded canvases after layout settles
      setTimeout(() => {{
        document.querySelectorAll('.jsroot-canvas[data-loaded]').forEach(div => resize(div));
      }}, 50);
    }};

    const {{ filter, regex, caseSensitive }} = getQueryParams();
    document.getElementById('filterInput').value = filter;
    document.getElementById('regexToggle').checked = regex;
    document.getElementById('caseToggle').checked = caseSensitive;
    if (filter) window.filterPlots();
  </script>
</body>
</html>
"""

    os.makedirs(os.path.dirname(output_html) or '.', exist_ok=True)
    with open(output_html, "w") as f:
        f.write(html)
    return os.path.abspath(output_html)
