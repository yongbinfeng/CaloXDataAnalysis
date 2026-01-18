import os
import re
from datetime import datetime
from zoneinfo import ZoneInfo


def generate_html(png_files, png_dir, plots_per_row=4, output_html="view_plots.html", title="PNG Plot Viewer", intro_text=""):
    html_dir = os.path.dirname(os.path.abspath(output_html))
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
