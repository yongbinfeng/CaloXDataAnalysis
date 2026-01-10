import os
from datetime import datetime
from zoneinfo import ZoneInfo
from pathlib import Path


def generate_html(png_files, png_dir, plots_per_row=4, output_html="view_plots.html", intro_text=""):
    """
    Generate an HTML file to view PNG plots with an optional introductory text.
    Supports **bold** via markdown-style syntax.
    """
    html_dir = os.path.dirname(os.path.abspath(output_html))
    png_dir_abs = os.path.abspath(png_dir)
    output_html_abs = os.path.abspath(output_html)

    real_png_files = [f for f in png_files if f != "NEWLINE"]
    png_paths = [os.path.join(png_dir, f) for f in real_png_files]
    rel_paths = [os.path.relpath(p, start=html_dir) for p in png_paths]

    path_map = dict(zip(real_png_files, rel_paths))

    timestamp = datetime.now(ZoneInfo("Europe/Zurich")
                             ).strftime("%B %d, %Y, %I:%M %p %Z")

    # Format the intro text:
    # 1. Handle Newlines
    # 2. Convert **bold** to <b> tags
    formatted_intro = intro_text.replace("\n", "<br>")
    import re
    formatted_intro = re.sub(r'\*\*(.*?)\*\*', r'<b>\1</b>', formatted_intro)

    html_header = f"""<!DOCTYPE html>
<html>
<head>
  <meta charset="UTF-8">
  <title>PNG Plot Viewer</title>
  <style>
    body {{
      font-family: -apple-system, BlinkMacSystemFont, "Segoe UI", Roboto, sans-serif;
      padding: 30px;
      margin: 0;
      line-height: 1.6;
      color: #333;
    }}
    h1 {{
      font-size: 24px;
      margin-bottom: 5px;
      color: #1a1a1a;
    }}
    .timestamp {{
      font-size: 14px;
      color: #666;
      margin-bottom: 20px;
    }}
    .intro-text {{
      background: #f8f9fa;
      border-left: 5px solid #007bff;
      padding: 20px;
      margin-bottom: 30px;
      font-size: 16px;
      border-radius: 0 4px 4px 0;
      box-shadow: 0 1px 3px rgba(0,0,0,0.1);
    }}
    .controls {{
      position: sticky;
      top: 0;
      background: white;
      padding: 15px 0;
      z-index: 10;
      border-bottom: 1px solid #eee;
      margin-bottom: 25px;
    }}
    input[type="text"] {{
      font-size: 16px;
      padding: 8px 12px;
      width: 350px;
      border: 1px solid #ccc;
      border-radius: 4px;
    }}
    .grid {{
      display: grid;
      grid-template-columns: repeat({plots_per_row}, 1fr);
      gap: 25px;
    }}
    .plot {{
      border: 1px solid #e1e4e8;
      padding: 12px;
      text-align: center;
      border-radius: 6px;
      background: #fff;
      transition: box-shadow 0.2s;
    }}
    .plot:hover {{
      box-shadow: 0 4px 12px rgba(0,0,0,0.1);
    }}
    .filename {{
      font-weight: 600;
      font-size: 13px;
      margin-bottom: 10px;
      word-break: break-all;
      color: #444;
    }}
    img {{
      max-width: 100%;
      height: auto;
      border-radius: 2px;
    }}
    @media (max-width: 1200px) {{ .grid {{ grid-template-columns: repeat(2, 1fr); }} }}
    @media (max-width: 600px) {{ .grid {{ grid-template-columns: 1fr; }} }}
  </style>
</head>
<body>

  <h1>PNG Plot Viewer</h1>
  <div class="timestamp">Generated on: {timestamp}</div>

  {f'<div class="intro-text">{formatted_intro}</div>' if intro_text else ''}

  <div class="controls">
    <input type="text" id="filterInput" placeholder="Filter by filename..." onkeyup="filterPlots()">
    <label style="margin-left:15px"><input type="checkbox" id="regexToggle" onchange="filterPlots()"> Regex</label>
    <label style="margin-left:10px"><input type="checkbox" id="caseToggle" onchange="filterPlots()"> Case</label>
    <button onclick="clearFilter()" style="margin-left:20px; padding: 8px 12px; cursor:pointer;">Clear</button>
  </div>

  <div id="plotContainer">
"""

    # ... [Rest of the body logic remains the same] ...
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

    # ... [Script section remains the same as previous version] ...
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
            const escaped = input.replace(/[-\\/\\^$+?.()|[\\]{{}}]/g, '\\\\$&').replace(/\\*/g, '.*').replace(/\\?/g, '.');
            const regex = new RegExp(escaped, caseSensitive ? "" : "i");
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

    os.makedirs(os.path.dirname(output_html), exist_ok=True)
    with open(output_html, "w") as f:
        f.write(html_header + html_body + html_footer)

    print(f"✅ HTML viewer generated at: {output_html_abs}")

    return output_html_abs
