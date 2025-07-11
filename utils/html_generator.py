import os
from datetime import datetime
from zoneinfo import ZoneInfo


def generate_html(png_files, png_dir, plots_per_row=4, output_html="view_plots.html"):
    """
    Generate an HTML file to view PNG plots in a grid layout with filter, regex, deep linking, and Geneva time.

    Parameters:
    - png_files (list of str): list of PNG filenames (e.g., ['a.png', 'b.png'])
    - png_dir (str or Path): path to the directory containing the PNG files
    - plots_per_row (int): number of plots to show per row
    - output_html (str): path to output HTML file
    """

    html_dir = os.path.dirname(os.path.abspath(output_html))
    png_dir_abs = os.path.abspath(png_dir)
    output_html_abs = os.path.abspath(output_html)

    png_paths = [os.path.join(png_dir, f) for f in png_files]
    rel_paths = [os.path.relpath(png_path, start=html_dir)
                 for png_path in png_paths]

    timestamp = datetime.now(ZoneInfo("Europe/Zurich")
                             ).strftime("%B %d, %Y, %I:%M %p %Z")

    html_header = f"""<!DOCTYPE html>
<html>
<head>
  <meta charset="UTF-8">
  <title>PNG Plot Viewer</title>
  <style>
    body {{
      font-family: sans-serif;
      padding: 20px;
      margin: 0;
    }}

    h1 {{
      font-size: 20px;
      margin-bottom: 5px;
    }}

    .timestamp {{
      font-size: 15px;
      font-weight: bold;
      color: #444;
      margin-bottom: 20px;
    }}

    .controls {{
      position: sticky;
      top: 0;
      background: white;
      padding: 10px 0;
      z-index: 10;
      border-bottom: 1px solid #ccc;
      margin-bottom: 15px;
    }}

    .controls label {{
      margin-right: 15px;
      font-size: 14px;
    }}

    input[type="text"] {{
      font-size: 16px;
      padding: 5px;
      width: 400px;
      margin-right: 10px;
    }}

    button {{
      font-size: 14px;
      padding: 5px 10px;
      margin-left: 10px;
      cursor: pointer;
    }}

    .grid {{
      display: grid;
      grid-template-columns: repeat({plots_per_row}, 1fr);
      gap: 20px;
    }}

    .plot {{
      border: 1px solid #ccc;
      padding: 5px;
      text-align: center;
    }}

    .filename {{
      font-weight: bold;
      font-size: 14px;
      margin-bottom: 5px;
      word-break: break-word;
    }}

    img {{
      max-width: 100%;
      height: auto;
      cursor: pointer;
      transition: transform 0.2s;
    }}

    img:hover {{
      transform: scale(1.02);
    }}

    @media (max-width: 1200px) {{
      .grid {{
        grid-template-columns: repeat(2, 1fr);
      }}
    }}

    @media (max-width: 600px) {{
      .grid {{
        grid-template-columns: 1fr;
      }}
      input[type="text"] {{
        width: 100%;
        margin-bottom: 10px;
      }}
    }}
  </style>
</head>
<body>

  <h1>PNG Plot Viewer</h1>
  <div class="timestamp">Generated on: {timestamp}</div>

  <div class="controls">
    <input type="text" id="filterInput" placeholder="Filter by filename..." onkeyup="filterPlots()">
    <label><input type="checkbox" id="regexToggle" onchange="filterPlots()"> Use Wildcard/Regex</label>
    <label><input type="checkbox" id="caseToggle" onchange="filterPlots()"> Case Sensitive</label>
    <button onclick="clearFilter()">Clear Filter</button>
  </div>

  <div id="plotContainer" class="grid">
"""

    html_body = ""
    for filename, rel_path in zip(png_files, rel_paths):
        html_body += f"""    <div class="plot" data-filename="{filename}">
      <div class="filename">{filename}</div>
      <a href="{rel_path}" target="_blank">
        <img src="{rel_path}" alt="{filename}">
      </a>
    </div>
"""

    html_footer = """  </div>

  <script>
    function getQueryParams() {
      const params = new URLSearchParams(window.location.search);
      return {
        filter: params.get("filter") || "",
        regex: params.get("regex") === "1",
        caseSensitive: params.get("case") === "1",
      };
    }

    function updateURLParams(filter, regex, caseSensitive) {
      const params = new URLSearchParams();
      if (filter) params.set("filter", filter);
      if (regex) params.set("regex", "1");
      if (caseSensitive) params.set("case", "1");
      const newUrl = `${window.location.pathname}?${params.toString()}`;
      history.replaceState(null, "", newUrl);
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
            const escaped = input.replace(/[-\\/\\^$+?.()|[\\]{{}}]/g, '\\\\$&')
                                 .replace(/\\*/g, '.*')
                                 .replace(/\\?/g, '.');
            const flags = caseSensitive ? "" : "i";
            const regex = new RegExp(escaped, flags);
            match = regex.test(name);
          } catch (e) {
            match = false;
          }
        } else {
          const nameToTest = caseSensitive ? name : name.toLowerCase();
          const inputToTest = caseSensitive ? input : input.toLowerCase();
          match = nameToTest.includes(inputToTest);
        }

        plot.style.display = match ? "" : "none";
      }
    }

    function clearFilter() {
      document.getElementById("filterInput").value = "";
      document.getElementById("regexToggle").checked = false;
      document.getElementById("caseToggle").checked = false;
      history.replaceState(null, "", window.location.pathname);
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
