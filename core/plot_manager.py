import os
from typing import Optional, Union, List, Dict, Tuple, Any
import ROOT
from configs.plot_style import PlotStyle, STYLE_2D_COLZ, STYLE_CER_SCI
from plotting.my_function import DrawHistos
from utils.html_generator import generate_html


class HistogramNotFoundError(Exception):
    """Raised when a histogram is not found in a ROOT file."""
    pass


class PlotManager:
    """
    Manages histogram plotting workflows with reduced code duplication.
    """

    _generated_html_registry: List[str] = []

    def __init__(
        self,
        rootdir: str,
        plotdir: str,
        htmldir: str,
        run_number: Optional[int] = None,
        default_style: Optional[PlotStyle] = None,
    ):
        """
        Initialize the PlotManager.

        Args:
            rootdir: Directory containing ROOT files
            plotdir: Base directory for output plots
            htmldir: Base directory for HTML output
            run_number: Run number for labeling (optional)
            default_style: Default PlotStyle to use (optional)
        """
        self.rootdir = rootdir
        self.plotdir = plotdir
        self.htmldir = htmldir
        self.run_number = run_number
        self.default_style = default_style or PlotStyle()

        # Internal state
        self._plots: List[str] = []
        self._open_files: Dict[str, ROOT.TFile] = {}
        self._current_outdir: Optional[str] = None

    def reset_plots(self) -> 'PlotManager':
        """Clear the collected plots list."""
        self._plots = []
        return self

    def set_output_dir(self, subdir: str) -> 'PlotManager':
        """Set the current output subdirectory."""
        self._current_outdir = os.path.join(self.plotdir, subdir)
        return self

    def get_output_dir(self) -> str:
        """Get the current output directory."""
        if self._current_outdir is None:
            return self.plotdir
        return self._current_outdir

    def _get_file(self, filename: str) -> ROOT.TFile:
        """Get or open a ROOT file (cached)."""
        filepath = os.path.join(self.rootdir, filename) if not os.path.isabs(
            filename) else filename
        if filepath not in self._open_files:
            if not os.path.exists(filepath):
                raise FileNotFoundError(f"ROOT file not found: {filepath}")
            self._open_files[filepath] = ROOT.TFile(filepath, "READ")
        return self._open_files[filepath]

    def get_histogram(
        self,
        filename: str,
        hist_name: str,
        required: bool = True
    ) -> Optional[ROOT.TH1]:
        """
        Retrieve a histogram from a ROOT file.

        Args:
            filename: Name of the ROOT file (relative to rootdir)
            hist_name: Name of the histogram
            required: If True, raise error if not found; if False, return None

        Returns:
            The histogram, or None if not found and required=False
        """
        infile = self._get_file(filename)
        hist = infile.Get(hist_name)

        if not hist:
            if required:
                raise HistogramNotFoundError(
                    f"Histogram '{hist_name}' not found in {filename}"
                )
            print(f"Warning: Histogram '{hist_name}' not found in {filename}")
            return None

        return hist

    def get_histogram_by_pattern(
        self,
        filename: str,
        varname: str,
        suffix: str = "",
        prefix: str = "hist_",
        required: bool = True
    ) -> Optional[ROOT.TH1]:
        """
        Retrieve a histogram using the common naming pattern: prefix_varname_suffix.

        Args:
            filename: Name of the ROOT file
            varname: Variable name
            suffix: Suffix (e.g., "inclusive", "electron")
            prefix: Prefix (default "hist_")
            required: If True, raise error if not found

        Returns:
            The histogram, or None if not found and required=False
        """
        hist_name = f"{prefix}{varname}_{suffix}" if suffix else f"{prefix}{varname}"
        return self.get_histogram(filename, hist_name, required=required)

    def add_plot(self, output_name: str, prepend: bool = False) -> 'PlotManager':
        """Add a plot to the collection."""
        plot_file = output_name if output_name.endswith(
            '.png') else f"{output_name}.png"
        if prepend:
            self._plots.insert(0, plot_file)
        else:
            self._plots.append(plot_file)
        return self

    def add_newline(self) -> 'PlotManager':
        """Add a newline marker in the plot list (for HTML layout)."""
        self._plots.append("NEWLINE")
        return self

    def plot_1d(
        self,
        hist: Union[ROOT.TH1, List[ROOT.TH1]],
        output_name: str,
        xlabel: str,
        xrange: Tuple[float, float],
        ylabel: str = "Events",
        yrange: Tuple[Optional[float], Optional[float]] = (None, None),
        legends: Union[str, List[str]] = "",
        style: Optional[PlotStyle] = None,
        prepend: bool = False,
        includeRunNumber: bool = True,
        **kwargs
    ) -> 'PlotManager':
        """
        Draw a 1D histogram plot.

        Args:
            hist: Histogram or list of histograms
            output_name: Output filename (without extension)
            xlabel: X-axis label
            xrange: (xmin, xmax) tuple
            ylabel: Y-axis label
            yrange: (ymin, ymax) tuple, None for auto
            legends: Legend labels
            style: PlotStyle to use (or use default)
            prepend: If True, add to beginning of plot list
            includeRunNumber: If True, include run number in the plot
            **kwargs: Additional arguments passed to DrawHistos

        Returns:
            self for chaining
        """
        hists = [hist] if isinstance(hist, ROOT.TH1) else hist
        style = style or self.default_style

        # Merge style kwargs with explicit kwargs
        plot_kwargs = style.to_kwargs()
        plot_kwargs.update(kwargs)

        DrawHistos(
            hists,
            legends,
            xrange[0], xrange[1],
            xlabel,
            yrange[0], yrange[1],
            ylabel,
            output_name,
            outdir=self.get_output_dir(),
            run_number=self.run_number if includeRunNumber else None,
            **plot_kwargs
        )

        return self.add_plot(output_name, prepend=prepend)

    def plot_2d(
        self,
        hist: ROOT.TH2,
        output_name: str,
        xlabel: str,
        xrange: Tuple[float, float],
        ylabel: str,
        yrange: Tuple[float, float],
        style: Optional[PlotStyle] = None,
        prepend: bool = False,
        **kwargs
    ) -> 'PlotManager':
        """
        Draw a 2D histogram plot.

        Args:
            hist: 2D histogram
            output_name: Output filename (without extension)
            xlabel: X-axis label
            xrange: (xmin, xmax) tuple
            ylabel: Y-axis label
            yrange: (ymin, ymax) tuple
            style: PlotStyle to use
            prepend: If True, add to beginning of plot list
            **kwargs: Additional arguments passed to DrawHistos

        Returns:
            self for chaining
        """
        style = style or STYLE_2D_COLZ

        plot_kwargs = style.to_kwargs()
        plot_kwargs['doth2'] = True
        plot_kwargs.update(kwargs)

        DrawHistos(
            [hist],
            "",
            xrange[0], xrange[1],
            xlabel,
            yrange[0], yrange[1],
            ylabel,
            output_name,
            outdir=self.get_output_dir(),
            run_number=self.run_number,
            **plot_kwargs
        )

        return self.add_plot(output_name, prepend=prepend)

    def plot_from_file(
        self,
        filename: str,
        hist_name: str,
        output_name: str,
        xlabel: str,
        xrange: Tuple[float, float],
        ylabel: str = "Events",
        yrange: Tuple[Optional[float], Optional[float]] = (None, None),
        style: Optional[PlotStyle] = None,
        required: bool = False,
        prepend: bool = False,
        **kwargs
    ) -> 'PlotManager':
        """
        Load a histogram from file and plot it.

        Combines get_histogram and plot_1d in one call.

        Returns:
            self for chaining (skips plot if histogram not found and required=False)
        """
        hist = self.get_histogram(filename, hist_name, required=required)
        if hist is None:
            return self

        return self.plot_1d(
            hist, output_name, xlabel, xrange,
            ylabel=ylabel, yrange=yrange,
            style=style, prepend=prepend, **kwargs
        )

    def plot_cer_sci(
        self,
        filename: str,
        varname_cer: str,
        varname_sci: str,
        suffix: str,
        output_name: str,
        xlabel: str,
        xrange: Tuple[float, float],
        ylabel: str = "Events",
        yrange: Tuple[Optional[float], Optional[float]] = (None, None),
        style: Optional[PlotStyle] = None,
        prepend: bool = False,
        **kwargs
    ) -> 'PlotManager':
        """
        Plot Cer and Sci histograms together (common pattern).

        Args:
            filename: ROOT file name
            varname_cer: Cerenkov variable name
            varname_sci: Scintillator variable name
            suffix: Histogram suffix
            output_name: Output filename
            xlabel, xrange, ylabel, yrange: Axis configuration
            style: PlotStyle (defaults to red/blue for Cer/Sci)

        Returns:
            self for chaining
        """
        hist_cer = self.get_histogram_by_pattern(
            filename, varname_cer, suffix, required=False)
        hist_sci = self.get_histogram_by_pattern(
            filename, varname_sci, suffix, required=False)

        hists = [h for h in [hist_cer, hist_sci] if h is not None]
        if not hists:
            return self

        legends = []
        if hist_cer:
            legends.append("Cer")
        if hist_sci:
            legends.append("Sci")

        style = style or STYLE_CER_SCI
        # Adjust colors if only one histogram
        if len(hists) == 1:
            style = PlotStyle(
                **{**style.__dict__, 'mycolors': [2] if hist_cer else [4]}
            )

        return self.plot_1d(
            hists, output_name, xlabel, xrange,
            ylabel=ylabel, yrange=yrange,
            legends=legends, style=style, prepend=prepend, **kwargs
        )

    def generate_html(
        self,
        output_path: str,
        plots_per_row: int = 4,
        title: str = "PNG Plot Viewer",
        intro_text: str = "",
        clear_plots: bool = True
    ) -> str:
        """
        Generate an HTML gallery from collected plots.

        Args:
            output_path: Path relative to htmldir
            plots_per_row: Number of plots per row
            title: Page title
            intro_text: Introduction text
            clear_plots: If True, reset plots list after generation

        Returns:
            Full path to generated HTML file
        """
        full_path = os.path.join(self.htmldir, output_path)

        if full_path not in PlotManager._generated_html_registry:
            PlotManager._generated_html_registry.append(full_path)

        generate_html(
            self._plots,
            self.get_output_dir(),
            plots_per_row=plots_per_row,
            output_html=full_path,
            title=title,
            intro_text=intro_text
        )

        if clear_plots:
            self._plots = []

        return full_path

    def close(self):
        """Close all open ROOT files."""
        for f in self._open_files.values():
            f.Close()
        self._open_files = {}

    def __enter__(self) -> 'PlotManager':
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        self.close()
        return False

    @staticmethod
    def print_html_summary():
        """Print all generated HTML files."""
        """Prints all HTML files generated during this session."""
        print("\n" + "="*50)
        print("📊 GENERATED HTML REPORTS SUMMARY")
        print("="*50)
        if not PlotManager._generated_html_registry:
            print("No HTML files were generated.")
        for path in PlotManager._generated_html_registry:
            print(f"✅ {path}")
        print("="*50 + "\n")
