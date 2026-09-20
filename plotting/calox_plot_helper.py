import ROOT
from typing import Optional, Union, List, Dict, Tuple, Any
from plotting.my_function import DrawHistos
from core.plot_manager import PlotManager
from utils.visualization import (DRS_XMIN_DISPLAY, DRS_XMAX_DISPLAY, DRS_W_REF,
                                  FERS_YMIN_DISPLAY, FERS_YMAX_DISPLAY, FERS_H_REF,
                                  get_drs_display_x)


def create_pave_text(
    x1: float = 0.20,
    y1: float = 0.65,
    x2: float = 0.60,
    y2: float = 0.90,
    font: int = 42,
    size: float = 0.04,
    align: int = 11
) -> ROOT.TPaveText:
    """
    Create a standard TPaveText with common settings.

    This is a helper to reduce the boilerplate for creating text boxes.
    """
    pave = ROOT.TPaveText(x1, y1, x2, y2, "NDC")
    pave.SetTextAlign(align)
    pave.SetFillColorAlpha(0, 0)
    pave.SetBorderSize(0)
    pave.SetTextFont(font)
    pave.SetTextSize(size)
    return pave


def create_board_info_pave(
    board_no: int,
    tower_x: int,
    tower_y: int,
    channel_info: Optional[Dict[str, Any]] = None,
    position: Tuple[float, float, float, float] = (0.20, 0.65, 0.60, 0.90)
) -> ROOT.TPaveText:
    """
    Create a TPaveText with standard board/tower/channel information.

    Args:
        board_no: Board number
        tower_x: Tower X index
        tower_y: Tower Y index
        channel_info: Optional dict with channel details (e.g., {"Cer": 5, "Sci": 6})
        position: (x1, y1, x2, y2) for the text box
    """
    pave = create_pave_text(*position)
    pave.AddText(f"Board: {board_no}")
    pave.AddText(f"Tower X: {tower_x}")
    pave.AddText(f"Tower Y: {tower_y}")
    if channel_info:
        for key, value in channel_info.items():
            pave.AddText(f"{key} Channel: {value}")
    return pave


def _pct_range_from_hists(hists, lo=0.05, hi=0.90):
    """Return (zmin, zmax) from the lo–hi percentile of non-zero TH2 bin contents."""
    vals = []
    for h in hists:
        if h is None:
            continue
        for bx in range(1, h.GetNbinsX() + 1):
            for by in range(1, h.GetNbinsY() + 1):
                v = h.GetBinContent(bx, by)
                if v != 0.0:
                    vals.append(v)
    if not vals:
        return None, None
    vals.sort()
    n = len(vals)
    return vals[max(0, int(lo * n))], vals[min(n - 1, int(hi * n))]


class BoardPlotHelper:
    """
    Helper class for plotting board visualization histograms (2D maps).

    This handles the common pattern of plotting Cer/Sci × HG/LG × various stats
    with consistent styling.
    """

    def __init__(
        self,
        manager: PlotManager,
        xrange: Optional[Tuple[float, float]] = None,
        yrange: Tuple[float, float] = (FERS_YMIN_DISPLAY, FERS_YMAX_DISPLAY),
        W_ref: Optional[int] = None,
        H_ref: int = FERS_H_REF,
        run_number: Optional[int] = None,
    ):
        # If a run number is given, take the run-dependent DRS X range / W_ref
        # (e.g. widened to [-18, 10] for run >= 1896) unless explicitly overridden.
        if run_number is not None:
            rx0, rx1, rwref = get_drs_display_x(run_number)
            if xrange is None:
                xrange = (rx0, rx1)
            if W_ref is None:
                W_ref = rwref
        if xrange is None:
            xrange = (DRS_XMIN_DISPLAY, DRS_XMAX_DISPLAY)
        if W_ref is None:
            W_ref = DRS_W_REF
        self.manager = manager
        self.xrange = xrange
        self.yrange = yrange
        self.W_ref = W_ref
        self.H_ref = H_ref

    def plot_board_map(
        self,
        hists: List[ROOT.TH2],
        output_name: str,
        extra_text: str = "",
        zmin: Optional[float] = None,
        zmax: Optional[float] = None,
        nTextDigits: int = 0,
        zlabel: str = "",
        **kwargs
    ) -> 'BoardPlotHelper':
        if zmin is None or zmax is None:
            auto_lo, auto_hi = _pct_range_from_hists(hists)
            if zmin is None:
                zmin = auto_lo
            if zmax is None:
                zmax = auto_hi
        pm = self.manager
        first_opt = "colz,text" if zlabel else "col,text"
        DrawHistos(
            hists, "",
            self.xrange[0], self.xrange[1], "X [cm]",
            self.yrange[0], self.yrange[1], "Y [cm]",
            output_name,
            dology=False,
            drawoptions=[first_opt, "col,text"],
            outdir=pm.get_output_dir(),
            doth2=True,
            W_ref=self.W_ref,
            H_ref=self.H_ref,
            extra_text=extra_text,
            run_number=pm.run_number,
            zmin=zmin,
            zmax=zmax,
            nTextDigits=nTextDigits,
            zlabel=zlabel,
            usePNG=not pm.use_jsroot,
            canvas_json_store=pm._canvas_jsons if pm.use_jsroot else None,
            **kwargs
        )
        pm.add_plot(output_name)
        return self

    def plot_cer_sci_pair(
        self,
        cer_hists: List[ROOT.TH2],
        sci_hists: List[ROOT.TH2],
        output_base: str,
        zmin: Optional[float] = None,
        zmax: Optional[float] = None,
        nTextDigits: int = 0,
        zlabel: str = "",
        **kwargs
    ) -> 'BoardPlotHelper':
        if zmin is None or zmax is None:
            auto_lo, auto_hi = _pct_range_from_hists(cer_hists + sci_hists)
            if zmin is None:
                zmin = auto_lo
            if zmax is None:
                zmax = auto_hi
        self.plot_board_map(
            cer_hists, f"{output_base}_Cer",
            extra_text="Cer", zmin=zmin, zmax=zmax, nTextDigits=nTextDigits,
            zlabel=zlabel, **kwargs
        )
        self.plot_board_map(
            sci_hists, f"{output_base}_Sci",
            extra_text="Sci", zmin=zmin, zmax=zmax, nTextDigits=nTextDigits,
            zlabel=zlabel, **kwargs
        )
        return self


# ---------------------------------------------------------------------------
# Per-fibre summary plot (one point per test fibre, colour per fibre type)
# ---------------------------------------------------------------------------

import os

FIBER_PALETTE = [
    ROOT.kBlue + 1, ROOT.kRed + 1, ROOT.kGreen + 2, ROOT.kOrange + 1,
    ROOT.kViolet + 2, ROOT.kCyan + 2, ROOT.kMagenta + 1, ROOT.kTeal + 3,
    ROOT.kAzure + 7, ROOT.kSpring + 6,
]
# ROOT marker styles: filled circle, square, triangle-up, diamond, triangle-down,
# open circle, open square, open triangle-up, star, open diamond
FIBER_MARKERS = [20, 21, 22, 33, 23, 24, 25, 26, 29, 27]
# marker used for the secondary series of plot_fiber_summary (a thin "x")
SECONDARY_MARKER = 5


def fiber_base_type(label):
    """Strip trailing _N suffix: 'Quartz400_1' → 'Quartz400', 'Sapphire' → 'Sapphire'."""
    import re
    return re.sub(r'_\d+$', '', label)


def plot_fiber_summary(values, ytitle, ymin, ymax, plot_name, pm, lumi,
                       secondary=None, logy=False, note=None):
    """One point per fiber (label -> (value, error)), one color+marker per base
    fiber type, with the fiber labels on the x axis.

    secondary: optional second label -> (value, error) map drawn as "x"
    markers in the same colours (e.g. out-of-time next to in-time).
    logy: draw the y axis logarithmic (ymin must then be > 0).
    note: optional short text drawn top-left, e.g. what the markers mean.

    Returns the plot name (without extension) so the caller can add it to pm,
    or None if there is nothing to draw.
    """
    if not values:
        return None
    secondary = secondary or {}

    # Sort by (base type, full label) so same-type channels are adjacent
    sorted_labels = sorted(values, key=lambda l: (fiber_base_type(l), l))
    types_ordered = list(dict.fromkeys(fiber_base_type(l) for l in sorted_labels))
    type_color    = {t: FIBER_PALETTE[i % len(FIBER_PALETTE)]
                     for i, t in enumerate(types_ordered)}
    type_marker   = {t: FIBER_MARKERS[i % len(FIBER_MARKERS)]
                     for i, t in enumerate(types_ordered)}

    # Group points by base type
    pts_by_type = {}
    for i, label in enumerate(sorted_labels):
        v, ev = values[label]
        pts_by_type.setdefault(fiber_base_type(label), []).append((i + 1, v, ev))

    n = len(sorted_labels)

    c = ROOT.TCanvas(f"c_{plot_name}", "", 1400, 650)
    c.SetBottomMargin(0.28)
    c.SetLeftMargin(0.09)
    c.SetRightMargin(0.04)
    c.SetTopMargin(0.08)
    if logy:
        c.SetLogy(True)

    # Frame histogram — only purpose is to carry x-axis bin labels
    frame = ROOT.TH1F(f"h_{plot_name}_frame", "", n, 0.5, n + 0.5)
    for i, label in enumerate(sorted_labels):
        frame.GetXaxis().SetBinLabel(i + 1, label)
    frame.GetXaxis().LabelsOption("v")
    frame.GetXaxis().SetLabelSize(0.042)
    frame.GetXaxis().SetTickLength(0.0)
    frame.GetYaxis().SetTitle(ytitle)
    frame.GetYaxis().SetTitleOffset(0.9)
    frame.GetYaxis().SetTitleSize(0.048)
    frame.SetMinimum(ymin)
    frame.SetMaximum(ymax)
    frame.SetStats(0)
    frame.Draw("AXIS")

    # One TGraphErrors per base type
    graphs = []
    leg = ROOT.TLegend(0.10, 0.84, 0.90, 0.90)
    leg.SetNColumns(min(len(types_ordered), 8))
    leg.SetBorderSize(0)
    leg.SetFillStyle(0)
    leg.SetTextSize(0.036)

    for ft in types_ordered:
        pts = pts_by_type[ft]
        gr  = ROOT.TGraphErrors(len(pts))
        col = type_color[ft]
        for j, (x, v, ev) in enumerate(pts):
            gr.SetPoint(j, x, v)
            gr.SetPointError(j, 0.0, ev)
        gr.SetMarkerStyle(type_marker[ft])
        gr.SetMarkerSize(1.3)
        gr.SetMarkerColor(col)
        gr.SetLineColor(col)
        gr.SetLineWidth(2)
        gr.Draw("P SAME")
        leg.AddEntry(gr, ft, "p")
        graphs.append(gr)
        pts2 = [(i + 1, *secondary[l]) for i, l in enumerate(sorted_labels)
                if fiber_base_type(l) == ft and l in secondary]
        if pts2:
            gr2 = ROOT.TGraphErrors(len(pts2))
            for j, (x, v, ev) in enumerate(pts2):
                gr2.SetPoint(j, x, v)
                gr2.SetPointError(j, 0.0, ev)
            gr2.SetMarkerStyle(SECONDARY_MARKER)
            gr2.SetMarkerSize(1.6)
            gr2.SetMarkerColor(col)
            gr2.SetLineColor(col)
            gr2.SetLineWidth(2)
            gr2.Draw("P SAME")
            graphs.append(gr2)

    leg.Draw()

    lat = ROOT.TLatex()
    lat.SetNDC()
    lat.SetTextSize(0.038)
    lat.SetTextFont(42)
    lat.SetTextAlign(31)   # right-align
    lat.DrawLatex(0.96, 0.935, lumi)
    if note:
        lat.SetTextAlign(11)
        lat.DrawLatex(0.10, 0.935, note)

    if pm.use_jsroot:
        pm._canvas_jsons[plot_name] = ROOT.TBufferJSON.ToJSON(c).Data()
    else:
        c.SaveAs(os.path.join(pm.get_output_dir(), f"{plot_name}.png"))

    # Keep objects alive until after SaveAs
    c._keep = [frame, leg, lat] + graphs

    return plot_name
