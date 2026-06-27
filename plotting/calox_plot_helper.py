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
