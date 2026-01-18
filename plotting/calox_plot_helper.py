import ROOT
from typing import Optional, Union, List, Dict, Tuple, Any
from plotting.my_function import DrawHistos
from core.plot_manager import PlotManager


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


class BoardPlotHelper:
    """
    Helper class for plotting board visualization histograms (2D maps).

    This handles the common pattern of plotting Cer/Sci × HG/LG × various stats
    with consistent styling.
    """

    def __init__(
        self,
        manager: PlotManager,
        xrange: Tuple[float, float] = (-14, 14),
        yrange: Tuple[float, float] = (-10, 10),
        W_ref: int = 1000,
        H_ref: int = 1100
    ):
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
        zmin: float = 0,
        zmax: Optional[float] = None,
        nTextDigits: int = 0,
        **kwargs
    ) -> 'BoardPlotHelper':
        """
        Plot a board visualization map (2D histogram with text).

        Args:
            hists: List of 2D histograms (typically [main, 3mm])
            output_name: Output filename
            extra_text: Text label (e.g., "Cer", "Sci")
            zmin, zmax: Z-axis range
            nTextDigits: Number of digits for text display
        """
        DrawHistos(
            hists, "",
            self.xrange[0], self.xrange[1], "iX",
            self.yrange[0], self.yrange[1], "iY",
            output_name,
            dology=False,
            drawoptions=["col,text", "col,text"],
            outdir=self.manager.get_output_dir(),
            doth2=True,
            W_ref=self.W_ref,
            H_ref=self.H_ref,
            extra_text=extra_text,
            run_number=self.manager.run_number,
            zmin=zmin,
            zmax=zmax,
            nTextDigits=nTextDigits,
            **kwargs
        )
        self.manager.add_plot(output_name)
        return self

    def plot_cer_sci_pair(
        self,
        cer_hists: List[ROOT.TH2],
        sci_hists: List[ROOT.TH2],
        output_base: str,
        zmin: float = 0,
        zmax: Optional[float] = None,
        nTextDigits: int = 0,
        **kwargs
    ) -> 'BoardPlotHelper':
        """
        Plot a Cer/Sci pair of board maps.

        Args:
            cer_hists: Cerenkov histograms [main, 3mm]
            sci_hists: Scintillator histograms [main, 3mm]
            output_base: Base output name (will append _Cer and _Sci)
        """
        self.plot_board_map(
            cer_hists, f"{output_base}_Cer",
            extra_text="Cer", zmin=zmin, zmax=zmax, nTextDigits=nTextDigits, **kwargs
        )
        self.plot_board_map(
            sci_hists, f"{output_base}_Sci",
            extra_text="Sci", zmin=zmin, zmax=zmax, nTextDigits=nTextDigits, **kwargs
        )
        return self
