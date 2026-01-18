from typing import Optional, Union, List, Dict, Tuple, Any
from dataclasses import dataclass


@dataclass
class PlotStyle:
    """Configuration for plot styling."""
    dology: bool = False
    dologx: bool = False
    dologz: bool = False
    drawoptions: Union[str, List[str]] = "HIST"
    mycolors: Optional[List[int]] = None
    addOverflow: bool = True
    addUnderflow: bool = True
    legendPos: Optional[List[float]] = None
    legendNCols: int = 1
    W_ref: int = 600
    H_ref: int = 600
    donormalize: bool = False
    extra_text: Optional[str] = None
    extraToDraw: Optional[Any] = None
    zmin: Optional[float] = None
    zmax: Optional[float] = None
    zlabel: str = "Events"
    linestyles: Optional[List[int]] = None

    def to_kwargs(self) -> Dict:
        """Convert to kwargs for DrawHistos."""
        kwargs = {
            'dology': self.dology,
            'dologx': self.dologx,
            'dologz': self.dologz,
            'drawoptions': self.drawoptions,
            'addOverflow': self.addOverflow,
            'addUnderflow': self.addUnderflow,
            'W_ref': self.W_ref,
            'H_ref': self.H_ref,
            'donormalize': self.donormalize,
            'zmin': self.zmin,
            'zmax': self.zmax,
            'zlabel': self.zlabel,
        }
        if self.mycolors is not None:
            kwargs['mycolors'] = self.mycolors
        if self.legendPos is not None:
            kwargs['legendPos'] = self.legendPos
        if self.legendNCols != 1:
            kwargs['legendNCols'] = self.legendNCols
        if self.extra_text is not None:
            kwargs['extra_text'] = self.extra_text
        if self.extraToDraw is not None:
            kwargs['extraToDraw'] = self.extraToDraw
        if self.linestyles is not None:
            kwargs['linestyles'] = self.linestyles
        return kwargs


# Pre-defined styles for common use cases
STYLE_1D_LINEAR = PlotStyle(dology=False, drawoptions="HIST")
STYLE_1D_LOG = PlotStyle(dology=True, drawoptions="HIST")
STYLE_2D_COLZ = PlotStyle(dology=False, drawoptions=[
                          "colz"], addOverflow=False, addUnderflow=False)
STYLE_CER = PlotStyle(dology=False, drawoptions="HIST", mycolors=[2])  # Red
STYLE_SCI = PlotStyle(dology=False, drawoptions="HIST", mycolors=[4])  # Blue
STYLE_CER_SCI = PlotStyle(dology=False, drawoptions="HIST", mycolors=[2, 4])
