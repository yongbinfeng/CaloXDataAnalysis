import ROOT
from plotting.my_function import DrawHistos


_X_SCALE = 1.2  # cm per bin in X
_Y_SCALE = 1.6  # cm per bin in Y
_NX = 28        # bins in X (-13.5 .. +13.5 half-integer centres)
_NY = 20        # bins in Y (-9.5 .. +9.5 half-integer centres)
_XMIN = -_NX / 2 * _X_SCALE   # -16.8 cm
_XMAX =  _NX / 2 * _X_SCALE   #  16.8 cm
_YMIN = -_NY / 2 * _Y_SCALE   # -16.0 cm
_YMAX =  _NY / 2 * _Y_SCALE   #  16.0 cm

# Full calorimeter face (FERS boards span the complete active area).
FERS_XMIN_DISPLAY = _XMIN   # -16.8 cm
FERS_XMAX_DISPLAY = _XMAX   #  16.8 cm
FERS_YMIN_DISPLAY = _YMIN   # -16.0 cm
FERS_YMAX_DISPLAY = _YMAX   #  16.0 cm
FERS_W_REF = 1270
FERS_H_REF = 1000

# DRS channels occupy only the central ±~4 cm in X; display in [-6, 6] cm.
# W_ref is scaled from the FERS reference so pixels/cm stays equal in X and Y.
DRS_XMIN_DISPLAY = -10.0
DRS_XMAX_DISPLAY =  10.0
DRS_W_REF = round(FERS_W_REF * (DRS_XMAX_DISPLAY - DRS_XMIN_DISPLAY) / (FERS_XMAX_DISPLAY - FERS_XMIN_DISPLAY))  # 454


def get_drs_display_x(run_number=None):
    """DRS board-map X display range and canvas W_ref (run-dependent).

    For run >= 1896 (PHASE_2 layout) the calorimeter extends further in -X, so
    widen X to [-18, 10] and scale W_ref proportionally so the cm-per-pixel
    (and thus the drawn cell size) stays the same as the default range.
    Returns (xmin, xmax, W_ref).
    """
    if run_number is not None and run_number >= 1896:
        xmin, xmax = -18.0, 10.0
    else:
        xmin, xmax = DRS_XMIN_DISPLAY, DRS_XMAX_DISPLAY
    w_ref = round(DRS_W_REF * (xmax - xmin)
                  / (DRS_XMAX_DISPLAY - DRS_XMIN_DISPLAY))
    return xmin, xmax, w_ref


def visualizeFERSBoards(fersboards, valuemaps=None, suffix="", gain="HG", quartzOnly=0):
    # quartzOnly: 0 use both quartz and plastic channels; 1 use only quartz channels; 2 use only plastic channels
    h2_Cer = ROOT.TH2D(f"h_FERSBoards_Cer_{suffix}", f"FERS Board Cer Channels for {suffix}",
                       _NX, _XMIN, _XMAX, _NY, _YMIN, _YMAX)
    h2_Sci = ROOT.TH2D(f"h_FERSBoards_Sci_{suffix}", f"FERS Board Sci Channels for {suffix}",
                       _NX, _XMIN, _XMAX, _NY, _YMIN, _YMAX)
    h2_Cer_3mm = ROOT.TH2D(f"h_FERSBoards_Cer_3mm_{suffix}", f"FERS Board Cer Channels for {suffix} (3mm)",
                           _NX, _XMIN, _XMAX, _NY * 4, _YMIN, _YMAX)
    h2_Sci_3mm = ROOT.TH2D(f"h_FERSBoards_Sci_3mm_{suffix}", f"FERS Board Sci Channels for {suffix} (3mm)",
                           _NX, _XMIN, _XMAX, _NY * 4, _YMIN, _YMAX)

    for FERSBoard in fersboards.values():
        for i_tower_x, i_tower_y in FERSBoard.get_list_of_towers():
            channel_Cer = FERSBoard.get_channel_by_tower(
                i_tower_x, i_tower_y, isCer=True)
            channel_Sci = FERSBoard.get_channel_by_tower(
                i_tower_x, i_tower_y, isCer=False)

            if quartzOnly == 1 and not channel_Cer.isQuartz:
                continue
            if quartzOnly == 2 and channel_Cer.isQuartz:
                continue

            channel_Cer_name = channel_Cer.get_channel_name(gain=gain)
            channel_Sci_name = channel_Sci.get_channel_name(gain=gain)

            # default encoding based on board and channel numbers
            cer_encoded = channel_Cer.board_no * 100 + channel_Cer.channel_no
            if cer_encoded == 0:
                cer_encoded = 0.001
            sci_encoded = channel_Sci.board_no * 100 + channel_Sci.channel_no
            if sci_encoded == 0:
                sci_encoded = 0.001

            if valuemaps and channel_Cer_name in valuemaps:
                cer_encoded = valuemaps[channel_Cer_name]
            if valuemaps and channel_Sci_name in valuemaps:
                sci_encoded = valuemaps[channel_Sci_name]
            if FERSBoard.is_3mm_size():
                h2_Cer_3mm.Fill(channel_Cer.i_tower_x * _X_SCALE, channel_Cer.i_tower_y * _Y_SCALE,
                                cer_encoded)
                h2_Sci_3mm.Fill(channel_Sci.i_tower_x * _X_SCALE, channel_Sci.i_tower_y * _Y_SCALE,
                                sci_encoded)
            else:
                h2_Cer.Fill(channel_Cer.i_tower_x * _X_SCALE, channel_Cer.i_tower_y * _Y_SCALE,
                            cer_encoded)
                h2_Sci.Fill(channel_Sci.i_tower_x * _X_SCALE, channel_Sci.i_tower_y * _Y_SCALE,
                            sci_encoded)

    h2_Cer.SetMarkerSize(0.60)
    h2_Sci.SetMarkerSize(0.60)
    h2_Cer_3mm.SetMarkerSize(0.60)
    h2_Sci_3mm.SetMarkerSize(0.60)

    return [h2_Cer, h2_Cer_3mm], [h2_Sci, h2_Sci_3mm]


def _drs_label(channel):
    """Numeric label drawn in a DRS channel-map cell.

    For bridge-numbered runs the digits encode bridge | board | channel, where
    the last two digits are the channel index on the board (group * 8 + channel,
    i.e. 0-31): e.g. Brg1 Board3 Group3 Channel7 -> 1331. Pre-bridge runs keep
    the old board * 100 + group * 10 + channel encoding.
    """
    if channel.bridge_no is not None:
        return (channel.bridge_no * 1000 + channel.board_no * 100
                + channel.group_no * 8 + channel.channel_no)
    return (channel.board_no * 100 + channel.group_no * 10
            + channel.channel_no)


def visualizeDRSBoards(drs_boards, valuemaps=None, suffix="", quartzOnly=0):
    h2_DRS_Cer = ROOT.TH2D(f"h_DRSBoards_Cer_{suffix}", f"DRS Board Channels for {suffix}",
                           _NX, _XMIN, _XMAX, _NY, _YMIN, _YMAX)
    h2_DRS_Cer_3mm = ROOT.TH2D(f"h_DRSBoards_Cer_3mm_{suffix}", f"DRS Board Channels for {suffix} (3mm)",
                               _NX, _XMIN, _XMAX, _NY * 4, _YMIN, _YMAX)

    h2_DRS_Sci = ROOT.TH2D(f"h_DRSBoards_Sci_{suffix}", f"DRS Board Channels for {suffix}",
                           _NX, _XMIN, _XMAX, _NY, _YMIN, _YMAX)
    h2_DRS_Sci_3mm = ROOT.TH2D(f"h_DRSBoards_Sci_3mm_{suffix}", f"DRS Board Channels for {suffix} (3mm)",
                               _NX, _XMIN, _XMAX, _NY * 4, _YMIN, _YMAX)

    for _, DRSBoard in drs_boards.items():
        for channel_Sci in DRSBoard.get_list_of_channels(isCer=False):
            ch_name = channel_Sci.get_channel_name()
            if valuemaps is not None:
                if ch_name not in valuemaps:
                    continue
                sci_encoded = valuemaps[ch_name]
            else:
                sci_encoded = _drs_label(channel_Sci)
                if sci_encoded == 0:
                    sci_encoded = 0.001
            if not channel_Sci.is6mm:
                h2_DRS_Sci_3mm.Fill(channel_Sci.i_tower_x * _X_SCALE, channel_Sci.i_tower_y * _Y_SCALE,
                                    sci_encoded)
            else:
                h2_DRS_Sci.Fill(channel_Sci.i_tower_x * _X_SCALE, channel_Sci.i_tower_y * _Y_SCALE,
                                sci_encoded)
        for channel_Cer in DRSBoard.get_list_of_channels(isCer=True):
            if quartzOnly == 1 and not channel_Cer.isQuartz:
                continue
            if quartzOnly == 2 and channel_Cer.isQuartz:
                continue
            ch_name = channel_Cer.get_channel_name()
            if valuemaps is not None:
                if ch_name not in valuemaps:
                    continue
                cer_encoded = valuemaps[ch_name]
            else:
                cer_encoded = _drs_label(channel_Cer)
                if cer_encoded == 0:
                    cer_encoded = 0.001
            if not channel_Cer.is6mm:
                h2_DRS_Cer_3mm.Fill(channel_Cer.i_tower_x * _X_SCALE, channel_Cer.i_tower_y * _Y_SCALE,
                                    cer_encoded)
            else:
                h2_DRS_Cer.Fill(channel_Cer.i_tower_x * _X_SCALE,
                                channel_Cer.i_tower_y * _Y_SCALE, cer_encoded)

    h2_DRS_Cer.SetMarkerSize(0.60)
    h2_DRS_Sci.SetMarkerSize(0.60)
    h2_DRS_Cer_3mm.SetMarkerSize(0.60)
    h2_DRS_Sci_3mm.SetMarkerSize(0.60)

    return [h2_DRS_Cer, h2_DRS_Cer_3mm], [h2_DRS_Sci, h2_DRS_Sci_3mm]


def makeEventDisplay(h2_6mm, h2_3mm, output_name, outdir, run_number, zmin, zmax, isCer=False, extraToDraw=None):
    extra_text = "Cer" if isCer else "Sci"
    DrawHistos([h2_6mm, h2_3mm], "", _XMIN, _XMAX, "X [cm]", _YMIN, _YMAX, "Y [cm]", output_name,
               dology=False, drawoptions=["col,text", "col,text"],
               outdir=outdir, doth2=True, W_ref=1270, H_ref=1000,
               extra_text=extra_text, run_number=run_number, zmin=zmin, zmax=zmax, extraToDraw=extraToDraw)
