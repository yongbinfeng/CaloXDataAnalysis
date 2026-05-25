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
            sci_encoded = channel_Sci.board_no * 100 + \
                channel_Sci.group_no * 10 + channel_Sci.channel_no
            if sci_encoded == 0:
                sci_encoded = 0.001
            if valuemaps and channel_Sci.get_channel_name() in valuemaps:
                sci_encoded = valuemaps[channel_Sci.get_channel_name()]
            if not channel_Sci.is6mm:
                h2_DRS_Sci_3mm.Fill(channel_Sci.i_tower_x * _X_SCALE, channel_Sci.i_tower_y * _Y_SCALE,
                                    sci_encoded)
            else:
                h2_DRS_Sci.Fill(channel_Sci.i_tower_x * _X_SCALE, channel_Sci.i_tower_y * _Y_SCALE,
                                sci_encoded)
        for channel_Cer in DRSBoard.get_list_of_channels(isCer=True):
            cer_encoded = channel_Cer.board_no * 100 + \
                channel_Cer.group_no * 10 + channel_Cer.channel_no
            if quartzOnly == 1 and not channel_Cer.isQuartz:
                continue
            if quartzOnly == 2 and channel_Cer.isQuartz:
                continue
            if cer_encoded == 0:
                cer_encoded = 0.001
            if valuemaps and channel_Cer.get_channel_name() in valuemaps:
                cer_encoded = valuemaps[channel_Cer.get_channel_name()]
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
