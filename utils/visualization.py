import ROOT
from plotting.my_function import DrawHistos


def visualizeFERSBoards(fersboards, valuemaps=None, suffix="", gain="HG", quartzOnly=0):
    # quartzOnly: 0 use both quartz and plastic channels; 1 use only quartz channels; 2 use only plastic channels
    xmax = 14
    xmin = -14
    ymax = 10
    ymin = -10
    h2_Cer = ROOT.TH2D(f"h_FERSBoards_Cer_{suffix}", f"FERS Board Cer Channels for {suffix}",
                       int(xmax - xmin), xmin, xmax, int(ymax - ymin), ymin, ymax)
    h2_Sci = ROOT.TH2D(f"h_FERSBoards_Sci_{suffix}", f"FERS Board Sci Channels for {suffix}",
                       int(xmax - xmin), xmin, xmax, int(ymax - ymin), ymin, ymax)
    h2_Cer_3mm = ROOT.TH2D(f"h_FERSBoards_Cer_3mm_{suffix}", f"FERS Board Cer Channels for {suffix} (3mm)",
                           int(xmax - xmin), xmin, xmax, int(ymax - ymin) * 4, ymin, ymax)
    h2_Sci_3mm = ROOT.TH2D(f"h_FERSBoards_Sci_3mm_{suffix}", f"FERS Board Sci Channels for {suffix} (3mm)",
                           int(xmax - xmin), xmin, xmax, int(ymax - ymin) * 4, ymin, ymax)

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
                h2_Cer_3mm.Fill(channel_Cer.i_tower_x, channel_Cer.i_tower_y,
                                cer_encoded)
                h2_Sci_3mm.Fill(channel_Sci.i_tower_x, channel_Sci.i_tower_y,
                                sci_encoded)
            else:
                h2_Cer.Fill(channel_Cer.i_tower_x, channel_Cer.i_tower_y,
                            cer_encoded)
                h2_Sci.Fill(channel_Sci.i_tower_x, channel_Sci.i_tower_y,
                            sci_encoded)

    h2_Cer.SetMarkerSize(0.60)
    h2_Sci.SetMarkerSize(0.60)
    h2_Cer_3mm.SetMarkerSize(0.60)
    h2_Sci_3mm.SetMarkerSize(0.60)

    return [h2_Cer, h2_Cer_3mm], [h2_Sci, h2_Sci_3mm]


def visualizeDRSBoards(drs_boards, valuemaps=None, suffix="", quartzOnly=0):
    xmax = 14
    xmin = -14
    ymax = 10
    ymin = -10

    h2_DRS_Cer = ROOT.TH2D(f"h_DRSBoards_Cer_{suffix}", f"DRS Board Channels for {suffix}",
                           int(xmax - xmin), xmin, xmax, int(ymax - ymin), ymin, ymax)
    h2_DRS_Cer_3mm = ROOT.TH2D(f"h_DRSBoards_Cer_3mm_{suffix}", f"DRS Board Channels for {suffix} (3mm)",
                               int(xmax - xmin), xmin, xmax, int(ymax - ymin) * 4, ymin, ymax)

    h2_DRS_Sci = ROOT.TH2D(f"h_DRSBoards_Sci_{suffix}", f"DRS Board Channels for {suffix}",
                           int(xmax - xmin), xmin, xmax, int(ymax - ymin), ymin, ymax)
    h2_DRS_Sci_3mm = ROOT.TH2D(f"h_DRSBoards_Sci_3mm_{suffix}", f"DRS Board Channels for {suffix} (3mm)",
                               int(xmax - xmin), xmin, xmax, int(ymax - ymin) * 4, ymin, ymax)

    for _, DRSBoard in drs_boards.items():
        for channel_Sci in DRSBoard.get_list_of_channels(isCer=False):
            sci_encoded = channel_Sci.board_no * 100 + \
                channel_Sci.group_no * 10 + channel_Sci.channel_no
            if sci_encoded == 0:
                sci_encoded = 0.001
            if valuemaps and channel_Sci.get_channel_name() in valuemaps:
                sci_encoded = valuemaps[channel_Sci.get_channel_name()]
            if not channel_Sci.is6mm:
                h2_DRS_Sci_3mm.Fill(channel_Sci.i_tower_x, channel_Sci.i_tower_y,
                                    sci_encoded)
            else:
                h2_DRS_Sci.Fill(channel_Sci.i_tower_x, channel_Sci.i_tower_y,
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
                h2_DRS_Cer_3mm.Fill(channel_Cer.i_tower_x, channel_Cer.i_tower_y,
                                    cer_encoded)
            else:
                h2_DRS_Cer.Fill(channel_Cer.i_tower_x,
                                channel_Cer.i_tower_y, cer_encoded)

    h2_DRS_Cer.SetMarkerSize(0.60)
    h2_DRS_Sci.SetMarkerSize(0.60)
    h2_DRS_Cer_3mm.SetMarkerSize(0.60)
    h2_DRS_Sci_3mm.SetMarkerSize(0.60)

    return [h2_DRS_Cer, h2_DRS_Cer_3mm], [h2_DRS_Sci, h2_DRS_Sci_3mm]


def makeEventDisplay(h2_6mm, h2_3mm, output_name, outdir, run_number, zmin, zmax, isCer=False, extraToDraw=None):
    xmax = 14
    xmin = -14
    ymax = 10
    ymin = -10
    W_ref = 1000
    H_ref = 1100
    if isCer:
        extra_text = "Cer"
    else:
        extra_text = "Sci"
    DrawHistos([h2_6mm, h2_3mm], "", xmin, xmax, "iX", ymin, ymax, "iY", output_name, dology=False, drawoptions=["col,text", "col,text"],
               outdir=outdir, doth2=True, W_ref=W_ref, H_ref=H_ref, extra_text=extra_text, run_number=run_number, zmin=zmin, zmax=zmax, extraToDraw=extraToDraw)
