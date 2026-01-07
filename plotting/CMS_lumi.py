import ROOT as rt

# CMS_lumi
#   Initiated by: Gautier Hamel de Monchenault (Saclay)
#   Translated in Python by: Joshua Hardenbrook (Princeton)
#   Updated by:   Dinko Ferencek (Rutgers)
#

cms_text = "CaloX"
cms_text_font = 61

write_extra_text = True
extra_text = "Preliminary"
extra_text_font = 52

# lumiTextSize     = 0.6
lumi_text_size = 0.72
lumi_text_offset = 0.2

cms_text_size = 0.75
cms_text_offset = 0.1

rel_pos_x = 0.085
rel_pos_y = 0.035
rel_extra_dy = 1.2

extra_over_cms_text_size = 0.72/0.75

lumi_13TeV = "Cosmic runs"
# lumi_13TeV = "200 pb^{-1}"
lumi_5TeV = "298 pb^{-1}"
lumi_8TeV = "19.7 fb^{-1}"
lumi_7TeV = "5.1 fb^{-1}"
lumi_sqrtS = ""

draw_logo = False


def CMS_lumi(pad,  i_period,  i_pos_x, plot_cms=True):
    out_of_frame = False
    if (i_pos_x/10 == 0):
        out_of_frame = True

    alignY_ = 3
    alignX_ = 2
    if (i_pos_x/10 == 0):
        alignX_ = 1
    if (i_pos_x == 0):
        alignY_ = 1
    if (i_pos_x/10 == 1):
        alignX_ = 1
    if (i_pos_x/10 == 2):
        alignX_ = 2
    if (i_pos_x/10 == 3):
        alignX_ = 3
    align_ = 10*alignX_ + alignY_

    H = pad.GetWh()
    W = pad.GetWw()
    l = pad.GetLeftMargin()
    t = pad.GetTopMargin()
    # print "top margin %f"%t
    r = pad.GetRightMargin()
    b = pad.GetBottomMargin()
    e = 0.025

    pad.cd()

    lumi_text = ""
    if (i_period == 1):
        lumi_text += lumi_7TeV
        lumi_text += " (7 TeV)"
    elif (i_period == 2):
        lumi_text += lumi_8TeV
        lumi_text += " (8 TeV)"

    elif (i_period == 3):
        lumi_text = lumi_8TeV
        lumi_text += " (8 TeV)"
        lumi_text += " + "
        lumi_text += lumi_7TeV
        lumi_text += " (7 TeV)"
    elif (i_period == 4):
        lumi_text += lumi_13TeV
        # lumiText += " (13 TeV)"
    elif (i_period == 5):
        lumi_text += lumi_5TeV
        lumi_text += " (5 TeV)"
    elif (i_period == 7):
        if (out_of_frame):
            lumi_text += "#scale[0.85]{"
        lumi_text += lumi_13TeV
        # lumiText += " (13 TeV)"
        lumi_text += " + "
        lumi_text += lumi_8TeV
        lumi_text += " (8 TeV)"
        lumi_text += " + "
        lumi_text += lumi_7TeV
        lumi_text += " (7 TeV)"
        if (out_of_frame):
            lumi_text += "}"
    elif (i_period == 12):
        lumi_text += "8 TeV"
    elif (i_period == 0):
        lumi_text += lumi_sqrtS

    # print lumiText

    latex = rt.TLatex()
    latex.SetNDC()
    latex.SetTextAngle(0)
    latex.SetTextColor(rt.kBlack)

    extra_text_size = extra_over_cms_text_size*cms_text_size

    latex.SetTextFont(42)
    latex.SetTextAlign(31)
    latex.SetTextSize(lumi_text_size*t)

    latex.DrawLatex(1-r, 1-t+lumi_text_offset*t, lumi_text)

    if (out_of_frame and plot_cms):
        latex.SetTextFont(cms_text_font)
        latex.SetTextAlign(11)
        latex.SetTextSize(cms_text_size*t)
        latex.DrawLatex(l, 1-t+lumi_text_offset*t, cms_text)

    pad.cd()

    posX_ = 0
    if (i_pos_x % 10 <= 1):
        posX_ = l + rel_pos_x*(1-l-r)
    elif (i_pos_x % 10 == 2):
        posX_ = l + 0.5*(1-l-r)
    elif (i_pos_x % 10 == 3):
        posX_ = 1-r - rel_pos_x*(1-l-r)

    posY_ = 1-t - rel_pos_y*(1-t-b)

    if (not out_of_frame):
        if (draw_logo):
            posX_ = l + 0.045*(1-l-r)*W/H
            posY_ = 1-t - 0.045*(1-t-b)
            xl_0 = posX_
            yl_0 = posY_ - 0.15
            xl_1 = posX_ + 0.15*H/W
            yl_1 = posY_
            CMS_logo = rt.TASImage("CMS-BW-label.png")
            pad_logo = rt.TPad("logo", "logo", xl_0, yl_0, xl_1, yl_1)
            pad_logo.Draw()
            pad_logo.cd()
            CMS_logo.Draw("X")
            pad_logo.Modified()
            pad.cd()
        elif plot_cms:
            latex.SetTextFont(cms_text_font)
            latex.SetTextSize(cms_text_size*t)
            latex.SetTextAlign(align_)
            latex.DrawLatex(posX_, posY_, cms_text)
        if (write_extra_text):
            latex.SetTextFont(extra_text_font)
            latex.SetTextAlign(align_)
            latex.SetTextSize(extra_text_size*t)
            latex.DrawLatex(posX_+0.1, posY_ - rel_extra_dy *
                            cms_text_size*t, extra_text)
    elif (write_extra_text):
        if (i_pos_x == 0):
            # posX_ =   l +  relPosX*(1-l-r)
            posX_ = l
            posY_ = 1-t+lumi_text_offset*t

        latex.SetTextFont(extra_text_font)
        latex.SetTextSize(extra_text_size*t)
        latex.SetTextAlign(align_)
        latex.DrawLatex(posX_ + 0.13, posY_, extra_text)

    pad.Update()
