#!/usr/bin/env python

import ROOT as r
from plottingUtils import Print, MakeCumu, SetBatch
from plotDetails import mapRanges, mapDMRanges, alphaTDict
from array import array
import math


###-------------------------------------------------------------------###
###-------------------------------------------------------------------###

r.gROOT.SetBatch(r.kTRUE)
r.gStyle.SetOptStat(0)
r.TH1.SetDefaultSumw2(1)
r.TH2.SetDefaultSumw2(1)


settings = {
    "model":    ["T2cc", "T2", "T2_4body", "T2tt"][-1],
    "version":  5,
    "mode":     ["JES", "ISR", "bTag", "LeptonVeto", "DeadECAL", "MHT_MET", "3jet"],
    "inclHT":   [False, True][1],
    "HTBins":   ["200_275", "275_325", "325_375", "375_475", "475_575", "575_675", "675_775", "775_875", "875_975", "975_1075", "1075"],
    "deltaM":   [False, True][1],
    "jMulti":   ["le3j", "ge4j", "ge2j"][:2],
    "bMulti":   ["eq0b", "eq1b", "eq2b", "eq3b", "ge0b"][:4],
    "text_plot":[False, True][1]
}

###-------------------------------------------------------------------###

def set_palette(name="", ncontours=30):
    """Set a color palette from a given RGB list
    stops, red, green and blue should all be lists of the same length
    see set_decent_colors for an example"""

    if name == "gray" or name == "grayscale":
        stops = [0.00, 0.34, 0.61, 0.84, 1.00]
        red   = [1.00, 0.84, 0.61, 0.34, 0.00]
        green = [1.00, 0.84, 0.61, 0.34, 0.00]
        blue  = [1.00, 0.84, 0.61, 0.34, 0.00]
    # elif name == "whatever":
        # (define more palettes)
    else:
        # default palette, looks cool
        stops = [0.00, 0.34, 0.61, 0.84, 1.00]
        red   = [0.50, 0.00, 0.87, 1.00, 0.51]
        green = [0.00, 0.81, 1.00, 0.20, 0.00]
        blue  = [0.81, 1.00, 0.12, 0.00, 0.00]

    stops_ = array('d', stops)
    red_ = array('d', red)
    green_ = array('d', green)
    blue_ = array('d', blue)

    npoints = len(stops_)
    r.TColor.CreateGradientColorTable(npoints, stops_, red_, green_, blue_, ncontours)
    r.gStyle.SetNumberContours(ncontours)

###-------------------------------------------------------------------###

def threeToTwo(h3):
    name = h3.GetName()
    h2 = r.TH2D(name+"_2D", h3.GetTitle(),
                h3.GetNbinsX(), h3.GetXaxis().GetXmin(), h3.GetXaxis().GetXmax(),
                h3.GetNbinsY(), h3.GetYaxis().GetXmin(), h3.GetYaxis().GetXmax(),
                )
    for iX in range(1, 1+h3.GetNbinsX()):
        for iY in range(1, 1+h3.GetNbinsY()):
            content = h3.GetBinContent(iX, iY, 1) + h3.GetBinContent(iX, iY, 2) + h3.GetBinContent(iX, iY, 0)
            h2.SetBinContent(iX, iY, content)
    h2.GetZaxis().SetTitle(h3.GetZaxis().GetTitle())

    #h2.RebinX(2)
    #h2.RebinY(2)

    return h2

###-------------------------------------------------------------------###


def GetHist(File=None, folder=None, hist=None, Norm=None, rebinX=None, rebinY=None):
    h = None

    for f in folder:
        directory = File.Get(f)
        a = directory.Get(hist)
        if h is None:
            h = a.Clone()
        else:
            h.Add(a)

    if rebinX:
        h.RebinX(rebinX)
    if rebinY:
        h.RebinY(rebinY)

    return h

###-------------------------------------------------------------------###


def getRootDirs(bMulti_="", jMulti_="", sitv_=False):

    dirs = []


    for ht in settings["HTBins"]:
        thisHT = ht.split("_")
        if len(thisHT)<2: thisHT.append(None)
        aT = alphaTDict["%s,%s" % (thisHT[0], thisHT[1])][0]
        aTString = "AlphaT%s" % (str(aT[0])+"_"+str(aT[1]) if aT[1] else str(aT[0]))
        dirs.append("smsScan_%s_%s_%s%s_%s" % (bMulti_, jMulti_,
                    "SITV_" if sitv_ else "", aTString, ht))

    return dirs

###-------------------------------------------------------------------###

def getOutFile(model="", htbin="", format="", jMulti_ = "", bMulti_ = "", mode_ = ""):

    if format == "txt":
        outName = "%s_%s_systOutput_%s_%s.txt" % (mode_, model, bMulti_, jMulti_)
    elif "pdf" in format:
        outName = "%s_%s_v%d_%s_%s_%s%s.pdf" % (mode_, model, settings['version'], bMulti_, jMulti_,
                                            htbin, "_dM" if settings["deltaM"] else "")

    return "./out/"+outName

###-------------------------------------------------------------------###


def getPointVal(hist=None, xval=0., yval=0):

    if hist is None:
        return 0

    bin = hist.FindBin(xval, yval)

    return math.fabs(hist.GetBinContent(bin))

###-------------------------------------------------------------------###

def deltaM(h2):

    minVal = 0.;
    maxVal = 0.;

    for iY in range(1, 1+h2.GetNbinsY()):
        for iX in range(1, 1+h2.GetNbinsX()):
            if h2.GetBinContent(iX, iY) > 0.:
                val = h2.GetXaxis().GetBinCenter(iX) - h2.GetYaxis().GetBinCenter(iY)
                if val>maxVal: maxVal=val
                if val<minVal: minVal=val

    nbins = int((int(maxVal)+10 - int(minVal))/h2.GetYaxis().GetBinWidth(5))*2

    h2_dM = r.TH2D(h2.GetName(), h2.GetTitle(),
                    h2.GetNbinsX(), h2.GetXaxis().GetXmin(), h2.GetXaxis().GetXmax(),
                    nbins, int(minVal), int(maxVal)+10)

    for iY in range(1, 1+h2.GetNbinsY()):
        for iX in range(1, 1+h2.GetNbinsX()):
            if h2.GetBinContent(iX, iY) > 0.:
                content = h2.GetBinContent(iX, iY)
                ybinVal = h2.GetXaxis().GetBinLowEdge(iX) - h2.GetYaxis().GetBinLowEdge(iY)
                print ybinVal
                ybin = h2.GetYaxis().FindBin(ybinVal)
                h2_dM.Fill(int(h2.GetXaxis().GetBinCenter(iX)), int(ybinVal), content)

    h2_dM.GetXaxis().SetTitle("mStop (GeV)")
    h2_dM.GetYaxis().SetTitle("deltaM (GeV)")

    #h2_dM.RebinX(2)

    return h2_dM

def get_hist_stat_vals(hist=None):

    min_ = 1.
    max_ = 0.
    avg_ = 0.
    count_ = 0
    
    for i in range(hist.GetNbinsX()*hist.GetNbinsY() + 1000):
        val = abs(hist.GetBinContent(i))

        if val == 1000.: continue

        if val > 0.:
            count_ += 1
            avg_ += val
            if val > max_:
                max_ = val
            if val < min_:
                min_ = val

    avg_ /= float(count_)

    return [min_, max_, avg_]


def copyHist(hist=None, name=""):
    xbins = hist.GetXaxis().GetNbins()
    xlow = hist.GetXaxis().GetBinLowEdge(1)
    xhigh = hist.GetXaxis().GetBinUpEdge(xbins)
    ybins = hist.GetYaxis().GetNbins()
    ylow = hist.GetYaxis().GetBinLowEdge(1)
    yhigh = hist.GetYaxis().GetBinUpEdge(ybins)

    nbins = xbins*ybins

    h = r.TH2D(name, name, xbins, xlow, xhigh, ybins, ylow, yhigh)

    for i in range(nbins+1000):
        val = hist.GetBinContent(i)
        err = hist.GetBinError(i)

        h.SetBinContent(i, val)
        h.SetBinError(i, err)

    return h


###-------------------------------------------------------------------###
###-------------------------------------------------------------------###

def make_syst_map_three(bMulti = "", jMulti = ""):

    set_palette()

    print ">>> Creating systs for", bMulti, jMulti

    xTitle = "m_{stop} (GeV)"

    if settings["deltaM"]:
        yTitle = "deltaM (GeV)"
    else:
        yTitle = "m_{LSP} (GeV)"

    total_syst_hist = r.TH2D()

    my_syst_hists = []

    for n_syst, mode in enumerate(settings["mode"]):

        if "T2cc" not in settings['model'] and mode == "3jet":
            continue

        print "    >>", mode

        if mode == "ISR":
            ins = "isr_"
        else:
            ins = ""

        cut_syst = False # set variable for a before/after type systematic
        if mode in ["LeptonVeto", "MHT_MET", "DeadECAL", "3jet"]:
            cut_syst = True

        # need to find a better way to do this!
        if mode in ["JES", "ISR"]:
            centalRootFile73 = r.TFile.Open("./rootFiles/%s_v%d/sigScan_%s_had_2012_73.7_bt0.0_MChi-1.0.root" %
                                            (settings["model"], settings["version"], settings["model"]))
            centalRootFile86 = r.TFile.Open("./rootFiles/%s_v%d/sigScan_%s_had_2012_86.7_bt0.0_MChi-1.0.root" %
                                            (settings["model"], settings["version"], settings["model"]))
            centalRootFile100 = r.TFile.Open("./rootFiles/%s_v%d/sigScan_%s_had_2012_100.0_bt0.0_MChi-1.0.root" %
                                            (settings["model"], settings["version"], settings["model"]))

            jesPlusRootFile73 = r.TFile.Open("./rootFiles/%s_v%d/%s/sigScan_%s_had_2012_73.7_%s+ve_bt0.0_MChi-1.0.root" %
                                            (settings["model"], settings["version"], mode, settings["model"], ins))
            jesPlusRootFile86 = r.TFile.Open("./rootFiles/%s_v%d/%s/sigScan_%s_had_2012_86.7_%s+ve_bt0.0_MChi-1.0.root" %
                                            (settings["model"], settings["version"], mode, settings["model"], ins))
            jesPlusRootFile100 = r.TFile.Open("./rootFiles/%s_v%d/%s/sigScan_%s_had_2012_100.0_%s+ve_bt0.0_MChi-1.0.root" %
                                             (settings["model"], settings["version"], mode, settings["model"], ins))

            jesNegRootFile73 = r.TFile.Open("./rootFiles/%s_v%d/%s/sigScan_%s_had_2012_73.7_%s-ve_bt0.0_MChi-1.0.root" %
                                            (settings["model"], settings["version"], mode, settings["model"], ins))
            jesNegRootFile86 = r.TFile.Open("./rootFiles/%s_v%d/%s/sigScan_%s_had_2012_86.7_%s-ve_bt0.0_MChi-1.0.root" %
                                            (settings["model"], settings["version"], mode, settings["model"], ins))
            jesNegRootFile100 = r.TFile.Open("./rootFiles/%s_v%d/%s/sigScan_%s_had_2012_100.0_%s-ve_bt0.0_MChi-1.0.root" %
                                            (settings["model"], settings["version"], mode, settings["model"], ins))
        elif cut_syst:
            centalRootFile73 = r.TFile.Open("./rootFiles/%s_v%d/sigScan_%s_had_2012_73.7_bt0.0_MChi-1.0.root" %
                                            (settings["model"], settings["version"], settings["model"]))
            centalRootFile86 = r.TFile.Open("./rootFiles/%s_v%d/sigScan_%s_had_2012_86.7_bt0.0_MChi-1.0.root" %
                                            (settings["model"], settings["version"], settings["model"]))
            centalRootFile100 = r.TFile.Open("./rootFiles/%s_v%d/sigScan_%s_had_2012_100.0_bt0.0_MChi-1.0.root" %
                                            (settings["model"], settings["version"], settings["model"]))

            jesPlusRootFile73 = r.TFile.Open("./rootFiles/%s_v%d/%s/sigScan_%s_had_2012_73.7_bt0.0_MChi-1.0.root" %
                                            (settings["model"], settings["version"], mode, settings["model"]))
            jesPlusRootFile86 = r.TFile.Open("./rootFiles/%s_v%d/%s/sigScan_%s_had_2012_86.7_bt0.0_MChi-1.0.root" %
                                            (settings["model"], settings["version"], mode, settings["model"]))
            jesPlusRootFile100 = r.TFile.Open("./rootFiles/%s_v%d/%s/sigScan_%s_had_2012_100.0_bt0.0_MChi-1.0.root" %
                                            (settings["model"], settings["version"], mode, settings["model"]))         

        else:
            centalRootFile73 = r.TFile.Open("./rootFiles/%s_v%d/sigScan_%s_had_2012_73.7_bt0.0_MChi-1.0.root" %
                                            (settings["model"], settings["version"], settings["model"]))
            centalRootFile86 = r.TFile.Open("./rootFiles/%s_v%d/sigScan_%s_had_2012_86.7_bt0.0_MChi-1.0.root" %
                                            (settings["model"], settings["version"], settings["model"]))
            centalRootFile100 = r.TFile.Open("./rootFiles/%s_v%d/sigScan_%s_had_2012_100.0_bt0.0_MChi-1.0.root" %
                                            (settings["model"], settings["version"], settings["model"]))

            jesPlusRootFile73 = r.TFile.Open("./rootFiles/%s_v%d/%s/sigScan_%s_had_2012_73.7_bt1.0_MChi-1.0.root" %
                                            (settings["model"], settings["version"], mode, settings["model"]))
            jesPlusRootFile86 = r.TFile.Open("./rootFiles/%s_v%d/%s/sigScan_%s_had_2012_86.7_bt1.0_MChi-1.0.root" %
                                            (settings["model"], settings["version"], mode, settings["model"]))
            jesPlusRootFile100 = r.TFile.Open("./rootFiles/%s_v%d/%s/sigScan_%s_had_2012_100.0_bt1.0_MChi-1.0.root" %
                                             (settings["model"], settings["version"], mode, settings["model"]))

            jesNegRootFile73 = r.TFile.Open("./rootFiles/%s_v%d/%s/sigScan_%s_had_2012_73.7_bt-1.0_MChi-1.0.root" %
                                            (settings["model"], settings["version"], mode, settings["model"]))
            jesNegRootFile86 = r.TFile.Open("./rootFiles/%s_v%d/%s/sigScan_%s_had_2012_86.7_bt-1.0_MChi-1.0.root" %
                                            (settings["model"], settings["version"], mode, settings["model"]))
            jesNegRootFile100 = r.TFile.Open("./rootFiles/%s_v%d/%s/sigScan_%s_had_2012_100.0_bt-1.0_MChi-1.0.root" %
                                            (settings["model"], settings["version"], mode, settings["model"]))


        HTList = []
        if settings["inclHT"]:
            HTList = ["incl"]
        else:
            print "fuck."
            exit()

        for htbin in HTList:

            processCrossSections = []
            cuts = []
            cutsJESNeg = []
            cutsJESPlus = []
            nocuts = []

            suf = ""
            if "275_" in htbin:
                suf = "73"
            elif "325_" in htbin:
                suf = "86"
            else:
                suf = "100"

            r.gROOT.SetBatch(r.kTRUE)
            c1 = Print(getOutFile(model=settings["model"], htbin=htbin, format="pdf", bMulti_ = bMulti, jMulti_ = jMulti, mode_ = mode))
            c1.DoPageNum = False

            r.gPad.SetRightMargin(0.15)
            r.gPad.SetLeftMargin(0.15)
            r.gPad.SetTopMargin(0.08)
            r.gPad.SetBottomMargin(0.15)

            if settings["model"] == "T2_4body":
                rebin_y_val = 2
            else:
                rebin_y_val = 1

            nocuts = GetHist(File=centalRootFile100, folder=["smsScan_before", ], hist="m0_m12_mChi_weight", Norm=None, rebinY=rebin_y_val)
            nocuts = threeToTwo(nocuts)

            if settings["inclHT"]:

                cutsHist = GetHist(File=centalRootFile73, folder=getRootDirs(bMulti_ = bMulti, jMulti_ = jMulti, sitv_ = True)[0:2], hist="m0_m12_mChi_weight", Norm=None, rebinY=rebin_y_val).Clone()
                cutsHist.Add(GetHist(File=centalRootFile86, folder=getRootDirs(bMulti_ = bMulti, jMulti_ = jMulti, sitv_ = True)[2:3], hist="m0_m12_mChi_weight", Norm=None, rebinY=rebin_y_val))
                cutsHist.Add(GetHist(File=centalRootFile100, folder=getRootDirs(bMulti_ = bMulti, jMulti_ = jMulti, sitv_ = True)[3:], hist="m0_m12_mChi_weight", Norm=None, rebinY=rebin_y_val))
                cutsHist = threeToTwo(cutsHist)

                cutsJESPlusHist = GetHist(File=jesPlusRootFile73, folder=getRootDirs(bMulti_ = bMulti, jMulti_ = jMulti, sitv_ = True)[0:2], hist="m0_m12_mChi_weight", Norm=None, rebinY=rebin_y_val).Clone()
                cutsJESPlusHist.Add(GetHist(File=jesPlusRootFile86, folder=getRootDirs(bMulti_ = bMulti, jMulti_ = jMulti, sitv_ = True)[2:3], hist="m0_m12_mChi_weight", Norm=None, rebinY=rebin_y_val))
                cutsJESPlusHist.Add(GetHist(File=jesPlusRootFile100, folder=getRootDirs(bMulti_ = bMulti, jMulti_ = jMulti, sitv_ = True)[3:], hist="m0_m12_mChi_weight", Norm=None, rebinY=rebin_y_val))
                cutsJESPlusHist = threeToTwo(cutsJESPlusHist)

                if not cut_syst:
                    cutsJESNegHist = GetHist(File=jesNegRootFile73, folder=getRootDirs(bMulti_ = bMulti, jMulti_ = jMulti, sitv_ = True)[0:2], hist="m0_m12_mChi_weight", Norm=None, rebinY=rebin_y_val).Clone()
                    cutsJESNegHist.Add(GetHist(File=jesNegRootFile86, folder=getRootDirs(bMulti_ = bMulti, jMulti_ = jMulti, sitv_ = True)[2:3], hist="m0_m12_mChi_weight", Norm=None, rebinY=rebin_y_val))
                    cutsJESNegHist.Add(GetHist(File=jesNegRootFile100, folder=getRootDirs(bMulti_ = bMulti, jMulti_ = jMulti, sitv_ = True)[3:], hist="m0_m12_mChi_weight", Norm=None, rebinY=rebin_y_val))
                    cutsJESNegHist = threeToTwo(cutsJESNegHist)
            
            else:
                print "fuck."
                exit()


            if settings["deltaM"]:
                nocuts = deltaM(nocuts)
                cutsHist = deltaM(cutsHist)
                cutsJESPlusHist = deltaM(cutsJESPlusHist)
                if not cut_syst: cutsJESNegHist = deltaM(cutsJESNegHist)

            # zero out the above-diagonal region
            l = [i for i in range(1001)]
            for a in l:
              xbinval = cutsHist.GetXaxis().GetBinCenter(a)
              for b in l:
                ybinval = cutsHist.GetYaxis().GetBinCenter(b)
                if xbinval - ybinval < 0. or a < 0.:
                  bin = cutsHist.FindBin(float(a), float(b))
                  cutsHist.SetBinContent(bin, 0.)
                  cutsJESPlusHist.SetBinContent(bin, 0.)
                  if not cut_syst: cutsJESNegHist.SetBinContent(bin, 0.)

            if mode in ["JES", "ISR"]:
                mini = 0.7
                maxi = 1.3
            elif mode in ["MHT_MET"]:
                mini = 0.65
                maxi = 1.35
            elif mode in ["DeadECAL"]:
                mini = 0.7
                maxi = 1.3
            elif mode in ["LeptonVeto"]:
                mini = 0.65
                maxi = 1.35
            else:
                mini = 0.96
                maxi = 1.04

            c1.canvas.SetLogz()
            offset = 1.1
            # c1.Print()
            c1.canvas.SetLogz(False)
            TotalEff = cutsHist.Clone()
            TotalEff.GetZaxis().SetTitle("Fraction of expected signal yield")
            TotalEff.GetZaxis().SetTitleOffset(offset)
            TotalEff.GetZaxis().SetTitleSize(0.05)
            TotalEff.GetXaxis().SetTitle(xTitle)

            TotalEff.SetTitle("Total Efficiency - Nominal")

            TotalEff.GetYaxis().SetTitleOffset(1.3)
            TotalEff.GetYaxis().SetTitleSize(0.05)
            TotalEff.GetYaxis().SetTitle(yTitle)
            TotalEff.Divide(nocuts)
            maxVal = TotalEff.GetMaximum()

            #TotalEff.Scale(100.)
            if settings["text_plot"] and settings["model"] in ["T2cc", "T2_4body"] :
                r.gStyle.SetPaintTextFormat("0.4f");
                TotalEff.SetMarkerSize(0.8)
                TotalEff.Draw("COLZ TEXT")
            else:
                TotalEff.Draw("COLZ")

            tot = 0.
            ctr = 0
            for i in range(TotalEff.GetNbinsX()*TotalEff.GetNbinsY()):
                val = TotalEff.GetBinContent(i)
                if val > 0.:
                    tot += val
                    ctr += 1

            num0 = r.TLatex(0.17, 0.85, "Average efficiency: %.3f%%" % (float(tot/ctr)*100))
            num0.SetNDC()
            num0.Draw("same")

            c1.Print()

            TotalEffPlus = cutsJESPlusHist.Clone()
            TotalEffPlus.GetZaxis().SetTitle("Relative change in efficiency")
            TotalEffPlus.GetZaxis().SetTitleOffset(offset)
            TotalEffPlus.SetTitle("Total Up Efficiency")
            TotalEffPlus.Divide(nocuts)
            TotalEffPlus.GetXaxis().SetTitle(xTitle)
            TotalEffPlus.GetYaxis().SetTitle(yTitle)
            TotalEffPlus.GetZaxis().SetTitleSize(0.05)

            TotalEffPlus.GetYaxis().SetTitleOffset(1.3)
            TotalEffPlus.GetYaxis().SetTitleSize(0.05)
            TotalEffPlus.SetMaximum(maxVal)
            if settings["text_plot"] and settings["model"] in ["T2cc", "T2_4body"] :
                r.gStyle.SetPaintTextFormat("0.4f");
                TotalEffPlus.SetMarkerSize(0.8)
                TotalEffPlus.Draw("COLZ TEXT")
            else:
                TotalEffPlus.Draw("COLZ")
            c1.Print()
            
            if not cut_syst:
                TotalEffNeg = cutsJESNegHist.Clone()
                TotalEffNeg.GetZaxis().SetTitle("Relative change in efficiency")
                TotalEffNeg.GetZaxis().SetTitleOffset(offset)
                TotalEffNeg.SetTitle("Total Down Efficiency")
                TotalEffNeg.Divide(nocuts)
                TotalEffNeg.GetZaxis().SetTitleSize(0.05)

                TotalEffNeg.GetXaxis().SetTitle(xTitle)
                TotalEffNeg.GetYaxis().SetTitle(yTitle)
                TotalEffNeg.GetYaxis().SetTitleOffset(1.3)
                TotalEffNeg.GetYaxis().SetTitleSize(0.05)
                TotalEffNeg.SetMaximum(maxVal)
                if settings["text_plot"] and settings["model"] in ["T2cc", "T2_4body"] :
                    r.gStyle.SetPaintTextFormat("0.4f");
                    TotalEffNeg.SetMarkerSize(0.8)
                    TotalEffNeg.Draw("COLZ TEXT")
                else:
                    TotalEffNeg.Draw("COLZ")
                c1.Print()

            EffOverJESPlus = TotalEffPlus.Clone()
            if cut_syst:
                EffOverJESPlus.SetTitle("Efficiency Change - %s" % mode)
            else:
                EffOverJESPlus.SetTitle("Up deltaEff (Up / Total)")
            EffOverJESPlus.Divide(TotalEff)
            EffOverJESPlus.Draw("COLZ")
            EffOverJESPlus.GetZaxis().SetRangeUser(0.5, 1.5)
            # c1.Print()

            # invert the above division for cuts
            if mode in ["MHT_MET", "DeadECAL"]:
                for i in range(EffOverJESPlus.GetNbinsX()*EffOverJESPlus.GetNbinsY()+1000):
                    # print EffOverJESPlus.GetBinContent(i)
                    if EffOverJESPlus.GetBinContent(i) > 0:
                        # print 1./EffOverJESPlus.GetBinContent(i), EffOverJESPlus.GetBinContent(i)
                        EffOverJESPlus.SetBinContent(i, 1./EffOverJESPlus.GetBinContent(i))

            if not cut_syst:
                EffOverJESNeg = TotalEffNeg.Clone()
                EffOverJESNeg.SetTitle("Down deltaEff (Down / Total)")
                EffOverJESNeg.Divide(TotalEff)
                EffOverJESNeg.Draw("COLZ")
                EffOverJESNeg.GetZaxis().SetRangeUser(0.5, 1.5)
                # c1.Print()

            r.gStyle.SetOptStat(0)

            if not cut_syst: EffOverJESNegClone = EffOverJESNeg.Clone()
            EffOverJESPlusClone = EffOverJESPlus.Clone()

            #force the zRange accross all bins
            if not cut_syst: 
                for bin in range(EffOverJESNegClone.GetNbinsX()*EffOverJESNegClone.GetNbinsY()+1000):
                  if EffOverJESNegClone.GetBinContent(bin) > 0:
                    if EffOverJESNegClone.GetBinContent(bin) < mini:
                        EffOverJESNegClone.SetBinContent(bin, mini)
                    if EffOverJESNegClone.GetBinContent(bin) > maxi:
                        EffOverJESNegClone.SetBinContent(bin, maxi)
            for bin in range(EffOverJESPlusClone.GetNbinsX()*EffOverJESPlusClone.GetNbinsY()+1000):
              if EffOverJESPlusClone.GetBinContent(bin) > 0:
                if EffOverJESPlusClone.GetBinContent(bin) < mini:
                    EffOverJESPlusClone.SetBinContent(bin, mini)
                if EffOverJESPlusClone.GetBinContent(bin) > maxi:
                    EffOverJESPlusClone.SetBinContent(bin, maxi)

            if mode not in ["MHT_MET", "DeadECAL"]:
                # centre distribution around 0., set all others to 1000               
                for bin in range(EffOverJESPlusClone.GetXaxis().GetNbins()*EffOverJESPlusClone.GetYaxis().GetNbins()+10000):
                  if not cut_syst:
                    if EffOverJESNegClone.GetBinContent(bin) > 0:
                      EffOverJESNegClone.SetBinContent(bin, EffOverJESNegClone.GetBinContent(bin)-1.)
                    else:
                      EffOverJESNegClone.SetBinContent(bin, -1000)
                  if EffOverJESPlusClone.GetBinContent(bin) > 0:
                    EffOverJESPlusClone.SetBinContent(bin, EffOverJESPlusClone.GetBinContent(bin)-1.)
                  else:
                    EffOverJESPlusClone.SetBinContent(bin, -1000)

                #force range
                if not cut_syst: 
                    EffOverJESNegClone.SetMinimum(mini-1.)
                    EffOverJESNegClone.SetMaximum(maxi-1.)
                EffOverJESPlusClone.SetMinimum(mini-1.)
                EffOverJESPlusClone.SetMaximum(maxi-1.)
            else:
                EffOverJESPlusClone.SetTitle("Cut Efficiency - %s" % mode)


            if settings["text_plot"] and settings["model"] in ["T2cc", "T2_4body"] :
                r.gStyle.SetPaintTextFormat("0.4f");
                EffOverJESPlusClone.SetMarkerSize(0.8)
                EffOverJESPlusClone.Draw("COLZ TEXT")
            else:
                EffOverJESPlusClone.Draw("COLZ")
            
            if mode in ["MHT_MET", "DeadECAL"]:
                stat_vals = get_hist_stat_vals(EffOverJESPlusClone)
                stat_vals_string = "Avg=%.3f, Min=%.3f, Max=%.3f"%(stat_vals[2], stat_vals[0], stat_vals[1])
                print "      >", stat_vals_string
                num = r.TLatex(0.16,0.8,stat_vals_string)
                num.SetNDC()
                num.Draw("same")

            c1.Print()

            if not cut_syst: 
                if settings["text_plot"] and settings["model"] in ["T2cc", "T2_4body"] :
                    r.gStyle.SetPaintTextFormat("0.4f");
                    EffOverJESNegClone.SetMarkerSize(0.8)
                    EffOverJESNegClone.Draw("COLZ TEXT")
                else:
                    EffOverJESNegClone.Draw("COLZ")
                c1.Print()

            if mode in ["MHT_MET", "DeadECAL", "3jet"]:
                c1.close()
                continue

            if mode in ["LeptonVeto"]:
                if settings["model"] in ["T2cc", "T2_4body"]:
                    c1.close()
                    continue

            # now calculate the overall systematics, point by point
            syst_hist = EffOverJESPlusClone.Clone()
            for i in range((EffOverJESPlusClone.GetNbinsX() * EffOverJESPlusClone.GetNbinsY()) + 1000):
                
                # skip null points
                if abs(EffOverJESPlusClone.GetBinContent(i)) == 1000.: continue

                if not cut_syst:
                    neg_syst = abs(EffOverJESNegClone.GetBinContent(i))
                pos_syst = abs(EffOverJESPlusClone.GetBinContent(i))

                if cut_syst:
                    this_syst = pos_syst
                else:
                    if pos_syst >= neg_syst:
                        this_syst = pos_syst
                    else:
                        this_syst = neg_syst

                # pick out the largest systematic - most conservative approach
                syst_hist.SetBinContent(i, this_syst)

            syst_hist.SetTitle("Systematics")
            syst_hist.SetMinimum(0.)

            if settings["text_plot"] and settings["model"] in ["T2cc", "T2_4body"]:
                r.gStyle.SetPaintTextFormat("0.4f");
                syst_hist.SetMarkerSize(0.8)
                syst_hist.Draw("COLZ TEXT")
            else:
                syst_hist.Draw("COLZ")

            stat_vals = get_hist_stat_vals(EffOverJESPlusClone)
            stat_vals_string = "Avg=%.3f, Min=%.3f, Max=%.3f"%(stat_vals[2], stat_vals[0], stat_vals[1])
            print "      >", stat_vals_string
            num = r.TLatex(0.16,0.8,stat_vals_string)
            num.SetNDC()
            num.Draw("same")

            c1.Print()

            c1.close()
            my_syst_hists.append(syst_hist.Clone())

    if len(my_syst_hists) == 0:
        return

    if settings["model"] == "T2cc":
        # open Yossof's pdf systematics file
        pdf_file = r.TFile.Open("envCvRelHist.root", 'READ')
        pdf_hist = pdf_file.Get("acc_cvRel_m0_m12")
        my_syst_hists.append(pdf_hist)

    total_syst_hist = r.TH2D("my_th2", "my_th2", my_syst_hists[0].GetNbinsX(), my_syst_hists[0].GetXaxis().GetBinLowEdge(1), my_syst_hists[0].GetXaxis().GetBinUpEdge(my_syst_hists[0].GetNbinsX()),
                                                my_syst_hists[0].GetNbinsY(), my_syst_hists[0].GetYaxis().GetBinLowEdge(1), my_syst_hists[0].GetYaxis().GetBinUpEdge(my_syst_hists[0].GetNbinsY()))

    # sum all point systs in quadrature and add to final histogram
    min_val = 1.
    max_val = 0.
    avg_val = 0.
    count = 0

    # loop over specific mass points and find corresponding bins, as pdf plot is different binning scheme
    stop_masses = []
    split_masses = []

    if settings["model"] == "T2cc":
        stop_masses = [100.+25.*i for i in range(11)]
        split_masses = [5., 10., 20., 30., 40., 60., 80.]
    elif settings["model"] == "T2_4body":
        stop_masses = [100.+25.*i for i in range(12)]
        split_masses = [10.*i for i in range(1,9)]
    else:
        for x_bin in range(1,my_syst_hists[0].GetNbinsX()+1):
            for y_bin in range(1,my_syst_hists[0].GetNbinsY()+1):
                if my_syst_hists[0].GetBinContent(x_bin,y_bin) > 0.:
                    stop_mass = my_syst_hists[0].GetXaxis().GetBinLowEdge(x_bin)
                    if stop_mass not in stop_masses:
                        stop_masses.append(stop_mass)
                    split_mass = my_syst_hists[0].GetXaxis().GetBinLowEdge(x_bin) - my_syst_hists[0].GetYaxis().GetBinLowEdge(y_bin)
                    if split_mass not in split_masses:
                        split_masses.append(split_mass)
        print stop_masses
        print split_masses

    for mstop in stop_masses:
        for mdiff in split_masses:
            mlsp = mstop - mdiff
            val = 0.

            this_bin = my_syst_hists[0].FindBin(mstop, mlsp)

            if my_syst_hists[0].GetBinContent(this_bin) <= 0.:
                print "Zero val bin in syst hist!"
                continue

            for s_hist in my_syst_hists:
                if "cvRel" in s_hist.GetName():
                    this_yossof_bin = s_hist.FindBin(mstop, mlsp)
                    hist_val = abs(s_hist.GetBinContent(this_yossof_bin)-1.)
                else:
                    hist_val = s_hist.GetBinContent(this_bin)

                val += hist_val*hist_val
            val = math.sqrt(val)

            if val > 0:
                count += 1
                avg_val += val
                if val > max_val:
                    max_val = val
                if val < min_val:
                    min_val = val
            
            total_syst_hist.SetBinContent(this_bin, val)

    if min_val == 1.:
        min_val = 0.

    out_name = getOutFile(model=settings["model"], htbin=htbin, format="pdf", bMulti_ = bMulti, jMulti_ = jMulti, mode_ = "total")

    c_total = r.TCanvas()
    r.gPad.SetRightMargin(0.15)
    total_syst_hist.SetTitle("Total Systematic")
    total_syst_hist.SetMaximum(0.35)
    if settings["text_plot"] and settings["model"] in ["T2cc", "T2_4body"] :
        r.gStyle.SetPaintTextFormat("0.4f");
        total_syst_hist.SetMarkerSize(0.8)
        total_syst_hist.Draw("COLZ TEXT")
    else:
        total_syst_hist.Draw("COLZ")

    try:
        avg_val /= float(count)
    except ZeroDivisionError:
        avg_val = 0.

    num = r.TLatex(0.15,0.8,"Avg=%.3f, Min=%.3f, Max=%.3f"%(avg_val, min_val, max_val))
    num.SetNDC()
    num.Draw("same")

    print "\t>> Total: Avg=%.3f, Min=%.3f, Max=%.3f"%(avg_val, min_val, max_val)

    c_total.Print(out_name)

    # write the output rootfile
    out_root_file = r.TFile(out_name.replace(".pdf", ".root"), "RECREATE")
    total_syst_hist.SetName(out_name.replace(".pdf","").replace("_v%d"%settings["version"], "").split("/")[-1])
    total_syst_hist.Write()
    
    out_root_file.Close()
    if settings["model"] == "T2cc":
        pdf_file.Close()

    del out_root_file

if __name__ == "__main__":

    if "list" not in str(type(settings["jMulti"])):
        print "settings['jMulti'] must be list."
        exit()
    if "list" not in str(type(settings["bMulti"])):
        print "settings['bMulti'] must be list."
        exit()

    for jM_ in settings["jMulti"]:
        for bM_ in settings["bMulti"]:
            make_syst_map_three(bMulti = bM_, jMulti = jM_)

