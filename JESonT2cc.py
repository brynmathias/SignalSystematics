#!/usr/bin/env python

import ROOT as r
from plottingstuff import *
from plottingUtils import Print, MakeCumu
import math

###-------------------------------------------------------------------###
###-------------------------------------------------------------------###

r.gStyle.SetOptStat(0)
r.TH1.SetDefaultSumw2(1)
r.TH2.SetDefaultSumw2(1)

settings = {
    "mode": ["JES", "ISR", "bTag"][2],
    "inclHT": [False, True][1],
    "HTBins": ["275_325", "325_375", "375_475", "475_575", "575_675", "675_775", "775_875", "875"],
    "deltaM": [False, True][0],
    "jMulti": ["le3j", "ge4j", "eq2j", "eq3j"][1],
    "bMulti": ["eq0b", "eq1b"][0]
}

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

    return h

###-------------------------------------------------------------------###


def getRootDirs():

    dirs = []
    bMulti = settings["bMulti"]
    jMulti = settings["jMulti"]

    ## convert individual string selections to lists for iteration
    if "str" in str(type(bMulti)):
        bMulti = [bMulti]
    if "str" in str(type(jMulti)):
        jMulti = [jMulti]

    for jM in jMulti:
        for bM in bMulti:
            for ht in settings["HTBins"]:
                dirs.append("smsScan_%s_%s_AlphaT55_%s" % (bM, jM, ht))

    return dirs

###-------------------------------------------------------------------###


def getOutFile(model="", htbin="", format=""):

    bMulti = settings["bMulti"]
    jMulti = settings["jMulti"]

    ## convert individual string selections to lists for iteration
    if "str" in str(type(bMulti)):
        bMulti = [bMulti]
    if "str" in str(type(jMulti)):
        jMulti = [jMulti]

    if format == "txt":
        outName = "%s_%s_systOutput_%s_%s.txt" % (settings["mode"], model, "_".join(bMulti), "_".join(jMulti))
    elif format == "pdf":
        outName = "%s_%s_%s_%s_%s%s.pdf" % (settings["mode"], model, "_".join(bMulti), "_".join(jMulti),
                                            htbin, "_dM" if settings["deltaM"] else "")

    return "./out/"+outName

###-------------------------------------------------------------------###
###-------------------------------------------------------------------###

model = "T2cc"

xTitle = "m_{stop} (GeV)"

if settings["deltaM"]:
    yTitle = "deltaM (GeV)"
else:
    yTitle = "m_{LSP} (GeV)"

if settings["mode"] == "ISR":
    ins = "isr_"
else:
    ins = ""

# need to find a better way to do this!
if settings["mode"] != "bTag":

    centalRootFile73 = r.TFile.Open("./rootFiles/%s/%s/sigScan_%s_had_2012_73.7_%sbt0.0_MChi-1.0.root" %
                                    (model+"_v5", settings["mode"], model, ins))
    centalRootFile86 = r.TFile.Open("./rootFiles/%s/%s/sigScan_%s_had_2012_86.7_%sbt0.0_MChi-1.0.root" %
                                    (model+"_v5", settings["mode"], model, ins))
    centalRootFile100 = r.TFile.Open("./rootFiles/%s/%s/sigScan_%s_had_2012_100.0_%sbt0.0_MChi-1.0.root" %
                                    (model+"_v5", settings["mode"], model, ins))

    jesPlusRootFile73 = r.TFile.Open("./rootFiles/%s/%s/sigScan_%s_had_2012_73.7_%s+ve_bt0.0_MChi-1.0.root" %
                                    (model+"_v5", settings["mode"], model, ins))
    jesPlusRootFile86 = r.TFile.Open("./rootFiles/%s/%s/sigScan_%s_had_2012_86.7_%s+ve_bt0.0_MChi-1.0.root" %
                                    (model+"_v5", settings["mode"], model, ins))
    jesPlusRootFile100 = r.TFile.Open("./rootFiles/%s/%s/sigScan_%s_had_2012_100.0_%s+ve_bt0.0_MChi-1.0.root" %
                                     (model+"_v5", settings["mode"], model, ins))

    jesNegRootFile73 = r.TFile.Open("./rootFiles/%s/%s/sigScan_%s_had_2012_73.7_%s-ve_bt0.0_MChi-1.0.root" %
                                    (model+"_v5", settings["mode"], model, ins))
    jesNegRootFile86 = r.TFile.Open("./rootFiles/%s/%s/sigScan_%s_had_2012_86.7_%s-ve_bt0.0_MChi-1.0.root" %
                                    (model+"_v5", settings["mode"], model, ins))
    jesNegRootFile100 = r.TFile.Open("./rootFiles/%s/%s/sigScan_%s_had_2012_100.0_%s-ve_bt0.0_MChi-1.0.root" %
                                    (model+"_v5", settings["mode"], model, ins))
else:
    centalRootFile73 = r.TFile.Open("./rootFiles/%s/%s/sigScan_%s_had_2012_73.7_%sbt0.0_MChi-1.0.root" %
                                    (model+"_v5", settings["mode"], model, ins))
    centalRootFile86 = r.TFile.Open("./rootFiles/%s/%s/sigScan_%s_had_2012_86.7_%sbt0.0_MChi-1.0.root" %
                                    (model+"_v5", settings["mode"], model, ins))
    centalRootFile100 = r.TFile.Open("./rootFiles/%s/%s/sigScan_%s_had_2012_100.0_%sbt0.0_MChi-1.0.root" %
                                    (model+"_v5", settings["mode"], model, ins))

    jesPlusRootFile73 = r.TFile.Open("./rootFiles/%s/%s/sigScan_%s_had_2012_73.7_%sbt1.0_MChi-1.0.root" %
                                    (model+"_v5", settings["mode"], model, ins))
    jesPlusRootFile86 = r.TFile.Open("./rootFiles/%s/%s/sigScan_%s_had_2012_86.7_%sbt1.0_MChi-1.0.root" %
                                    (model+"_v5", settings["mode"], model, ins))
    jesPlusRootFile100 = r.TFile.Open("./rootFiles/%s/%s/sigScan_%s_had_2012_100.0_%sbt1.0_MChi-1.0.root" %
                                     (model+"_v5", settings["mode"], model, ins))

    jesNegRootFile73 = r.TFile.Open("./rootFiles/%s/%s/sigScan_%s_had_2012_73.7_%sbt-1.0_MChi-1.0.root" %
                                    (model+"_v5", settings["mode"], model, ins))
    jesNegRootFile86 = r.TFile.Open("./rootFiles/%s/%s/sigScan_%s_had_2012_86.7_%sbt-1.0_MChi-1.0.root" %
                                    (model+"_v5", settings["mode"], model, ins))
    jesNegRootFile100 = r.TFile.Open("./rootFiles/%s/%s/sigScan_%s_had_2012_100.0_%sbt-1.0_MChi-1.0.root" %
                                    (model+"_v5", settings["mode"], model, ins))

HTList = []
if settings["inclHT"]:
    HTList = ["incl"]
else:
    HTList = settings["HTBins"]
    oF = open(getOutFile(model=model, format="txt"), 'w')

for htbin in HTList:

    processCrossSections = []
    cuts = []
    cutsJESNeg = []
    cutsJESPlus = []
    # cutsJESRan = []
    nocuts = []

    suf = ""
    if "275_" in htbin:
        suf = "73"
    elif "325_" in htbin:
        suf = "86"
    else:
        suf = "100"

    c1 = Print(getOutFile(model=model, htbin=htbin, format="pdf"))
    c1.DoPageNum = False

    r.gPad.SetRightMargin(0.175)
    r.gPad.SetLeftMargin(0.15)
    r.gPad.SetTopMargin(0.08)
    r.gPad.SetBottomMargin(0.15)

    nocuts = GetHist(File=centalRootFile100, folder=["smsScan_before", ], hist="m0_m12_mChi_noweight", Norm=None, rebinX=4)
    nocuts = threeToTwo(nocuts)

    if settings["inclHT"]:

        cutsHist = GetHist(File=centalRootFile73, folder=getRootDirs()[0:1], hist="m0_m12_mChi_noweight", Norm=None).Clone()
        cutsHist.Add(GetHist(File=centalRootFile86, folder=getRootDirs()[1:2], hist="m0_m12_mChi_noweight", Norm=None))
        cutsHist.Add(GetHist(File=centalRootFile100, folder=getRootDirs()[2:], hist="m0_m12_mChi_noweight", Norm=None))

        cutsJESPlusHist = GetHist(File=jesPlusRootFile73, folder=getRootDirs()[0:1], hist="m0_m12_mChi_noweight", Norm=None).Clone()
        cutsJESPlusHist.Add(GetHist(File=jesPlusRootFile86, folder=getRootDirs()[1:2], hist="m0_m12_mChi_noweight", Norm=None))
        cutsJESPlusHist.Add(GetHist(File=jesPlusRootFile100, folder=getRootDirs()[2:], hist="m0_m12_mChi_noweight", Norm=None))

        cutsJESNegHist = GetHist(File=jesNegRootFile73, folder=getRootDirs()[0:1], hist="m0_m12_mChi_noweight", Norm=None).Clone()
        cutsJESNegHist.Add(GetHist(File=jesNegRootFile86, folder=getRootDirs()[1:2], hist="m0_m12_mChi_noweight", Norm=None))
        cutsJESNegHist.Add(GetHist(File=jesNegRootFile100, folder=getRootDirs()[2:], hist="m0_m12_mChi_noweight", Norm=None))
    else:
        d = [getRootDirs()[i] for i, x in enumerate(getRootDirs()) if htbin in x]
        cutsHist = GetHist(File=eval("centalRootFile%s" % suf), folder=d, hist="m0_m12_mChi_noweight", Norm=None).Clone()
        cutsJESPlusHist = GetHist(File=eval("jesPlusRootFile%s" % suf), folder=d, hist="m0_m12_mChi_noweight", Norm=None).Clone()
        cutsJESNegHist = GetHist(File=eval("jesNegRootFile%s" % suf), folder=d, hist="m0_m12_mChi_noweight", Norm=None).Clone()

    # convert to TH2 plots
    cutsHist = threeToTwo(cutsHist)
    cutsJESPlusHist = threeToTwo(cutsJESPlusHist)
    cutsJESNegHist = threeToTwo(cutsJESNegHist)

    if settings["deltaM"]:
        nocuts = deltaM(nocuts)
        cutsHist = deltaM(cutsHist)
        cutsJESPlusHist = deltaM(cutsJESPlusHist)
        cutsJESNegHist = deltaM(cutsJESNegHist)

    l = [i for i in range(301)]
    for a in l:
      xbinval = cutsHist.GetXaxis().GetBinCenter(a)
      for b in l:
        ybinval = cutsHist.GetYaxis().GetBinCenter(b)
        if xbinval - ybinval < 0. or a < 0.:
          bin = cutsHist.FindBin(float(a), float(b))
          cutsHist.SetBinContent(bin, 0.)
          cutsJESPlusHist.SetBinContent(bin, 0.)
          cutsJESNegHist.SetBinContent(bin, 0.)

    if settings["mode"] == "JES":
        mini = 0.85
        maxi = 1.15
    else:
        mini = 0.75
        maxi = 1.25

    c1.canvas.SetLogz()
    offset = 1.1
    # c1.Print()
    c1.canvas.SetLogz(False)
    TotalEff = cutsHist.Clone()
    TotalEff.GetZaxis().SetTitle("Fraction of expected signal yield")
    TotalEff.GetZaxis().SetTitleOffset(offset)
    TotalEff.GetZaxis().SetTitleSize(0.05)
    TotalEff.GetXaxis().SetTitle(xTitle)

    TotalEff.SetTitle("Total Efficiency")

    TotalEff.GetYaxis().SetTitleOffset(1.3)
    TotalEff.GetYaxis().SetTitleSize(0.05)
    TotalEff.GetYaxis().SetTitle(yTitle)
    TotalEff.Divide(nocuts)
    maxVal = TotalEff.GetMaximum()

    #TotalEff.Scale(100.)
    #r.gStyle.SetPaintTextFormat("0.2f %%");
    #TotalEff.SetMarkerSize(1.4)
    #r.gStyle.SetPalette(3)
    #TotalEff.Draw("COLZ TEXT40")
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
    TotalEffPlus.Draw("COLZ")
    c1.Print()

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
    TotalEffNeg.Draw("COLZ")
    c1.Print()

    EffOverJESPlus = TotalEffPlus.Clone()
    EffOverJESPlus.SetTitle("Up deltaEff (Up / Total)")
    EffOverJESPlus.Divide(TotalEff)
    EffOverJESPlus.Draw("COLZ")
    EffOverJESPlus.GetZaxis().SetRangeUser(0.5, 1.5)
    #c1.Print()

    EffOverJESNeg = TotalEffNeg.Clone()
    EffOverJESNeg.SetTitle("Down deltaEff (Down / Total)")
    EffOverJESNeg.Divide(TotalEff)
    EffOverJESNeg.Draw("COLZ")
    EffOverJESNeg.GetZaxis().SetRangeUser(0.5, 1.5)
    #c1.Print()

    r.gStyle.SetOptStat(0)

    EffOverJESNegClone = EffOverJESNeg.Clone()
    EffOverJESPlusClone = EffOverJESPlus.Clone()

    #force the zRange accross all bins
    for bin in range(EffOverJESNegClone.GetNbinsX()*EffOverJESNegClone.GetNbinsY()):
      if EffOverJESNegClone.GetBinContent(bin) > 0:
        if EffOverJESNegClone.GetBinContent(bin) < mini:
            EffOverJESNegClone.SetBinContent(bin, mini)
        if EffOverJESNegClone.GetBinContent(bin) > maxi:
            EffOverJESNegClone.SetBinContent(bin, maxi)
    for bin in range(EffOverJESPlusClone.GetNbinsX()*EffOverJESPlusClone.GetNbinsY()):
      if EffOverJESPlusClone.GetBinContent(bin) > 0:
        if EffOverJESPlusClone.GetBinContent(bin) < mini:
            EffOverJESPlusClone.SetBinContent(bin, mini)
        if EffOverJESPlusClone.GetBinContent(bin) > maxi:
            EffOverJESPlusClone.SetBinContent(bin, maxi)

    # centre distribution around 0., set all others to 1000
    for bin in range(EffOverJESPlusClone.GetXaxis().GetNbins()*EffOverJESPlusClone.GetYaxis().GetNbins()+10000):
      if EffOverJESNegClone.GetBinContent(bin) > 0:
        EffOverJESNegClone.SetBinContent(bin, EffOverJESNegClone.GetBinContent(bin)-1.)
      else:
        EffOverJESNegClone.SetBinContent(bin, -1000)
      if EffOverJESPlusClone.GetBinContent(bin) > 0:
        EffOverJESPlusClone.SetBinContent(bin, EffOverJESPlusClone.GetBinContent(bin)-1.)
      else:
        EffOverJESPlusClone.SetBinContent(bin, -1000)

    #force range
    EffOverJESNegClone.SetMinimum(mini-1.)
    EffOverJESNegClone.SetMaximum(maxi-1.)
    EffOverJESPlusClone.SetMinimum(mini-1.)
    EffOverJESPlusClone.SetMaximum(maxi-1.)

    EffOverJESPlusClone.Draw("COLZ")
    #EffOverJESPlusClone.GetZaxis().SetRangeUser(-.14, .14)
    c1.Print()

    EffOverJESNegClone.Draw("COLZ")
    #EffOverJESNegClone.GetZaxis().SetRangeUser(-.14, .14)
    c1.Print()

    oneDJesMinus = r.TH1D("oneDJesMinus", "", 500, 0., 0.5)
    oneDJesPlus = r.TH1D("oneDJesPlus", "", 500, 0., 0.5)
    # oneDJesRan = r.TH1D("oneDJesRan","JES variation for signal efficiency 1D Projection",1000,-5.,5.)
    nEvents = r.TH1D("totEv", "totEv", 2500, 0, 2500)
    minf = 0.9
    maxf = 1.1
    fit = r.TF1("Gaussian", "gaus", minf, maxf)

    ## loop over all bins
    totalBins = TotalEff.GetXaxis().GetNbins()*TotalEff.GetYaxis().GetNbins()
    for bin in range(totalBins):
        ##get deviations from "no change"
        contentMinus = math.fabs(EffOverJESNeg.GetBinContent(bin)-1.)
        errMinus = EffOverJESNeg.GetBinError(bin)
        contentPlus = math.fabs(EffOverJESPlus.GetBinContent(bin)-1.)
        errPlus = EffOverJESPlus.GetBinError(bin)
        #print contentMinus, contentPlus
        # contentRan =   EffOverJESRan.GetBinContent(bin)
        content = nocuts.GetBinContent(bin)
        ## skip if no original events, smaller than .01
        if content == 0:
            continue
        if EffOverJESPlus.GetBinContent(bin) < 0.01:
            continue
        if EffOverJESPlus.GetBinContent(bin) < 0. or EffOverJESPlus.GetBinContent(bin) > 100.:
           nEvents.Fill(math.fabs(cutsJESPlusHist.GetBinContent(bin)-cutsHist.GetBinContent(bin)))
        ## fill the variation in 1d
        if contentMinus > 0.:
            #oneDJesMinus.Fill(math.fabs(contentMinus))
            this = oneDJesMinus.FindBin(contentMinus)
            oneDJesMinus.SetBinContent(this, contentMinus)
            oneDJesMinus.SetBinError(this, errMinus)
        if contentPlus > 0.:
            #oneDJesPlus.Fill(math.fabs(contentPlus))
            this = oneDJesPlus.FindBin(contentPlus)
            oneDJesPlus.SetBinContent(this, contentPlus)
            oneDJesPlus.SetBinError(this, errPlus)

    # with the 1d distro of variations from 1., plot a normalised cumulative distro
    scalValPlus = oneDJesPlus.Integral()
    oneDJesPlus = MakeCumu(oneDJesPlus)
    oneDJesPlus.Scale(1./scalValPlus)
    oneDJesPlus.SetTitle("Up deltaEff")
    oneDJesPlus.Draw("hist")
    oneDJesPlus.GetXaxis().SetTitle("Relative change in efficiency")
    oneDJesPlus.GetXaxis().SetTitleSize(0.05)
    c1.Print()

    r.gStyle.SetOptStat(0)
    scalValMinus = oneDJesMinus.Integral()
    oneDJesMinus = MakeCumu(oneDJesMinus)
    oneDJesMinus.Scale(1./scalValMinus)
    oneDJesMinus.SetTitle("Down deltaEff")
    oneDJesMinus.Draw("hist")
    oneDJesMinus.GetXaxis().SetTitle("Relative change in efficiency")
    oneDJesMinus.GetXaxis().SetTitleSize(0.05)
    c1.Print()

    JesTotal = oneDJesPlus.Clone()
    JesTotal.Add(oneDJesMinus)
    JesTotal.Scale(1./2.)
    JesTotal.SetName("JESTotal")
    bin68 = 0
    for bin in range(JesTotal.GetNbinsX()):
      if JesTotal.GetBinContent(bin) <= 0.68:
        bin68 = bin
    JesTotClone = r.TH1D(JesTotal)
    for bin in range(JesTotal.GetNbinsX()):
      if bin > bin68:
        JesTotClone.SetBinContent(bin, 0.)
        JesTotClone.SetBinError(bin, 0.)
    JesTotal.SetTitle("Total deltaEff (Sum of Up and Down)")
    JesTotal.Draw("h")
    JesTotClone.SetFillColor(r.kRed)
    JesTotClone.Draw("sameh")

    num = r.TLatex(0.4, 0.3, "68%% of events below %.3f" % (JesTotClone.GetBinLowEdge(bin68)))
    num.SetNDC()
    num.Draw("same")

    c1.Print()

    print "\n\nSYST: %f\n\n" % JesTotClone.GetBinLowEdge(bin68)
    if not settings["inclHT"]:
      oF.write("%s\t\t%.3f +\\- %.3f\n" % (htbin, JesTotClone.GetBinLowEdge(bin68), JesTotClone.GetBinError(bin68)))

    c1.close()

if not settings["inclHT"]:
  oF.close()
