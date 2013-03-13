#!/usr/bin/env python
import sys
import os
import ROOT as r
from plottingstuff import *
from plottingUtils import Print, MakeCumu
import math

r.gStyle.SetOptStat(0)

r.TH1.SetDefaultSumw2(1)
r.TH2.SetDefaultSumw2(1)

settings = {
    "mode":["JES", "ISR"][1],
    "inclHT":[False, True][0],
    "HTBins":["275_325", "325_375", "375_475", "475_575", "575_675", "675_775", "775_875", "875"],
    "deltaM":[False, True][1],
    "jMulti":["le3j", "ge4j", "eq2j", "eq3j"][0],
    "bMulti":["eq0b","eq1b","eq2b"][0]
}


subDirListHigh = [
    "smsScan_%s_%s_AlphaT55_375_475"%(settings["bMulti"], settings["jMulti"]),
    "smsScan_%s_%s_AlphaT55_475_575"%(settings["bMulti"], settings["jMulti"]),
    "smsScan_%s_%s_AlphaT55_575_675"%(settings["bMulti"], settings["jMulti"]),
    "smsScan_%s_%s_AlphaT55_675_775"%(settings["bMulti"], settings["jMulti"]),
    "smsScan_%s_%s_AlphaT55_775_875"%(settings["bMulti"], settings["jMulti"]),
    "smsScan_%s_%s_AlphaT55_875"%(settings["bMulti"], settings["jMulti"]),
]

def threeToTwo(h3) :
    name = h3.GetName()
    binsz = h3.GetNbinsZ()
    # print binsz
    h2 = r.TH2D(name+"_2D",h3.GetTitle(),
                h3.GetNbinsX(), h3.GetXaxis().GetXmin(), h3.GetXaxis().GetXmax(),
                h3.GetNbinsY(), h3.GetYaxis().GetXmin(), h3.GetYaxis().GetXmax(),
                )
                
    for iX in range(1, 1+h3.GetNbinsX()) :
        for iY in range(1, 1+h3.GetNbinsY()) :
            content = h3.GetBinContent(iX, iY, 1) + h3.GetBinContent(iX, iY, 2)+ h3.GetBinContent(iX, iY, 0)
            h2.SetBinContent(iX, iY, content)
    h2.GetZaxis().SetTitle(h3.GetZaxis().GetTitle())
    return h2

def GetHist(File = None, folder = None, hist = None, Norm = None, rebinX = None, rebinY = None):
    h = None
    print File.GetName()
    for f in folder:
        print f
        directory = File.Get(f)
        print hist
        a = directory.Get(hist)
        if h is None:
            h = a.Clone()
        else: h.Add(a)

    return h

def copyHist(hist=None, name=""):
    xbins = hist.GetXaxis().GetNbins()
    xlow = hist.GetXaxis().GetBinLowEdge(1)
    xhigh = hist.GetXaxis().GetBinUpEdge(xbins)
    ybins = hist.GetYaxis().GetNbins()
    ylow = hist.GetYaxis().GetBinLowEdge(1)
    yhigh = hist.GetYaxis().GetBinUpEdge(ybins)

    nbins = xbins*ybins

    h = r.TH2D(name, name, xbins, xlow, xhigh, ybins, ylow, yhigh)

    for i in range(nbins):
        val = hist.GetBinContent(i)
        err = hist.GetBinError(i)

        h.SetBinContent(i, val)
        h.SetBinError(i, err)


    return h

model = "T2cc"
ins = "isr_"
ins = ""

centralRootFile100 = r.TFile.Open("./rootFiles/T2cc_v5/sigScan_%s_had_2012_100.0_%sbt0.0_MChi-1.0.root"%(model, ins))
centralRootFile87 = r.TFile.Open("./rootFiles/T2cc_v5/sigScan_%s_had_2012_86.7_%sbt0.0_MChi-1.0.root"%(model, ins))
centralRootFile73 = r.TFile.Open("./rootFiles/T2cc_v5/sigScan_%s_had_2012_73.7_%sbt0.0_MChi-1.0.root"%(model, ins))

print "./rootFiles/T2cc_v5/sigScan_%s_had_2012_73.7_%sbt0.0_MChi-1.0.root"%(model, ins)

c1 = r.TCanvas()

r.gPad.SetRightMargin(0.23)
r.gPad.SetLeftMargin(0.15)
r.gPad.SetTopMargin(0.08)
r.gPad.SetBottomMargin(0.15)

nocuts = GetHist(File = centralRootFile100,folder = ["smsScan_before",],hist = "m0_m12_mChi_noweight", Norm = None ,rebinX= 1)
nocuts = threeToTwo(nocuts)

cutsJESPlusHist = (GetHist(File = centralRootFile100,folder = subDirListHigh,hist = "m0_m12_mChi_noweight", Norm = None ,rebinX= 1))
cutsJESPlusHist.Add(GetHist(File = centralRootFile87,folder = ["smsScan_%s_%s_AlphaT55_325_375"%(settings["bMulti"], settings["jMulti"])],hist = "m0_m12_mChi_noweight", Norm = None ,rebinX= 1))
cutsJESPlusHist.Add(GetHist(File = centralRootFile73,folder = ["smsScan_%s_%s_AlphaT55_275_325"%(settings["bMulti"], settings["jMulti"])],hist = "m0_m12_mChi_noweight", Norm = None ,rebinX= 1))

### RUN FOR
# eq2j
# eq3j
# le3j
# ge4j


cutsJESPlusHist = threeToTwo(cutsJESPlusHist)

eff = cutsJESPlusHist.Clone()
eff.Divide(nocuts)
eff.RebinY(2)

eff.GetXaxis().SetTitle("mStop (GeV)")
eff.GetYaxis().SetTitle("mLSP (GeV)")
eff.GetZaxis().SetTitle("Efficiency")
eff.GetZaxis().SetRangeUser(0.,0.005)

eff.SetTitleSize(0.05,"x")
eff.SetTitleOffset(1.2,"x")
eff.SetTitleSize(0.05,"y")
eff.SetTitleOffset(1.2,"y")
eff.SetTitleSize(0.05,"z")
eff.SetTitleOffset(1.6,"z")

eff.SetTitle("Total Efficiency %s %s"%(settings["bMulti"], settings["jMulti"]))

eff.Draw("COLZ")

#nocuts_copy = copyHist(nocuts)
#cutsJESPlusHist_copy = copyHist(cutsJESPlusHist)
#
#eff = r.TEfficiency(cutsJESPlusHist_copy, nocuts_copy)

c1.Print("EffMap_T2cc_%s_%s.pdf"%(settings["bMulti"], settings["jMulti"]))



