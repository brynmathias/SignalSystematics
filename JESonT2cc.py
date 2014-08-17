#!/usr/bin/env python

import ROOT as r
# from plottingstuff import *
from plottingUtils import Print, MakeCumu, SetBatch
from plotDetails import mapRanges, mapDMRanges, alphaTDict
import math

###-------------------------------------------------------------------###
###-------------------------------------------------------------------###

r.gROOT.SetBatch(r.kTRUE)
r.gStyle.SetOptStat(0)
r.TH1.SetDefaultSumw2(1)
r.TH2.SetDefaultSumw2(1)


settings = {
    "model":    ["T2cc", "T2"][0],
    "version":  20,
    "mode":     ["JES", "ISR", "bTag"][1],
    "inclHT":   [False, True][1],
    "HTBins":   ["200_275", "275_325", "325_375", "375_475", "475_575", "575_675", "675_775", "775_875", "875_975", "975_1075", "1075"],
    "deltaM":   [False, True][0],
    "jMulti":   ["le3j", "ge4j", "eq2j", "eq3j", "ge2j"][-1],
    "bMulti":   ["eq0b", "eq1b", "ge0b"][-1]
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


def getRootDirs(bMulti_="", jMulti_="", sitv_=False):

    dirs = []
    if not bMulti_:
        bMulti = settings["bMulti"]
    else:
        bMulti = bMulti_
        
    if not jMulti_:
        jMulti = settings["jMulti"]
    else:
        jMulti = jMulti_

    ## convert individual string selections to lists for iteration
    if "str" in str(type(bMulti)):
        bMulti = [bMulti]
    if "str" in str(type(jMulti)):
        jMulti = [jMulti]

    for jM in jMulti:
        for bM in bMulti:
            for ht in settings["HTBins"]:
                thisHT = ht.split("_")
                if len(thisHT)<2: thisHT.append(None)
                aT = alphaTDict["%s,%s" % (thisHT[0], thisHT[1])][0]
                aTString = "AlphaT%s" % (str(aT[0])+"_"+str(aT[1]) if aT[1] else str(aT[0]))
                dirs.append("smsScan_%s_%s_%s%s_%s" % (bM, jM,
                            "SITV_" if sitv_ else "", aTString, ht))
    # print dirs
    return dirs

###-------------------------------------------------------------------###


def getOutFile(model="", htbin="", format="", jMulti_ = "", bMulti_ = ""):

    ## convert individual string selections to lists for iteration
    if "str" in str(type(bMulti_)):
        bMulti_ = [bMulti_]
    if "str" in str(type(jMulti_)):
        jMulti_ = [jMulti_]

    if format == "txt":
        outName = "%s_%s_systOutput_%s_%s.txt" % (settings["mode"], model, "_".join(bMulti_), "_".join(jMulti_))
    elif "pdf" in format:
        outName = "%s_%s_%s_%s_%s%s.pdf" % (settings["mode"], model, "_".join(bMulti_), "_".join(jMulti_),
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

###-------------------------------------------------------------------###
###-------------------------------------------------------------------###

def make_syst_map(central_files = [], positive_files = [], negative_files = [], bMulti = "", jMulti = "", syst_mode = ""):
    
    xTitle = "m_{stop} (GeV)"

    if settings["deltaM"]:
        yTitle = "deltaM (GeV)"
    else:
        yTitle = "m_{LSP} (GeV)"

    for mode in settings["mode"]:

        # need to find a better way to do this!
        if mode != "bTag":
            centalRootFile73 = r.TFile.Open("./rootFiles/%s_v%d/sigScan_%s_had_2012_73.7_bt0.0_MChi-1.0.root" %
                                            (settings["model"], settings["version"], settings["model"]))
            centalRootFile86 = r.TFile.Open("./rootFiles/%s_v%d/sigScan_%s_had_2012_86.7_bt0.0_MChi-1.0.root" %
                                            (settings["model"], settings["version"], settings["model"]))
            centalRootFile100 = r.TFile.Open("./rootFiles/%s_v%d/sigScan_%s_had_2012_100.0_bt0.0_MChi-1.0.root" %
                                            (settings["model"], settings["version"], settings["model"]))

            jesPlusRootFile73 = r.TFile.Open("./rootFiles/%s_v%d/%s/sigScan_%s_had_2012_73.7_%s+ve_bt0.0_MChi-1.0.root" %
                                            (settings["model"], settings["version"], settings["mode"], settings["model"], ins))
            jesPlusRootFile86 = r.TFile.Open("./rootFiles/%s_v%d/%s/sigScan_%s_had_2012_86.7_%s+ve_bt0.0_MChi-1.0.root" %
                                            (settings["model"], settings["version"], settings["mode"], settings["model"], ins))
            jesPlusRootFile100 = r.TFile.Open("./rootFiles/%s_v%d/%s/sigScan_%s_had_2012_100.0_%s+ve_bt0.0_MChi-1.0.root" %
                                             (settings["model"], settings["version"], settings["mode"], settings["model"], ins))

            jesNegRootFile73 = r.TFile.Open("./rootFiles/%s_v%d/%s/sigScan_%s_had_2012_73.7_%s-ve_bt0.0_MChi-1.0.root" %
                                            (settings["model"], settings["version"], settings["mode"], settings["model"], ins))
            jesNegRootFile86 = r.TFile.Open("./rootFiles/%s_v%d/%s/sigScan_%s_had_2012_86.7_%s-ve_bt0.0_MChi-1.0.root" %
                                            (settings["model"], settings["version"], settings["mode"], settings["model"], ins))
            jesNegRootFile100 = r.TFile.Open("./rootFiles/%s_v%d/%s/sigScan_%s_had_2012_100.0_%s-ve_bt0.0_MChi-1.0.root" %
                                            (settings["model"], settings["version"], settings["mode"], settings["model"], ins))
        else:
            centalRootFile73 = r.TFile.Open("./rootFiles/%s_v%d/sigScan_%s_had_2012_73.7_bt0.0_MChi-1.0.root" %
                                            (settings["model"], settings["version"], settings["model"]))
            centalRootFile86 = r.TFile.Open("./rootFiles/%s_v%d/sigScan_%s_had_2012_86.7_bt0.0_MChi-1.0.root" %
                                            (settings["model"], settings["version"], settings["model"]))
            centalRootFile100 = r.TFile.Open("./rootFiles/%s_v%d/sigScan_%s_had_2012_100.0_bt0.0_MChi-1.0.root" %
                                            (settings["model"], settings["version"], settings["model"]))

            jesPlusRootFile73 = r.TFile.Open("./rootFiles/%s_v%d/%s/sigScan_%s_had_2012_73.7_bt1.0_MChi-1.0.root" %
                                            (settings["model"], settings["version"], settings["mode"], settings["model"]))
            jesPlusRootFile86 = r.TFile.Open("./rootFiles/%s_v%d/%s/sigScan_%s_had_2012_86.7_bt1.0_MChi-1.0.root" %
                                            (settings["model"], settings["version"], settings["mode"], settings["model"]))
            jesPlusRootFile100 = r.TFile.Open("./rootFiles/%s_v%d/%s/sigScan_%s_had_2012_100.0_bt1.0_MChi-1.0.root" %
                                             (settings["model"], settings["version"], settings["mode"], settings["model"]))

            jesNegRootFile73 = r.TFile.Open("./rootFiles/%s_v%d/%s/sigScan_%s_had_2012_73.7_bt-1.0_MChi-1.0.root" %
                                            (settings["model"], settings["version"], settings["mode"], settings["model"]))
            jesNegRootFile86 = r.TFile.Open("./rootFiles/%s_v%d/%s/sigScan_%s_had_2012_86.7_bt-1.0_MChi-1.0.root" %
                                            (settings["model"], settings["version"], settings["mode"], settings["model"]))
            jesNegRootFile100 = r.TFile.Open("./rootFiles/%s_v%d/%s/sigScan_%s_had_2012_100.0_bt-1.0_MChi-1.0.root" %
                                            (settings["model"], settings["version"], settings["mode"], settings["model"]))

        if mode == "ISR":
            ins = "isr_"
        else:
            ins = ""

        HTList = []
        if settings["inclHT"]:
            HTList = ["incl"]
        else:
            HTList = settings["HTBins"]
            oF = open(getOutFile(model=settings["model"], format="txt", bMulti_ = bMulti, jMulti_ = jMulti), 'w')

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

            c1 = Print(getOutFile(model=settings["model"], htbin=htbin, format="pdf"))
            c1.DoPageNum = False

            r.gPad.SetRightMargin(0.175)
            r.gPad.SetLeftMargin(0.15)
            r.gPad.SetTopMargin(0.08)
            r.gPad.SetBottomMargin(0.15)

            nocuts = GetHist(File=centalRootFile100, folder=["smsScan_before", ], hist="m0_m12_mChi_weight", Norm=None, rebinX=4)
            nocuts = threeToTwo(nocuts)

            if settings["inclHT"]:

                cutsHist = GetHist(File=central_files[0], folder=getRootDirs(bMulti_ = bMulti, jMulti_ = jMulti, sitv_ = True)[0:2], hist="m0_m12_mChi_weight", Norm=None).Clone()
                cutsHist.Add(GetHist(File=central_files[1], folder=getRootDirs(bMulti_ = bMulti, jMulti_ = jMulti, sitv_ = True)[2:3], hist="m0_m12_mChi_weight", Norm=None))
                cutsHist.Add(GetHist(File=central_files[2], folder=getRootDirs(bMulti_ = bMulti, jMulti_ = jMulti, sitv_ = True)[3:], hist="m0_m12_mChi_weight", Norm=None))

                cutsJESPlusHist = GetHist(File=positive_files[0], folder=getRootDirs(bMulti_ = bMulti, jMulti_ = jMulti, sitv_ = True)[0:2], hist="m0_m12_mChi_weight", Norm=None).Clone()
                cutsJESPlusHist.Add(GetHist(File=positive_files[1], folder=getRootDirs(bMulti_ = bMulti, jMulti_ = jMulti, sitv_ = True)[2:3], hist="m0_m12_mChi_weight", Norm=None))
                cutsJESPlusHist.Add(GetHist(File=positive_files[2], folder=getRootDirs(bMulti_ = bMulti, jMulti_ = jMulti, sitv_ = True)[3:], hist="m0_m12_mChi_weight", Norm=None))

                cutsJESNegHist = GetHist(File=negative_files[0], folder=getRootDirs(bMulti_ = bMulti, jMulti_ = jMulti, sitv_ = True)[0:2], hist="m0_m12_mChi_weight", Norm=None).Clone()
                cutsJESNegHist.Add(GetHist(File=negative_files[1], folder=getRootDirs(bMulti_ = bMulti, jMulti_ = jMulti, sitv_ = True)[2:3], hist="m0_m12_mChi_weight", Norm=None))
                cutsJESNegHist.Add(GetHist(File=negative_files[2], folder=getRootDirs(bMulti_ = bMulti, jMulti_ = jMulti, sitv_ = True)[3:], hist="m0_m12_mChi_weight", Norm=None))

            else:
                print "fuck."
                exit()

            # convert to TH2 plots
            cutsHist = threeToTwo(cutsHist)
            cutsJESPlusHist = threeToTwo(cutsJESPlusHist)
            cutsJESNegHist = threeToTwo(cutsJESNegHist)

            if settings["deltaM"]:
                nocuts = deltaM(nocuts)
                cutsHist = deltaM(cutsHist)
                cutsJESPlusHist = deltaM(cutsJESPlusHist)
                cutsJESNegHist = deltaM(cutsJESNegHist)

            # zero out the above-diagonal region
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

            if mode == "JES":
                mini = 0.9
                maxi = 1.1
            else:
                mini = 0.70
                maxi = 1.30

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
            # c1.Print()

            EffOverJESNeg = TotalEffNeg.Clone()
            EffOverJESNeg.SetTitle("Down deltaEff (Down / Total)")
            EffOverJESNeg.Divide(TotalEff)
            EffOverJESNeg.Draw("COLZ")
            EffOverJESNeg.GetZaxis().SetRangeUser(0.5, 1.5)
            # c1.Print()

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
            c1.Print()

            EffOverJESNegClone.Draw("COLZ")
            c1.Print()

            if not settings["inclHT"]:
              oF.write("%s\t\t%.3f +\\- %.3f\n" % (htbin, JesTotClone.GetBinLowEdge(bin68), JesTotClone.GetBinError(bin68)))

            c1.close()

        if not settings["inclHT"]:
          oF.close()


if __name__ == "__main__":

    for jM_ in settings["jMulti"]:
        for bM_ in settings["bMulti"]:
            make_syst_map(bMulti = bM_, jMulti = jM_)

