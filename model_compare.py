#!/usr/bin/env python
import ROOT as r
import signalUtils as sutils
import plotting_classes as pCla
from plotDetails import alphaTDict


settings = {
    "model": ["T2cc", "T2_4body", "T2bw_0p75", "T2bw_0p25", "T2tt"][-1],
    "HTBins":   ["200_275", "275_325", "325_375", "375_475", "475_575", "575_675", "675_775", "775_875", "875_975", "975_1075", "1075"],
    "jMulti":   ["le3j", "ge4j", "ge2j"][:2],
    "bMulti":   ["eq0b", "eq1b", "eq2b", "eq3b", "ge0b"][:2],
}

version_dict = {
    "T2cc":32,
    "T2_4body":7
}

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

def cut_off_end_bin(hist = None):
    new_hist = r.TH2D("new_hist", "new_hist", 13 ,75.0, 400.0, 70, 10.0, 360.0)

    # for bin in range(hist.GetNbinsX()*hist.GetNbinsY()+100):
    #     xbin, ybin, zbin = r.Long(0.), r.Long(0.), r.Long(0.)
    #     hist.GetBinXYZ(bin, xbin, ybin, zbin)
    
    for xbin in range(1, hist.GetNbinsX()+10):
        for ybin in range(1, hist.GetNbinsY()+10):

            val = hist.GetBinContent(xbin, ybin)

            xbinval = hist.GetXaxis().GetBinCenter(xbin)
            ybinval = hist.GetYaxis().GetBinCenter(ybin)

            # print xbinval

            if xbinval in [125., 175., 225., 275., 325.]:
                ybinval -= 5.
                # print xbinval, ybinval
                pass

            new_hist.Fill(xbinval, ybinval, val)
            # new_hist.SetBinContent(xbin, ybin, val)

    return new_hist

def threeToTwo(h3) :
    name = h3.GetName()
    binsz = h3.GetNbinsZ()

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
    for f in folder:
        directory = File.Get(f)
        a = directory.Get(hist)
        if h is None:
            h = a.Clone()
        else: h.Add(a)

    if rebinX:
        h.RebinX(rebinX)
    if rebinY:
        h.RebinY(rebinY)

    return h

def compare_versions(bMulti = None, jMulti = None):
    models = ["T2cc", "T2_4body"]

    print bMulti, jMulti

    eff_maps = {}

    cballs = r.TCanvas()
    cballs.SetGridy(1)

    for v in models:

        centralRootFile100 = r.TFile.Open("./rootFiles/%s_v%s/sigScan_%s_%s_2012_100.0_bt0.0_MChi-1.0.root"%(v, version_dict[v], v, "had"))
        centralRootFile87 = r.TFile.Open("./rootFiles/%s_v%s/sigScan_%s_%s_2012_86.7_bt0.0_MChi-1.0.root"%(v, version_dict[v], v, "had"))
        centralRootFile73 = r.TFile.Open("./rootFiles/%s_v%s/sigScan_%s_%s_2012_73.3_bt0.0_MChi-1.0.root"%(v, version_dict[v], v, "had"))
        nocuts = GetHist(File=centralRootFile100, folder=["smsScan_before", ], hist="m0_m12_mChi_weight", Norm=None, rebinY=1)
        nocuts = threeToTwo(nocuts)

        cutsHist = GetHist(File=centralRootFile73, folder=getRootDirs(bMulti_ = bMulti, jMulti_ = jMulti, sitv_ = True)[0:2], hist="m0_m12_mChi_weight", Norm=None, rebinY=1).Clone()
        cutsHist.Add(GetHist(File=centralRootFile87, folder=getRootDirs(bMulti_ = bMulti, jMulti_ = jMulti, sitv_ = True)[2:3], hist="m0_m12_mChi_weight", Norm=None, rebinY=1))
        cutsHist.Add(GetHist(File=centralRootFile100, folder=getRootDirs(bMulti_ = bMulti, jMulti_ = jMulti, sitv_ = True)[3:], hist="m0_m12_mChi_weight", Norm=None, rebinY=1))
        cutsHist = threeToTwo(cutsHist)

        if v == "T2_4body":
            nocuts = cut_off_end_bin(nocuts)
            cutsHist = cut_off_end_bin(cutsHist)

        if v == "T2cc":
            # nocuts.RebinY(2)
            # cutsHist.RebinY(2)
            pass

        # print v, nocuts.GetNbinsX(), nocuts.GetXaxis().GetBinLowEdge(1), nocuts.GetXaxis().GetBinUpEdge(nocuts.GetNbinsX()), nocuts.GetNbinsY(), nocuts.GetYaxis().GetBinLowEdge(1), nocuts.GetYaxis().GetBinUpEdge(nocuts.GetNbinsY())

        eff_maps[v] = pCla.effMap(cutsHist, nocuts)

        # nocuts.Draw("colztext")
        # cballs.Print("out/nocuts_%s_%s_%s.pdf" % (v, bMulti, jMulti))
        # cutsHist.Draw("colztext")
        # cballs.Print("out/cutsHist_%s_%s_%s.pdf" % (v, bMulti, jMulti))
        # eff_maps[v]._hist.Draw("colztext")
        # cballs.Print("out/effHist_%s_%s_%s.pdf" % (v, bMulti, jMulti))


    compare_map = pCla.effMap(eff_maps[models[0]]._hist, eff_maps[models[1]]._hist)

    c1 = r.TCanvas()

    r.gStyle.SetOptStat(0)

    compare_map._hist.Draw("colztext")
    compare_map._hist.GetZaxis().SetRangeUser(0., 3.)
    compare_map._hist.SetTitle("Eff compare (%s/%s) - %s %s" % (models[0], models[1], bMulti, jMulti))
    c1.Print("out/model_compare_%s_%s_%s.pdf" % ("_vs_".join(models), bMulti, jMulti))

    r.gStyle.SetOptStat("neMRou")
    compare_1d_dm10 = r.TH1D("compare_dm10", "compare dM=10 (%s/%s)- %s %s" % (models[0], models[1], bMulti, jMulti),300, 0., 3.)
    compare_1d_dm20 = r.TH1D("compare_dm20", "compare dM=20 (%s/%s)- %s %s" % (models[0], models[1], bMulti, jMulti),300, 0., 3.)
    
    for i in range(compare_map._hist.GetNbinsX()*compare_map._hist.GetNbinsY()+500):
        xbin, ybin, zbin = r.Long(0.), r.Long(0.), r.Long(0.)
        compare_map._hist.GetBinXYZ(i, xbin, ybin, zbin)
        xbinval = compare_map._hist.GetXaxis().GetBinLowEdge(xbin)
        ybinval = compare_map._hist.GetYaxis().GetBinLowEdge(ybin)

        val = compare_map._hist.GetBinContent(i)
        dmass = xbinval-ybinval
        if val:
            # if xbinval == 250: print xbinval, ybinval, xbinval-ybinval
            if dmass <12:
                # print xbinval, ybinval, dmass
                compare_1d_dm10.Fill(val)
            # if dmass > 10 and dmass <= 20:
            #     print "20split:", xbinval, ybinval, dmass
            #     compare_1d_dm20.Fill(val)



    compare_1d_dm10.Draw("hist")
    r.gStyle.SetOptStat("neMRou")
    c1.Print("out/model_compare_1d_%s_%s_%s.pdf(" % ("_vs_".join(models), bMulti, jMulti)) 

    compare_1d_dm20.Draw("hist")
    c1.Print("out/model_compare_1d_%s_%s_%s.pdf)" % ("_vs_".join(models), bMulti, jMulti))     

if __name__ == "__main__":
    for jM in settings["jMulti"]:
        for bM in settings["bMulti"]:
            compare_versions(bM, jM)