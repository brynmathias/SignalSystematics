#!/usr/bin/env python
import ROOT as r
import signalUtils as sutils
import plotting_classes as pCla
from plotDetails import alphaTDict


settings = {
	"model": ["T2cc", "T2_4body", "T2bw_0p75", "T2bw_0p25", "T2tt"][1],
	"HTBins":   ["200_275", "275_325", "325_375", "375_475", "475_575", "575_675", "675_775", "775_875", "875_975", "975_1075", "1075"],
    "compare": ["eff", "syst"][0],
    "jMulti":   ["le3j", "ge4j", "ge2j"][:2],
    "bMulti":   ["eq0b", "eq1b", "eq2b", "eq3b", "ge0b"][:2],
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

def compare_versions(bMulti = None, jMulti = None):
	
	# comparison is: [num, denom]
	versions = ["15", "9"]

	print bMulti, jMulti

	eff_maps = {}

	for v in versions:

		if settings['compare'] == "eff":
			centralRootFile100 = r.TFile.Open("./rootFiles/%s_v%s/sigScan_%s_%s_2012_100.0_bt0.0_MChi-1.0.root"%(settings["model"], v, settings["model"], "had"))
			centralRootFile87 = r.TFile.Open("./rootFiles/%s_v%s/sigScan_%s_%s_2012_86.7_bt0.0_MChi-1.0.root"%(settings["model"], v, settings["model"], "had"))
			# if v=="9":
			# 	centralRootFile73 = r.TFile.Open("./rootFiles/%s_v%s/sigScan_%s_%s_2012_73.7_bt0.0_MChi-1.0.root"%(settings["model"], v, settings["model"], "had"))
			# elif v=="15":
			# 	centralRootFile73 = r.TFile.Open("./rootFiles/%s_v%s/sigScan_%s_%s_2012_73.3_bt0.0_MChi-1.0.root"%(settings["model"], v, settings["model"], "had"))
			centralRootFile73 = r.TFile.Open("./rootFiles/%s_v%s/sigScan_%s_%s_2012_73.3_bt0.0_MChi-1.0.root"%(settings["model"], v, settings["model"], "had"))
			nocuts = sutils.GetHist(File=centralRootFile100, folder=["smsScan_before", ], hist="m0_m12_mChi_weight", Norm=None, rebinY=1)
			nocuts = sutils.threeToTwo(nocuts)

			cutsHist = sutils.GetHist(File=centralRootFile73, folder=getRootDirs(bMulti_ = bMulti, jMulti_ = jMulti, sitv_ = True)[0:2], hist="m0_m12_mChi_weight", Norm=None, rebinY=1).Clone()
			cutsHist.Add(sutils.GetHist(File=centralRootFile87, folder=getRootDirs(bMulti_ = bMulti, jMulti_ = jMulti, sitv_ = True)[2:3], hist="m0_m12_mChi_weight", Norm=None, rebinY=1))
			cutsHist.Add(sutils.GetHist(File=centralRootFile100, folder=getRootDirs(bMulti_ = bMulti, jMulti_ = jMulti, sitv_ = True)[3:], hist="m0_m12_mChi_weight", Norm=None, rebinY=1))
			cutsHist = sutils.threeToTwo(cutsHist)
			eff_maps[v] = pCla.effMap(cutsHist, nocuts)

		elif settings['compare'] == "syst":
			rootfile = r.TFile.Open("./rootfiles/%s_v%s/syst.root" % (settings["model"], v))
			syst_hist = rootfile.Get("total_T2cc_%s_%s_incl" % (bMulti, jMulti))
			eff_maps[v] = pCla.effMap(syst_hist)

	compare_map = pCla.effMap(eff_maps[versions[0]]._hist, eff_maps[versions[1]]._hist)

	c1 = r.TCanvas()

	r.gStyle.SetOptStat(0)

	compare_map._hist.Draw("colztext")
	compare_map._hist.GetZaxis().SetRangeUser(0.7, 1.3)
	compare_map._hist.SetTitle("%s_%s" % (bMulti, jMulti))
	c1.Print("out/version_compare_%s_%s_%s_%s.pdf" % (settings["model"], "_vs_".join(versions), bMulti, jMulti))

	r.gStyle.SetOptStat("neMRou")
	compare_1d = r.TH1D("compare", "compare",80, 0.6, 1.4)
	for i in range(compare_map._hist.GetNbinsX()*compare_map._hist.GetNbinsY()+500):
		val = compare_map._hist.GetBinContent(i)
		if val: compare_1d.Fill(val)
	compare_1d.Draw("hist")
	r.gStyle.SetOptStat("neMRou")
	compare_1d.SetTitle("%s_%s" % (bMulti, jMulti))
	c1.Print("out/version_compare_1d_%s_%s_%s_%s.pdf" % (settings["model"], "_vs_".join(versions), bMulti, jMulti))	

if __name__ == "__main__":
	for jM in settings["jMulti"]:
		for bM in settings["bMulti"]:
			compare_versions(bM, jM)
