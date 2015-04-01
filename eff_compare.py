import ROOT as r
import signalUtils as sutils
import plotting_classes as pcla
from model_versions import versions
from itertools import product
from copy import deepcopy

r.gStyle.SetOptStat(0)
r.gROOT.SetBatch(r.kTRUE)

r.TH1.SetDefaultSumw2(1)
r.TH2.SetDefaultSumw2(1)

settings = {
    "model1":["T2cc", "T1ttcc", "T2tt", "T2bb", "T2_4body", "T2bw_0p25", "T2bw_0p75"][0], # this should be T2cc!
    "version1":35,
    "model2":["T2cc", "T1ttcc", "T2tt", "T2bb", "T2_4body", "T2bw_0p25", "T2bw_0p75"][0],
    "version2":38,
    "HTBins":["200_275", "275_325", "325_375", "375_475", "475_575", "575_675", "675_775", "775_875", "875_975", "975_1075", "1075"],
    "deltaM":[False, True][0],
    "jMulti":["le3j", "ge4j", "ge2j", "eq2j", "eq3j", "eq4j", "ge5j"][:3],
    "bMulti":["eq0b","eq1b","eq2b","eq3b","ge0b"],
    "SITV":[False, True][1],
    "run_mode":["eff_maps", "sitv_acceptance", "leptVeto_acceptance", "btag_compare"][0],
    "combine_output":[False, True][0],
    "text_plot":[False, True][1],
    "sele":["had","muon"][0]
}

def getRootDirs(bMulti_="", jMulti_="", htbins = [], sitv_=False, aT = []):

    dirs = []

    for ht in htbins:
        thisHT = ht.split("_")
        if len(thisHT)<2: thisHT.append(None)
        # aT = pdeets.alphaTDict["%s,%s" % (thisHT[0], thisHT[1])][0]
        aTString = "AlphaT%s" % (str(aT[0])+"_"+str(aT[1]) if aT[1] else str(aT[0]))
        dirs.append("smsScan_%s_%s_%s%s_%s" % (bMulti_, jMulti_,
                    "SITV_" if sitv_ else "", aTString, ht))
    # print dirs
    return dirs

def harvest_effs(eff = None):
    
    effs = {}
    hist = eff._hist.Clone()

    for xbin in range(1, hist.GetNbinsX()+10):
        for ybin in range(1, hist.GetNbinsY()+10):
            content = hist.GetBinContent(xbin, ybin) 
            if content > 0.:
                # print hist.GetXaxis().GetBinCenter(xbin), hist.GetYaxis().GetBinCenter(ybin), content
                x = hist.GetXaxis().GetBinCenter(xbin)
                y = hist.GetYaxis().GetBinCenter(ybin)

                if x not in effs:
                    effs[x] = {}

                effs[x][y] = content

    return effs

def dict_printer(dicto = {}, indent = 1):
  
  print "{ (%d keys)\n" % len(dicto)
  for key in dicto:
    print "\t"*indent, "'%s': " % key,
    if dict == type(dicto[key]):
      dict_printer(dicto[key], indent+1)
    else:
      print dicto[key]
  print "\t"*indent, "}\n"

def get_eff_map(files={}, bMulti="", jMulti="", aT_ = [0.55, None]):

    nocuts = sutils.GetHist(File = files["hi"],folder = ["smsScan_before",],hist = "m0_m12_mChi_weight", Norm = None, rebinY=1)
    nocuts = sutils.threeToTwo(nocuts)

    cutsJESPlusHist = sutils.GetHist(File=files["lo"], folder=sutils.getRootDirs(bMulti_ = bMulti, jMulti_ = jMulti, htbins = settings['HTBins'], sitv_ = settings["SITV"])[0:2], hist="m0_m12_mChi_weight", Norm=None, rebinY=1).Clone()
    cutsJESPlusHist.Add(sutils.GetHist(File=files["mid"], folder=sutils.getRootDirs(bMulti_ = bMulti, jMulti_ = jMulti, htbins = settings['HTBins'], sitv_ = settings["SITV"])[2:3], hist="m0_m12_mChi_weight", Norm=None, rebinY=1))
    cutsJESPlusHist.Add(sutils.GetHist(File=files["hi"], folder=sutils.getRootDirs(bMulti_ = bMulti, jMulti_ = jMulti, htbins = settings['HTBins'], sitv_ = settings["SITV"])[3:], hist="m0_m12_mChi_weight", Norm=None, rebinY=1))
    # cutsJESPlusHist = sutils.GetHist(File=files["hi"], folder=getRootDirs(bMulti_ = bMulti, jMulti_ = jMulti, htbins = settings['HTBins'], sitv_ = settings["SITV"], aT = aT_)[3:], hist="m0_m12_mChi_weight", Norm=None, rebinY=1)
    cutsJESPlusHist = sutils.threeToTwo(cutsJESPlusHist)

    return pcla.effMap(cutsJESPlusHist, nocuts)

def compare_effs(effs1 = {}, effs2 = {}, out_string = ""):
    """do comparison of eff dicts here"""

    c1 = r.TCanvas()
    # c1.Divide(1, 3)
    compare_2d      = r.TH2D("comp_2d", "comp_2d; m_{stop} (GeV); m_{LSP} (GeV);", 13, 62.5, 387.5, 78, -2.5, 387.5)
    # compare_2d      = r.TH2D("comp_2d", "comp_2d; m_{stop} (GeV); m_{LSP} (GeV);", 38, 87.5, 1037.5, 39, -12.5, 962.5)
    compare_1d      = r.TH1D("comp_1d", "comp_1d;Eff change", 100, 0., 2.)
    compare_1d_dm10 = r.TH1D("comp_1d_dm10", "comp_1d_dm10;Eff change", 100, 0.7, 1.3)

    # sutils.set_palette()

    r.gPad.SetRightMargin(0.21)
    r.gPad.SetLeftMargin(0.15)
    r.gPad.SetTopMargin(0.08)
    r.gPad.SetBottomMargin(0.15)

    # deepcopy of effs1, but we'll replace with effs
    compare = deepcopy(effs1)

    key1 = effs1.keys()
    key2 = effs2.keys()
    key1.sort()
    key2.sort()
    #print key1
    #print key2


    for stopmass in compare:
        # check mstop is in model2 plane also
        if stopmass not in effs2:
            continue
        for lspmass in compare[stopmass]:
            # check mlsp is in model2 plane also
            if lspmass not in effs2[stopmass]:
                continue
            compare[stopmass][lspmass] = float(effs1[stopmass][lspmass] / effs2[stopmass][lspmass])
            compare_2d.Fill(float(stopmass), float(lspmass), compare[stopmass][lspmass])
            compare_1d.Fill(compare[stopmass][lspmass])
            if float(stopmass) - float(lspmass) == 10.:
                compare_1d_dm10.Fill(compare[stopmass][lspmass])


    # c1.cd(1)
    p1 = r.TPad("p1", "p1", .02, .36, .98, 0.98)
    p1.SetBorderSize(12)
    p1.SetRightMargin(0.15)
    p1.SetLeftMargin(0.10)
    p1.SetTopMargin(0.08)
    # p1.SetBottomMargin(0.15)
    p1.Draw()
    p2 = r.TPad("p2", "p2", .02, .02, .48, 0.35)
    p2.SetBorderSize(12)
    p2.SetTopMargin(0.08)
    p2.SetBottomMargin(0.08)
    p2.Draw()
    p3 = r.TPad("p3", "p3", .52, .02, .98, 0.35)
    p3.SetBorderSize(12)
    p3.SetTopMargin(0.08)
    p3.SetBottomMargin(0.08)
    p3.Draw()    

    p1.cd()
    r.gStyle.SetPaintTextFormat("0.4f")
    compare_2d.SetMarkerSize(.8)
    compare_2d.SetMarkerColor(r.kWhite)
    compare_2d.Draw("colztext")
    # compare_2d.GetZaxis().SetRangeUser(0.9, 1.25)
    compare_2d.GetZaxis().SetTitle("Eff change")
    compare_2d.SetTitle("%s_v%d / %s_v%d" % (settings['model1'], settings['version1'],
                                                settings['model2'], settings['version2']))
    # c1.Print(out_string.replace("compare", "compare_2d"))

    # c1.cd(2)
    p2.cd()
    compare_1d.Draw("hist")
    compare_1d.SetLabelSize(0.04, "X")
    compare_1d.SetTitleSize(0.04, "X") 
    compare_1d.SetTitleOffset(0.9, "X")
    compare_1d.SetTitle("%s_v%d / %s_v%d - all masses" % (settings['model1'], settings['version1'],
                                                            settings['model2'], settings['version2']))
    # c1.Print(out_string.replace("compare", "compare_1d"))

    # c1.cd(3)
    p3.cd()
    compare_1d_dm10.Draw("hist")
    compare_1d_dm10.SetTitle("%s_v%d / %s_v%d - dm10" % (settings['model1'], settings['version1'],
                                                            settings['model2'], settings['version2']))
    compare_1d_dm10.SetLabelSize(0.04, "X")
    compare_1d_dm10.SetTitleSize(0.04, "X")
    compare_1d_dm10.SetTitleOffset(0.9, "X")
    # c1.Print(out_string.replace("compare", "compare_1d_dm10"))
    # c1.cd()
    c1.Print(out_string)



if __name__ == "__main__":

    print ">>> Running comparison for %s vs %s" % (settings['model1'], settings['model2'])

    infiles_1 = {}
    infiles_1["hi"] = r.TFile.Open("./rootFiles/%s_v%d/sigScan_%s_%s_2012_100.0_bt0.0_MChi-1.0.root"%(settings["model1"], settings['version1'], settings["model1"], settings["sele"]))
    infiles_1["mid"] = r.TFile.Open("./rootFiles/%s_v%d/sigScan_%s_%s_2012_86.7_bt0.0_MChi-1.0.root"%(settings["model1"], settings['version1'], settings["model1"], settings["sele"]))
    infiles_1["lo"] = r.TFile.Open("./rootFiles/%s_v%d/sigScan_%s_%s_2012_73.3_bt0.0_MChi-1.0.root"%(settings["model1"], settings['version1'], settings["model1"], settings["sele"]))
    infiles_2 = {}
    infiles_2["hi"] = r.TFile.Open("./rootFiles/%s_v%d/sigScan_%s_%s_2012_100.0_bt0.0_MChi-1.0.root"%(settings["model2"], settings['version2'], settings["model2"], settings["sele"]))
    infiles_2["mid"] = r.TFile.Open("./rootFiles/%s_v%d/sigScan_%s_%s_2012_86.7_bt0.0_MChi-1.0.root"%(settings["model2"], settings['version2'], settings["model2"], settings["sele"]))
    infiles_2["lo"] = r.TFile.Open("./rootFiles/%s_v%d/sigScan_%s_%s_2012_73.3_bt0.0_MChi-1.0.root"%(settings["model2"], settings['version2'], settings["model2"], settings["sele"]))

    for jM in settings['jMulti']:
        for bM in settings['bMulti']:

            print " >>", bM, jM

            effs_1 = harvest_effs(get_eff_map(infiles_1, bM, jM)) # num
            effs_2 = harvest_effs(get_eff_map(infiles_2, bM, jM)) # denom

            compare_effs(effs_1, effs_2, "out/eff_compare_%s_v%d_vs_%s_v%d_%s_%s.pdf" % (
                            settings['model1'], settings['version1'], 
                            settings['model2'], settings['version2'], bM, jM))
