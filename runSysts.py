#!/usr/bin/env python

import ROOT as r
import signalUtils as sutils
import numpy as np
import math
from copy import deepcopy
from model_versions import versions
from plotDetails import alphaTDict
from plotting_classes import systMap, multiPagePDF

###-------------------------------------------------------------------###
###-------------------------------------------------------------------###

r.gROOT.SetBatch(r.kTRUE)
r.gStyle.SetOptStat(0)
r.TH1.SetDefaultSumw2(1)
r.TH2.SetDefaultSumw2(1)


settings = {
    "model":    ["T2cc", "T2", "T2_4body", "T2tt", "T2bw_0p25", "T2bw_0p75"][2],
    "mode":     ["JES", "ISR", "bTag", "LeptonVeto", "DeadECAL", "MHT_MET", "3jet"][:1],
    "systTests":["JES", "ISR", "bTag", "LeptonVeto", "DeadECAL", "MHT_MET", "3jet", "Lumi", "PDF"][:1],
    "HTBins":   ["200_275", "275_325", "325_375", "375_475", "475_575", "575_675", "675_775", "775_875", "875_975", "975_1075", "1075"][3:],
    "deltaM":   [False, True][0],
    "SITV":     [False, True][1],
    "jMulti":   ["le3j", "ge4j", "ge2j"][:1],
    "bMulti":   ["eq0b", "eq1b", "eq2b", "eq3b", "ge0b"][:1],
    "text_plot":[False, True][1],
    "smooth":   [False, True][0],
    "s_iters":  10,
}


def whiteListOfPoints(respect = True):
    if respect:
        return [250., 240.]
    else:
        return None

# settings['systTests'].pop(3)

# flat systs currently applied flat across plane AND cats
flat_systs = {
    "T2cc":{
        "MHT_MET":  0.02,
        "DeadECAL": 0.02,
        "3jet":     0.04, # estimated from most significant deviation
    },
    "T2_4body":{
        "MHT_MET":  0.02,
        "DeadECAL": 0.02,
        "3jet":     0.04, # copied from T2cc
    },
    "T2bw_0p75":{
        "MHT_MET":  0.02,
        "DeadECAL": 0.02,
    },
    "T2bw_0p25":{
        "MHT_MET":  0.02,
        "DeadECAL": 0.02,
    },
    "T2tt":{
        "MHT_MET":  0.02,
        "DeadECAL": 0.02,
    },
}

# add in flat 4.4% systematic for luminosity
for model in flat_systs:
    flat_systs[model]["Lumi"] = 0.023

###-------------------------------------------------------------------###
def syst_smooth(eff = None, err = None, iterations = 1):
    """smooth the eff"""

    new_hist = eff.Clone()
    n = 0
    # note: ONLY WORKS FOR T2CC AND T24BODY because of binning assumptions
    while n < iterations:
        for xbin in range(1, eff.GetNbinsX()+100):
            for ybin in range(1, eff.GetNbinsY()+100):
                val = eff.GetBinContent(xbin, ybin)
                if val <= 0.: continue
                vals = []
                errs = []
                for xtmp in range(-2,3):
                    tmp_val = eff.GetBinContent( xbin+xtmp, ybin+xtmp*5)
                    if tmp_val > 0.:
                        if err:
                            tmp_err = float(err.GetBinContent( xbin+xtmp, ybin+xtmp*5)/tmp_val)
                        else:
                            tmp_err = 0.
                        vals.append(tmp_val)
                        if tmp_err > 0.:
                            errs.append(float(1./math.pow(tmp_err, 2)))
                        else:
                            errs.append(0.)
                
                # get weighted average considering errs
                err_sum = np.sum(errs)
                if err_sum:
                    for j in range(len(errs)):
                        errs[j]/=err_sum

                    ave_val = np.average(vals, weights=errs)
                else:
                    # should only be in here in err == None
                    print "Non-weighted avg used"
                    ave_val = np.average(vals)

                new_hist.SetBinContent(xbin, ybin, ave_val)
        n+=1

        eff = new_hist.Clone()
    return new_hist

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
        outName = "%s_%s_%s_%s_%s%s.pdf" % (mode_, model, bMulti_, jMulti_,
                                            htbin, "_dM" if settings["deltaM"] else "")

    return "./out/"+outName


###-------------------------------------------------------------------###
###-------------------------------------------------------------------###

def make_syst_map_three(bMulti = "", jMulti = ""):
    '''make systematic studies for a given jet+bjet bin'''

    print ">> Run systs:", bMulti, jMulti

    total_syst_hist = r.TH2D()

    my_syst_hists = []
    my_syst_errhists = []

    # if "LeptonVeto" in settings["systTests"]:
    #     settings["systTests"].pop(3)
    #     # pass

    # loop through all given systematic tests
    for n_syst, mode in enumerate(settings["mode"]):

        # only consider 3jet test for T2cc
        if "T2cc" not in settings['model'] and mode == "3jet":
            print "Warning: 3jet test only available for T2cc. Skipping..."
            continue

        # don't consider LeptonVeto systematic for T2cc
        if settings['model'] == "T2cc" and mode == "LeptonVeto":
            continue

        print "    >>", mode

        # set variable for a before/after type systematic
        if mode in ["LeptonVeto", "MHT_MET", "DeadECAL", "3jet"]:
            cut_syst = True
        else:
            cut_syst = False
        

        ### get input files ###

        if mode == "ISR":
            ins = "isr_"
        else:
            ins = ""

        if mode == "LeptonVeto":
            centalRootFile73 = r.TFile.Open("./rootFiles/%s_v%d/LeptonVeto/sigScan_%s_had_2012_73.3_bt0.0_GRLeptVeto_MChi-1.0.root" %
                                            (settings["model"], versions[settings['model']], settings["model"]))
            centalRootFile86 = r.TFile.Open("./rootFiles/%s_v%d/LeptonVeto/sigScan_%s_had_2012_86.7_bt0.0_GRLeptVeto_MChi-1.0.root" %
                                            (settings["model"], versions[settings['model']], settings["model"]))
            centalRootFile100 = r.TFile.Open("./rootFiles/%s_v%d/LeptonVeto/sigScan_%s_had_2012_100.0_bt0.0_GRLeptVeto_MChi-1.0.root" %
                                            (settings["model"], versions[settings['model']], settings["model"]))
        else:
            centalRootFile73 = r.TFile.Open("./rootFiles/%s_v%d/sigScan_%s_had_2012_73.3_bt0.0_MChi-1.0.root" %
                                            (settings["model"], versions[settings['model']], settings["model"]))
            centalRootFile86 = r.TFile.Open("./rootFiles/%s_v%d/sigScan_%s_had_2012_86.7_bt0.0_MChi-1.0.root" %
                                            (settings["model"], versions[settings['model']], settings["model"]))
            centalRootFile100 = r.TFile.Open("./rootFiles/%s_v%d/sigScan_%s_had_2012_100.0_bt0.0_MChi-1.0.root" %
                                            (settings["model"], versions[settings['model']], settings["model"]))


        # need to find a better way to do this!
        if mode in ["JES", "ISR"]:
            jesPlusRootFile73 = r.TFile.Open("./rootFiles/%s_v%d/%s/sigScan_%s_had_2012_73.3_%s+ve_bt0.0_MChi-1.0.root" %
                                            (settings["model"], versions[settings['model']], mode, settings["model"], ins))
            jesPlusRootFile86 = r.TFile.Open("./rootFiles/%s_v%d/%s/sigScan_%s_had_2012_86.7_%s+ve_bt0.0_MChi-1.0.root" %
                                            (settings["model"], versions[settings['model']], mode, settings["model"], ins))
            jesPlusRootFile100 = r.TFile.Open("./rootFiles/%s_v%d/%s/sigScan_%s_had_2012_100.0_%s+ve_bt0.0_MChi-1.0.root" %
                                             (settings["model"], versions[settings['model']], mode, settings["model"], ins))

            jesNegRootFile73 = r.TFile.Open("./rootFiles/%s_v%d/%s/sigScan_%s_had_2012_73.3_%s-ve_bt0.0_MChi-1.0.root" %
                                            (settings["model"], versions[settings['model']], mode, settings["model"], ins))
            jesNegRootFile86 = r.TFile.Open("./rootFiles/%s_v%d/%s/sigScan_%s_had_2012_86.7_%s-ve_bt0.0_MChi-1.0.root" %
                                            (settings["model"], versions[settings['model']], mode, settings["model"], ins))
            jesNegRootFile100 = r.TFile.Open("./rootFiles/%s_v%d/%s/sigScan_%s_had_2012_100.0_%s-ve_bt0.0_MChi-1.0.root" %
                                            (settings["model"], versions[settings['model']], mode, settings["model"], ins))
        elif cut_syst:
            if mode != "LeptonVeto":
                jesPlusRootFile73 = r.TFile.Open("./rootFiles/%s_v%d/%s/sigScan_%s_had_2012_73.3_bt0.0_MChi-1.0.root" %
                                                (settings["model"], versions[settings['model']], mode, settings["model"]))
                jesPlusRootFile86 = r.TFile.Open("./rootFiles/%s_v%d/%s/sigScan_%s_had_2012_86.7_bt0.0_MChi-1.0.root" %
                                                (settings["model"], versions[settings['model']], mode, settings["model"]))
                jesPlusRootFile100 = r.TFile.Open("./rootFiles/%s_v%d/%s/sigScan_%s_had_2012_100.0_bt0.0_MChi-1.0.root" %
                                                (settings["model"], versions[settings['model']], mode, settings["model"]))         
            else:
                jesPlusRootFile73 = r.TFile.Open("./rootFiles/%s_v%d/%s/sigScan_%s_had_2012_73.3_bt0.0_GLeptVeto_MChi-1.0.root" %
                                                (settings["model"], versions[settings['model']], mode, settings["model"]))
                jesPlusRootFile86 = r.TFile.Open("./rootFiles/%s_v%d/%s/sigScan_%s_had_2012_86.7_bt0.0_GLeptVeto_MChi-1.0.root" %
                                                (settings["model"], versions[settings['model']], mode, settings["model"]))
                jesPlusRootFile100 = r.TFile.Open("./rootFiles/%s_v%d/%s/sigScan_%s_had_2012_100.0_bt0.0_GLeptVeto_MChi-1.0.root" %
                                                (settings["model"], versions[settings['model']], mode, settings["model"]))         
        else:
            jesPlusRootFile73 = r.TFile.Open("./rootFiles/%s_v%d/%s/sigScan_%s_had_2012_73.3_bt1.0_MChi-1.0.root" %
                                            (settings["model"], versions[settings['model']], mode, settings["model"]))
            jesPlusRootFile86 = r.TFile.Open("./rootFiles/%s_v%d/%s/sigScan_%s_had_2012_86.7_bt1.0_MChi-1.0.root" %
                                            (settings["model"], versions[settings['model']], mode, settings["model"]))
            jesPlusRootFile100 = r.TFile.Open("./rootFiles/%s_v%d/%s/sigScan_%s_had_2012_100.0_bt1.0_MChi-1.0.root" %
                                             (settings["model"], versions[settings['model']], mode, settings["model"]))

            jesNegRootFile73 = r.TFile.Open("./rootFiles/%s_v%d/%s/sigScan_%s_had_2012_73.3_bt-1.0_MChi-1.0.root" %
                                            (settings["model"], versions[settings['model']], mode, settings["model"]))
            jesNegRootFile86 = r.TFile.Open("./rootFiles/%s_v%d/%s/sigScan_%s_had_2012_86.7_bt-1.0_MChi-1.0.root" %
                                            (settings["model"], versions[settings['model']], mode, settings["model"]))
            jesNegRootFile100 = r.TFile.Open("./rootFiles/%s_v%d/%s/sigScan_%s_had_2012_100.0_bt-1.0_MChi-1.0.root" %
                                            (settings["model"], versions[settings['model']], mode, settings["model"]))


        ### get histograms ###

        nocuts = sutils.GetHist(File=centalRootFile100, folder=["smsScan_before", ], hist="m0_m12_mChi_weight", Norm=None, rebinY=1)
        nocuts = sutils.threeToTwo(nocuts)

        nocuts_noweight = sutils.GetHist(File=centalRootFile100, folder=["smsScan_before", ], hist="m0_m12_mChi_noweight", Norm=None, rebinY=1)
        nocuts_noweight = sutils.threeToTwo(nocuts_noweight)


        cutsHist = sutils.GetHist(File=centalRootFile73, folder=getRootDirs(bMulti_ = bMulti, jMulti_ = jMulti, sitv_ = settings['SITV'] if mode != "LeptonVeto" else False)[0:2], hist="m0_m12_mChi_weight", Norm=None, rebinY=1).Clone()
        cutsHist.Add(sutils.GetHist(File=centalRootFile86, folder=getRootDirs(bMulti_ = bMulti, jMulti_ = jMulti, sitv_ = settings['SITV'] if mode != "LeptonVeto" else False)[2:3], hist="m0_m12_mChi_weight", Norm=None, rebinY=1))
        cutsHist.Add(sutils.GetHist(File=centalRootFile100, folder=getRootDirs(bMulti_ = bMulti, jMulti_ = jMulti, sitv_ = settings['SITV'] if mode != "LeptonVeto" else False)[3:], hist="m0_m12_mChi_weight", Norm=None, rebinY=1))
        cutsHist = sutils.threeToTwo(cutsHist)

        cutsJESPlusHist = sutils.GetHist(File=jesPlusRootFile73, folder=getRootDirs(bMulti_ = bMulti, jMulti_ = jMulti, sitv_ = settings['SITV'] if mode != "LeptonVeto" else False)[0:2], hist="m0_m12_mChi_weight", Norm=None, rebinY=1).Clone()
        cutsJESPlusHist.Add(sutils.GetHist(File=jesPlusRootFile86, folder=getRootDirs(bMulti_ = bMulti, jMulti_ = jMulti, sitv_ = settings['SITV'] if mode != "LeptonVeto" else False)[2:3], hist="m0_m12_mChi_weight", Norm=None, rebinY=1))
        cutsJESPlusHist.Add(sutils.GetHist(File=jesPlusRootFile100, folder=getRootDirs(bMulti_ = bMulti, jMulti_ = jMulti, sitv_ = settings['SITV'] if mode != "LeptonVeto" else False)[3:], hist="m0_m12_mChi_weight", Norm=None, rebinY=1))
        cutsJESPlusHist = sutils.threeToTwo(cutsJESPlusHist)

        if not cut_syst:
            cutsJESNegHist = sutils.GetHist(File=jesNegRootFile73, folder=getRootDirs(bMulti_ = bMulti, jMulti_ = jMulti, sitv_ = settings['SITV'])[0:2], hist="m0_m12_mChi_weight", Norm=None, rebinY=1).Clone()
            cutsJESNegHist.Add(sutils.GetHist(File=jesNegRootFile86, folder=getRootDirs(bMulti_ = bMulti, jMulti_ = jMulti, sitv_ = settings['SITV'])[2:3], hist="m0_m12_mChi_weight", Norm=None, rebinY=1))
            cutsJESNegHist.Add(sutils.GetHist(File=jesNegRootFile100, folder=getRootDirs(bMulti_ = bMulti, jMulti_ = jMulti, sitv_ = settings['SITV'])[3:], hist="m0_m12_mChi_weight", Norm=None, rebinY=1))
            cutsJESNegHist = sutils.threeToTwo(cutsJESNegHist)

        cutsHist_noweight = sutils.GetHist(File=centalRootFile73, folder=getRootDirs(bMulti_ = bMulti, jMulti_ = jMulti, sitv_ = settings['SITV'] if mode != "LeptonVeto" else False)[0:2], hist="m0_m12_mChi_noweight", Norm=None, rebinY=1).Clone()
        cutsHist_noweight.Add(sutils.GetHist(File=centalRootFile86, folder=getRootDirs(bMulti_ = bMulti, jMulti_ = jMulti, sitv_ = settings['SITV'] if mode != "LeptonVeto" else False)[2:3], hist="m0_m12_mChi_noweight", Norm=None, rebinY=1))
        cutsHist_noweight.Add(sutils.GetHist(File=centalRootFile100, folder=getRootDirs(bMulti_ = bMulti, jMulti_ = jMulti, sitv_ = settings['SITV'] if mode != "LeptonVeto" else False)[3:], hist="m0_m12_mChi_noweight", Norm=None, rebinY=1))
        cutsHist_noweight = sutils.threeToTwo(cutsHist_noweight)

        cutsJESPlusHist_noweight = sutils.GetHist(File=jesPlusRootFile73, folder=getRootDirs(bMulti_ = bMulti, jMulti_ = jMulti, sitv_ = settings['SITV'] if mode != "LeptonVeto" else False)[0:2], hist="m0_m12_mChi_noweight", Norm=None, rebinY=1).Clone()
        cutsJESPlusHist_noweight.Add(sutils.GetHist(File=jesPlusRootFile86, folder=getRootDirs(bMulti_ = bMulti, jMulti_ = jMulti, sitv_ = settings['SITV'] if mode != "LeptonVeto" else False)[2:3], hist="m0_m12_mChi_noweight", Norm=None, rebinY=1))
        cutsJESPlusHist_noweight.Add(sutils.GetHist(File=jesPlusRootFile100, folder=getRootDirs(bMulti_ = bMulti, jMulti_ = jMulti, sitv_ = settings['SITV'] if mode != "LeptonVeto" else False)[3:], hist="m0_m12_mChi_noweight", Norm=None, rebinY=1))
        cutsJESPlusHist_noweight = sutils.threeToTwo(cutsJESPlusHist_noweight)

        if not cut_syst:
            cutsJESNegHist_noweight = sutils.GetHist(File=jesNegRootFile73, folder=getRootDirs(bMulti_ = bMulti, jMulti_ = jMulti, sitv_ = settings['SITV'])[0:2], hist="m0_m12_mChi_noweight", Norm=None, rebinY=1).Clone()
            cutsJESNegHist_noweight.Add(sutils.GetHist(File=jesNegRootFile86, folder=getRootDirs(bMulti_ = bMulti, jMulti_ = jMulti, sitv_ = settings['SITV'])[2:3], hist="m0_m12_mChi_noweight", Norm=None, rebinY=1))
            cutsJESNegHist_noweight.Add(sutils.GetHist(File=jesNegRootFile100, folder=getRootDirs(bMulti_ = bMulti, jMulti_ = jMulti, sitv_ = settings['SITV'])[3:], hist="m0_m12_mChi_noweight", Norm=None, rebinY=1))
            cutsJESNegHist_noweight = sutils.threeToTwo(cutsJESNegHist_noweight)

        if settings["deltaM"]:
            nocuts = sutils.deltaM(nocuts)
            cutsHist = sutils.deltaM(cutsHist)
            cutsJESPlusHist = sutils.deltaM(cutsJESPlusHist)
            nocuts_noweight = sutils.deltaM(nocuts_noweight)
            cutsHist_noweight = sutils.deltaM(cutsHist_noweight)
            cutsJESPlusHist_noweight = sutils.deltaM(cutsJESPlusHist_noweight)
            if not cut_syst:
                cutsJESNegHist = sutils.deltaM(cutsJESNegHist)
                cutsJESNegHist_noweight = sutils.deltaM(cutsJESNegHist_noweight)

        # c5 = r.TCanvas()
        # cutsHist.Draw("colz")
        # # cutsHist.GetXaxis().SetRangeUser(350., 500.)
        # # cutsHist.GetYaxis().SetRangeUser(0., 400.)
        # c5.Print("out/cutHist.pdf")

        # cutsJESPlusHist.Draw("colz")
        # # cutsJESPlusHist.GetXaxis().SetRangeUser(350., 500.)
        # # cutsJESPlusHist.GetYaxis().SetRangeUser(0., 400.)
        # c5.Print("out/cutPosHist.pdf")

        # cutsJESPlusHist.Divide(cutsHist)
        # cutsJESPlusHist.Draw("colz")
        # cutsJESPlusHist.GetZaxis().SetRangeUser(0.95, 1.05)
        # c5.Print("out/outDivide.pdf")

        # create systMap object
        my_systMap = systMap(
            central = cutsHist,
            central_noweight = cutsHist_noweight,
            up = cutsJESPlusHist,
            up_noweight = cutsJESPlusHist_noweight,
            down = cutsJESNegHist if not cut_syst else None,
            down_noweight = cutsJESNegHist_noweight if not cut_syst else None,
            nocuts = nocuts,
            nocuts_noweight = nocuts_noweight,
            test = mode,
            model = settings['model'],
            )

        # make systMap output pdf file
        my_systMap.print_all("%s_%s" % (bMulti, jMulti), plotText = settings['text_plot'])
        

        # if test mode in list, then add to total systematic
        # only add if point-by-point value is used
        if mode in settings["systTests"] and mode in ["JES", "ISR", "bTag", "LeptonVeto"]:
            my_syst_hists.append((deepcopy(my_systMap._syst._hist), deepcopy(my_systMap._syst._errHist)))
            # my_syst_errhists.append(deepcopy(my_systMap._syst._errHist))

        del my_systMap

        centalRootFile100.Close()
        centalRootFile86.Close()
        centalRootFile73.Close()

        jesPlusRootFile100.Close()
        jesPlusRootFile86.Close()
        jesPlusRootFile73.Close()

        if not cut_syst:
            jesNegRootFile100.Close()
            jesNegRootFile86.Close()
            jesNegRootFile73.Close()

    ### Now should have an array of systematic test histograms to combined into a total syst histo ###

    if len(my_syst_hists) == 0:
        return

    if "PDF" in settings['systTests']:
        # open Yossof's pdf systematics file
        pdf_file = r.TFile.Open("%s_systematics_shifted.root" % settings['model'], 'READ')
        pdf_hist = pdf_file.Get("new_pdf")
        my_syst_hists.append(deepcopy(pdf_hist))
        pdf_file.Close()


    # create an empty hist to fill with total systematics
    total_syst_hist = r.TH2D("my_th2", "my_th2", my_syst_hists[0][0].GetNbinsX(), my_syst_hists[0][0].GetXaxis().GetBinLowEdge(1), my_syst_hists[0][0].GetXaxis().GetBinUpEdge(my_syst_hists[0][0].GetNbinsX()),
                                                my_syst_hists[0][0].GetNbinsY(), my_syst_hists[0][0].GetYaxis().GetBinLowEdge(1), my_syst_hists[0][0].GetYaxis().GetBinUpEdge(my_syst_hists[0][0].GetNbinsY()))
    total_syst_errhist = r.TH2D("my_th2_err", "my_th2_err", my_syst_hists[0][0].GetNbinsX(), my_syst_hists[0][0].GetXaxis().GetBinLowEdge(1), my_syst_hists[0][0].GetXaxis().GetBinUpEdge(my_syst_hists[0][0].GetNbinsX()),
                                                my_syst_hists[0][0].GetNbinsY(), my_syst_hists[0][0].GetYaxis().GetBinLowEdge(1), my_syst_hists[0][0].GetYaxis().GetBinUpEdge(my_syst_hists[0][0].GetNbinsY()))

    # sum all point systs in quadrature and add to final histogram
    min_val = 1.
    max_val = 0.
    avg_val = 0.
    count = 0

    print "My syst hists:", my_syst_hists

    for bin in range(1, my_syst_hists[0][0].GetNbinsX() * my_syst_hists[0][0].GetNbinsY() + 1000):
        if my_syst_hists[0][0].GetBinContent(bin) <= 0.:
            continue
        mass_val = 0.
        mass_err = 0.
        this_hist_vals = []
        
        # loop through every systematic
        for s_hist in my_syst_hists:
            # this_s_hist = s_hist
            # this_err_s_hist = s_hist.Clone()
            if 'tuple' in str(type(s_hist)):
                this_s_hist = s_hist[0]
                this_err_s_hist = s_hist[1]
            else:
                this_s_hist = s_hist
            # slightly safer method for pdf systematic
            if "new_pdf" in this_s_hist.GetName() or "acc_cvRel_m0_m12" in this_s_hist.GetName():
                # subtract 1. to get a change around zero
                hist_val = this_s_hist.GetBinContent(bin)
                hist_err = 0.
                if hist_val > 0.:
                    hist_val -= 1.
                    # skip dodgy strip
                    if hist_val < -0.2:
                        hist_val = 0.
                xbin, ybin, zbin = r.Long(0.), r.Long(0.), r.Long(0.)
                this_s_hist.GetBinXYZ(bin, xbin, ybin, zbin)
                if this_s_hist.GetXaxis().GetBinCenter(xbin) != my_syst_hists[0][0].GetXaxis().GetBinCenter(xbin):
                    print "Somethings wrong with pdf hist x binning. 1"
                    

                    print this_s_hist.GetXaxis().GetBinLowEdge(1), this_s_hist.GetXaxis().GetBinUpEdge(this_s_hist.GetNbinsX()), this_s_hist.GetNbinsX()
                    print my_syst_hists[0][0].GetXaxis().GetBinLowEdge(1), my_syst_hists[0][0].GetXaxis().GetBinUpEdge(my_syst_hists[0][0].GetNbinsX()), my_syst_hists[0][0].GetNbinsX()
                    print my_syst_hists[0][0].GetYaxis().GetBinLowEdge(1), my_syst_hists[0][0].GetYaxis().GetBinUpEdge(my_syst_hists[0][0].GetNbinsY()), my_syst_hists[0][0].GetNbinsY()

                    exit()
                if this_s_hist.GetYaxis().GetBinCenter(ybin) != my_syst_hists[0][0].GetYaxis().GetBinCenter(ybin):
                    print "Somethings wrong with pdf hist y binning. 2"
                    print this_s_hist.GetYaxis().GetBinCenter(ybin), my_syst_hists[0][0].GetYaxis().GetBinCenter(ybin)

                    print this_s_hist.GetXaxis().GetBinLowEdge(1), this_s_hist.GetXaxis().GetBinUpEdge(this_s_hist.GetNbinsX()), this_s_hist.GetNbinsX()
                    print my_syst_hists[0][0].GetXaxis().GetBinLowEdge(1), my_syst_hists[0][0].GetXaxis().GetBinUpEdge(my_syst_hists[0][0].GetNbinsX()), my_syst_hists[0][0].GetNbinsX()

                    print this_s_hist.GetYaxis().GetBinLowEdge(1), this_s_hist.GetYaxis().GetBinUpEdge(this_s_hist.GetNbinsY()), this_s_hist.GetNbinsY()
                    print my_syst_hists[0][0].GetYaxis().GetBinLowEdge(1), my_syst_hists[0][0].GetYaxis().GetBinUpEdge(my_syst_hists[0][0].GetNbinsY()), my_syst_hists[0][0].GetNbinsY()
                    exit()

            else:
                hist_val = this_s_hist.GetBinContent(bin)
                hist_err = this_err_s_hist.GetBinContent(bin)
                # print hist_err
                # skip any points that are exactly 666 as this is artificial
                if abs(hist_val) == 666:
                    continue

            mass_val += hist_val*hist_val
            mass_err += 4*hist_err
            # this_hist_vals.append(hist_val)

        # now add the flat contributions from DeadECAL, MHT_MET, "Lumi"
        for flat_test in ["MHT_MET", "DeadECAL", "Lumi", "3jet"]:
            if flat_test not in settings['systTests']:
                continue
            
            try:
                flat_val = flat_systs[settings["model"]][flat_test]
            except KeyError:
                flat_val = 0.

            mass_val += flat_val*flat_val
            mass_err += 0.

        mass_val = math.sqrt(mass_val)
        mass_err = math.sqrt(mass_err)

        # fill total syst histogram with quad summed systematic        
        total_syst_hist.SetBinContent(bin, mass_val)
        total_syst_errhist.SetBinContent(bin, mass_err)

    if settings['smooth']:            
        # apply smoothing algorithm to the final plot
        total_syst_hist = sutils.syst_smooth(eff=total_syst_hist, err=None, iterations=settings['s_iters']).Clone()

    out_name = getOutFile(model=settings["model"], htbin="incl", format="pdf", bMulti_ = bMulti, jMulti_ = jMulti, mode_ = "total")

    syst_vals = []
    for i in range(total_syst_hist.GetNbinsX()*total_syst_hist.GetNbinsY()+1000):
        if total_syst_hist.GetBinContent(i)>0.:
            syst_vals.append(total_syst_hist.GetBinContent(i))

    c_total = r.TCanvas()
    r.gPad.SetRightMargin(0.15)
    total_syst_hist.SetTitle("Total Systematic")
    total_syst_hist.SetMaximum(0.4)
    if settings["text_plot"] and settings["model"] in ["T2cc", "T2_4body"] :
        r.gStyle.SetPaintTextFormat("0.4f");
        total_syst_hist.SetMarkerSize(0.8)
        total_syst_hist.Draw("COLZ TEXT")
    else:
        total_syst_hist.Draw("COLZ")

    # get stat vals
    if len(syst_vals)>0:
        tot_mean = np.mean(syst_vals)
        tot_min = np.min(syst_vals)
        tot_max = np.max(syst_vals)
    else:
        tot_mean = np.nan
        tot_min = np.nan
        tot_max = np.nan

    total_stat_string = "Avg=%.3f, Min=%.3f, Max=%.3f"%(tot_mean, tot_min, tot_max)
    # print "\t>> Total: %s" % total_stat_string
    num = r.TLatex(0.15,0.8, total_stat_string)
    num.SetNDC()
    num.Draw("same")

    c_total.Print(out_name)

    if settings['text_plot']:
        total_syst_errhist.Draw("colztext")    
    else:
        total_syst_errhist.Draw("colz")
    c_total.Print(out_name.replace("total", "error_total"))

    # write the output rootfile
    out_root_file = r.TFile(out_name.replace(".pdf", ".root"), "RECREATE")
    total_syst_hist.SetName(out_name.replace(".pdf","").replace("_v%d"%versions[settings['model']], "").split("/")[-1])
    total_syst_hist.Write()

    total_syst_errhist.Write()
    
    out_root_file.Close()

    del out_root_file
    return total_stat_string

if __name__ == "__main__":

    if "list" not in str(type(settings["jMulti"])):
        print "settings['jMulti'] must be list."
        exit()
    if "list" not in str(type(settings["bMulti"])):
        print "settings['bMulti'] must be list."
        exit()

    total_syst_txt = "*"*32
    total_syst_txt += "\n*** Total systematic summary ***\n"
    total_syst_txt += "*"*32
    total_syst_txt += "\n\n"

    # loop through all bjet and jet mutliplicity combinations
    for jM_ in settings["jMulti"]:
        for bM_ in settings["bMulti"]:
            total_syst_txt += "* %s %s *\n" % (bM_, jM_)
            total_syst_txt += make_syst_map_three(bMulti = bM_, jMulti = jM_)
            total_syst_txt += "\n\n"

    print "\n\n", total_syst_txt
