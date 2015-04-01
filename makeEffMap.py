#!/usr/bin/env python

import itertools
import ROOT as r
from array import array
from model_versions import versions
from plotDetails import mapRanges, mapDMRanges, alphaTDict

r.gStyle.SetOptStat(0)
r.gROOT.SetBatch(r.kTRUE)

r.TH1.SetDefaultSumw2(1)
r.TH2.SetDefaultSumw2(1)

settings = {
    "model":["T2cc", "T1ttcc", "T2tt", "T2bb", "T2_4body", "T2bw_0p25", "T2bw_0p75", "pmssm"][0],
    "HTBins":["200_275", "275_325", "325_375", "375_475", "475_575", "575_675", "675_775", "775_875", "875_975", "975_1075", "1075"],
    "deltaM":[False, True][0],
    "jMulti":["le3j", "ge4j", "ge2j", "eq2j", "eq3j", "eq4j", "ge5j"][:2],
    "bMulti":["eq0b","eq1b","eq2b","eq3b","ge0b"][:-1],
    "SITV":[False, True][1],
    "run_mode":["eff_maps", "sitv_acceptance", "leptVeto_acceptance", "btag_compare", "minbias_acceptance"][0],
    "combine_output":[False, True][0],
    "text_plot":[False, True][1],
    "sele":["had","muon"][1]
}

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
                if settings['sele'] == "had":
                    aTString = "AlphaT%s" % (str(aT[0])+"_"+str(aT[1]) if aT[1] else str(aT[0]))

                    out_string = "smsScan_%s_%s_%s%s_%s" % (bM, jM,
                            "SITV_" if sitv_ else "", aTString, ht)
                else:
                    aTString = "NoAlphaT"
                    # out_string = "smsScan_%s_%s_%s%s_%s" % (bM, jM,
                    #         aTString, "_SITV" if sitv_ else "", ht)
                    out_string = "smsScan_%s_%s%s_%s_%s" % (bM, jM,
                            "_SITV" if sitv_ else "", aTString, ht)

                dirs.append(out_string)

    return dirs

def getOutFile(model="", htbin="", format="", bMulti_="", jMulti_="", sitv_=False, prefix=""):

    if not bMulti_:
        bMulti = settings["bMulti"]
    else:
        bMulti = bMulti_
        
    if not jMulti_:
        jMulti = settings["jMulti"]
    else:
        jMulti = jMulti_

    ## convert individual string selections to lists for iteration
    if "list" not in str( type(bMulti) ):
        bMulti = [bMulti]
    if "list" not in str( type(jMulti) ):
        jMulti = [jMulti]

    if format == "txt":
        outName = "%s_v%d_%s_systOutput_%s_%s.txt" % (model, versions[settings['model']], settings["run_mode"], "_".join(bMulti), "_".join(jMulti))
    elif "pdf" in format:
        outName = "%s_%s_%s_%s_%s%s%s.pdf" % (model, settings['sele'], settings['run_mode'], "_".join(bMulti), "_".join(jMulti),
                                            "_SITV" if sitv_ else "",
                                            "_dM" if settings["deltaM"] else "")

    if prefix:
        outName = prefix + "_" + outName

    return "./out/"+outName

def getScanMax(hist=None):
    nbins = hist.GetNbinsX()*hist.GetNbinsY()

    maxCont = 0.
    for i in range(nbins+1000):
        val = hist.GetBinContent(i)
        if val>maxCont:
            maxCont = val
    
    # return(0.001)
    # return max val rounded to 4 sig figs  
    return round(maxCont,4)

def getScanMin(hist=None):
    nbins = hist.GetNbinsX()*hist.GetNbinsY()

    minCont = 1.
    for i in range(nbins+1000):
        val = hist.GetBinContent(i)
        if val == 0: continue
        if val<minCont:
            minCont = val 

    return round(minCont,3)

def deltaM(h2):

    minVal = 100.;
    maxVal = 0.;

    # get the splitting range
    for iY in range(1, 1+h2.GetNbinsY()):
        for iX in range(1, 1+h2.GetNbinsX()):
            if h2.GetBinContent(iX, iY) > 0.:
                val = h2.GetXaxis().GetBinCenter(iX) - h2.GetYaxis().GetBinCenter(iY)
                if val>maxVal:
                    maxVal=val
                if val<minVal:
                    minVal=val

    nbins = int((int(maxVal)+10. - int(minVal))/h2.GetYaxis().GetBinWidth(5))

    h2_dM = r.TH2D(h2.GetName(), h2.GetTitle(),
                    h2.GetNbinsX(), h2.GetXaxis().GetXmin(), h2.GetXaxis().GetXmax(),
                    nbins, int(minVal), int(maxVal)+10.)

    for iY in range(1, 1+h2.GetNbinsY()):
        for iX in range(1, 1+h2.GetNbinsX()):
            if h2.GetBinContent(iX, iY) > 0.:
                content = h2.GetBinContent(iX, iY)
                ybinVal = h2.GetXaxis().GetBinCenter(iX) - h2.GetYaxis().GetBinCenter(iY)
                ybin = h2.GetYaxis().FindBin(ybinVal)
                h2_dM.Fill(int(h2.GetXaxis().GetBinCenter(iX)), int(ybinVal), content)

    h2_dM.GetXaxis().SetTitle("mStop (GeV)")
    h2_dM.GetYaxis().SetTitle("deltaM (GeV)")

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

    try:
        avg_ /= float(count_)
    except ZeroDivisionError:
        avg_ = 0.

    return [min_, max_, avg_]


def make_eff_map_plot(file73 = None, file87 = None, file100 = None, bMulti = "", jMulti = "", sitv=False, no_cuts=False, file_output=False):

    print ">>> Making eff map for", jMulti, bMulti

    c1 = r.TCanvas()

    set_palette()

    r.gPad.SetRightMargin(0.21)
    r.gPad.SetLeftMargin(0.15)
    r.gPad.SetTopMargin(0.08)
    r.gPad.SetBottomMargin(0.15)

    if settings["deltaM"]:
        xRange = mapDMRanges[settings["model"]][0]
        yRange = mapDMRanges[settings["model"]][1]
        suf = "_dM"
    else:
        xRange = mapRanges[settings["model"]][0]
        yRange = mapRanges[settings["model"]][1]

    rebin_y_val = 1

    nocuts = GetHist(File = file100,folder = ["smsScan_before",],hist = "m0_m12_mChi_weight", Norm = None, rebinY=rebin_y_val)
    nocuts = threeToTwo(nocuts)

    nocuts.GetXaxis().SetRangeUser(xRange[0], xRange[1])
    nocuts.GetYaxis().SetRangeUser(yRange[0], yRange[1])

    if settings["text_plot"]:
        r.gStyle.SetPaintTextFormat("0.4f");
        nocuts.SetMarkerSize(0.8)
        nocuts.Draw("COLZ TEXT")
    else:
        nocuts.Draw("COLZ")

    if no_cuts: return nocuts

    cutsJESPlusHist = GetHist(File=file73, folder=getRootDirs(bMulti_ = bMulti, jMulti_ = jMulti, sitv_ = sitv)[0:2], hist="m0_m12_mChi_weight", Norm=None, rebinY=rebin_y_val).Clone()
    cutsJESPlusHist.Add(GetHist(File=file87, folder=getRootDirs(bMulti_ = bMulti, jMulti_ = jMulti, sitv_ = sitv)[2:3], hist="m0_m12_mChi_weight", Norm=None, rebinY=rebin_y_val))
    cutsJESPlusHist.Add(GetHist(File=file100, folder=getRootDirs(bMulti_ = bMulti, jMulti_ = jMulti, sitv_ = sitv)[3:], hist="m0_m12_mChi_weight", Norm=None, rebinY=rebin_y_val))
    cutsJESPlusHist = threeToTwo(cutsJESPlusHist)

    # cutsJESPlusHist.Draw("COLZTEXT")
    # c1.Print("test.pdf")  

    if settings["deltaM"]:
        nocuts = deltaM(nocuts).Clone()
        cutsJESPlusHist = deltaM(cutsJESPlusHist).Clone()

    if file_output:
        c1.Print("out/%s_v%d_nocutsEff.pdf" % (settings['model'], versions[settings['model']]))

    eff = cutsJESPlusHist.Clone()
    eff.Divide(nocuts)

    eff.GetXaxis().SetTitle("mStop (GeV)")

    if settings["deltaM"]:
        eff.GetYaxis().SetTitle("deltaM (GeV)")
    else:
        eff.GetYaxis().SetTitle("mLSP (GeV)")

    eff.GetXaxis().SetRangeUser(xRange[0], xRange[1])
    eff.GetYaxis().SetRangeUser(yRange[0], yRange[1])     
    eff.GetZaxis().SetTitle("Fraction of expected signal yield")
    if settings['sele'] == "muon":
        pass
    else:
        eff.GetZaxis().SetRangeUser(0.,getScanMax(eff))

    eff.SetTitleSize(0.05,"x")
    eff.SetTitleOffset(1.2,"x")
    eff.SetTitleSize(0.05,"y")
    eff.SetTitleOffset(1.2,"y")
    eff.SetTitleSize(0.05,"z")
    eff.SetTitleOffset(1.6,"z")

    eff.SetLabelSize(0.04, "z")

    eff.SetTitle("Total Efficiency - %s %s"%(bMulti, jMulti))

    if settings["text_plot"]:
        r.gStyle.SetPaintTextFormat("0.4f");
        eff.SetMarkerSize(0.8)
        eff.Draw("COLZ TEXT")
    else:
        eff.Draw("COLZ")

    if not settings["deltaM"]:
        diag = r.TF1("diag", "x", 0., 2000.)
        diag.Draw("same")
        diag.SetLineStyle(2)

    stat_vals = get_hist_stat_vals(eff)
    num = r.TLatex(0.2,0.85,"#scale[0.6]{Avg=%.3f,Min=%.3f,Max=%.3f}"%(stat_vals[2], stat_vals[0], stat_vals[1]))
    num.SetNDC()
    num.Draw("same")

    if file_output:
        c1.Print( getOutFile(model=settings["model"], format="pdf", bMulti_=bMulti, jMulti_=jMulti, sitv_=sitv) )
    else:
        return eff

def make_var_eff_map_plot(raw_eff = None, variation_eff = None, bMulti = "", jMulti = "", file_output = False, test_name = ""):

    c1 = r.TCanvas()

    sitv_ratio = raw_eff.Clone()
    sitv_ratio.Divide(variation_eff)

    sitv_ratio.SetTitle("%s - %s %s" % (test_name, bMulti, jMulti))

    if settings["text_plot"]:
        r.gStyle.SetPaintTextFormat("0.4f");
        sitv_ratio.SetMarkerSize(0.8)
        draw_string = "colztext"
    else:
        draw_string = "colz"
    sitv_ratio.Draw(draw_string)

    if settings['model'] in ['T2cc', 'T2_4body']:
        sitv_ratio.GetZaxis().SetRangeUser( 0.6, 1. )
    else:
        sitv_ratio.GetZaxis().SetRangeUser( 0.6, 1. )

    if file_output:
        c1.Print( getOutFile(model=settings["model"], format="pdf", bMulti_=bMulti, jMulti_=jMulti, prefix=test_name) )
    else:
        return sitv_ratio


def print_plot_list(list = [], combine = True, add_string = ""):

    if len(list)<1:
        print ">>> Error in print_plot_list."
        print "  > Plot list empty."
        return

    c1 = r.TCanvas()

    for n, plot in enumerate(list):
        if n==0: suf = "("
        elif n==len(list)-1: suf = ")"
        else: suf = ""

        if settings["text_plot"]:
            r.gStyle.SetPaintTextFormat("0.4f");
            plot.SetMarkerSize(0.8)
            plot.Draw("COLZ TEXT")
        else:
            plot.Draw("COLZ")

        plot.SetTitleSize(0.05,"x")
        plot.SetTitleOffset(1.2,"x")
        plot.SetTitleSize(0.05,"y")
        plot.SetTitleOffset(1.2,"y")
        plot.SetTitleSize(0.05,"z")
        plot.SetTitleOffset(1.6,"z")

        plot.SetLabelSize(0.04, "z")

        c1.Print( "out/%s_%s%s.pdf%s" % (settings['model'], settings['run_mode'],
                                                "_" if add_string else "", suf) )

def divide_eff_maps(num = None, denom = None, title = "", file_output = False):

    c1 = r.TCanvas()

    eff_change = num.Clone()
    eff_change.Divide(denom)

    eff_change.SetTitle("Btag Compare - " + title)

    eff_change.GetZaxis().SetRangeUser(0., getScanMax(eff_change))

    eff_change.Draw("colz")

    if file_output:
        print "Chris - enable individual file printing for btag compare!"
    else:
        return eff_change

if __name__ == "__main__":
    model = settings["model"]

    centralRootFile100 = r.TFile.Open("./rootFiles/%s_v%d/sigScan_%s_%s_2012_100.0_bt0.0_MChi-1.0.root"%(settings["model"], versions[settings['model']], settings["model"], settings["sele"]))
    centralRootFile87 = r.TFile.Open("./rootFiles/%s_v%d/sigScan_%s_%s_2012_86.7_bt0.0_MChi-1.0.root"%(settings["model"], versions[settings['model']], settings["model"], settings["sele"]))
    centralRootFile73 = r.TFile.Open("./rootFiles/%s_v%d/sigScan_%s_%s_2012_73.3_bt0.0_MChi-1.0.root"%(settings["model"], versions[settings['model']], settings["model"], settings["sele"]))

    consts = {  'file73':centralRootFile73,
                'file87':centralRootFile87,
                'file100':centralRootFile100,
                'sitv':settings['SITV']}
    
    plots_to_print = []

    if 'list' in str(type(settings['bMulti'])) and 'list' in str(type(settings['jMulti'])):
        variables_to_iterate = ['bMulti', 'jMulti']
        variable_values = [settings[x] for x in variables_to_iterate ]
        argument_set = list(itertools.product(*variable_values))

        arg_list = [ dict(zip(variables_to_iterate, x)) for x in argument_set ]

        

        for kwargs in arg_list:
            kwargs.update(consts)

            # switch on individual plot printing (for documentation)
            if not settings["combine_output"]:
                kwargs["file_output"] = True

            if settings["run_mode"] == "eff_maps":
                plots_to_print.append( make_eff_map_plot(**kwargs) )
            elif settings["run_mode"] == "sitv_acceptance":
                kwargs["file_output"] = False
                
                kwargs["sitv"] = False
                no_sitv_eff = make_eff_map_plot(**kwargs)

                kwargs["sitv"] = True
                sitv_eff = make_eff_map_plot(**kwargs)

                plots_to_print.append( make_var_eff_map_plot(raw_eff = sitv_eff, variation_eff = no_sitv_eff, bMulti = kwargs["bMulti"], jMulti = kwargs["jMulti"],
                                                                file_output = False if settings["combine_output"] else True, test_name = settings["run_mode"]) )
            
            elif settings["run_mode"] == "leptVeto_acceptance":
                kwargs["file_output"] = False

                nom_eff = make_eff_map_plot(**kwargs)
                
                kwargs['file100'] = r.TFile.Open("./rootFiles/%s_v%d/LeptonVeto/sigScan_%s_%s_2012_100.0_bt0.0_MChi-1.0_noLeptVeto.root"%(settings["model"], versions[settings['model']], settings["model"], settings["sele"]))
                kwargs['file87'] = r.TFile.Open("./rootFiles/%s_v%d/LeptonVeto/sigScan_%s_%s_2012_86.7_bt0.0_MChi-1.0_noLeptVeto.root"%(settings["model"], versions[settings['model']], settings["model"], settings["sele"]))
                kwargs['file73'] = r.TFile.Open("./rootFiles/%s_v%d/LeptonVeto/sigScan_%s_%s_2012_73.3_bt0.0_MChi-1.0_noLeptVeto.root"%(settings["model"], versions[settings['model']], settings["model"], settings["sele"]))

                no_leptVeto_eff = make_eff_map_plot(**kwargs)
                
                plots_to_print.append( make_var_eff_map_plot(raw_eff = nom_eff, variation_eff = no_leptVeto_eff, bMulti = kwargs["bMulti"], jMulti = kwargs["jMulti"],
                                                                file_output = False if settings["combine_output"] else True, test_name = settings["run_mode"]) )

            elif settings["run_mode"] == "minbias_acceptance":
                kwargs["file_output"] = False

                nom_eff = make_eff_map_plot(**kwargs)

                kwargs['file100'] = r.TFile.Open("./rootFiles/%s_v%d/MinBias/sigScan_%s_%s_2012_100.0_bt0.0_MChi-1.0.root"%(settings["model"], versions[settings['model']], settings["model"], settings["sele"]))
                kwargs['file87'] = r.TFile.Open("./rootFiles/%s_v%d/MinBias/sigScan_%s_%s_2012_86.7_bt0.0_MChi-1.0.root"%(settings["model"], versions[settings['model']], settings["model"], settings["sele"]))
                kwargs['file73'] = r.TFile.Open("./rootFiles/%s_v%d/MinBias/sigScan_%s_%s_2012_73.3_bt0.0_MChi-1.0.root"%(settings["model"], versions[settings['model']], settings["model"], settings["sele"]))

                minbias_eff = make_eff_map_plot(**kwargs)

                plots_to_print.append( make_var_eff_map_plot(raw_eff = minbias_eff, variation_eff = nom_eff, bMulti = kwargs["bMulti"], jMulti = kwargs["jMulti"],
                                                                file_output = False if settings["combine_output"] else True, test_name = settings["run_mode"]) )                
            
            else:
                pass

        if settings["run_mode"] == "eff_maps":
            plots_to_print.append( make_eff_map_plot(file73 = centralRootFile73, file87 = centralRootFile87, file100 = centralRootFile100, no_cuts=True) )

    if settings["combine_output"] and settings["run_mode"] != "btag_compare":
        print_plot_list(list = plots_to_print)
    
    if settings["run_mode"] == "btag_compare":

        for jM in ["le3j", "ge4j", "ge2j"]:
            plots_to_print = []

            zero_eff = make_eff_map_plot(file73 = centralRootFile73, file87 = centralRootFile87, file100 = centralRootFile100,
                                            jMulti = jM, bMulti = "eq0b")
            one_eff = make_eff_map_plot(file73 = centralRootFile73, file87 = centralRootFile87, file100 = centralRootFile100,
                                            jMulti = jM, bMulti = "eq1b")
            two_eff = make_eff_map_plot(file73 = centralRootFile73, file87 = centralRootFile87, file100 = centralRootFile100,
                                            jMulti = jM, bMulti = "eq2b")

            c1 = r.TCanvas()

            div_one_zero = divide_eff_maps(num = one_eff, denom = zero_eff, title = "eq1b/eq0b, %s" % jM)
            div_two_one = divide_eff_maps(num = two_eff, denom = one_eff, title = "eq2b/eq1b, %s" % jM)

            div_one_zero.Draw("colz")
            div_one_zero.GetZaxis().SetRangeUser(0.,.8)
            c1.Print("out/%s_%s_one_zero_%s.pdf" % (settings['model'], settings['run_mode'], jM))

            div_two_one.Draw("colz")
            div_two_one.GetZaxis().SetRangeUser(0.,.8)
            c1.Print("out/%s_%s_two_one_%s.pdf" % (settings['model'], settings['run_mode'], jM))

