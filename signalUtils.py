import ROOT as r
from array import array
import plotDetails as pdeets
import numpy as np
import math

r.gROOT.SetBatch(r.kTRUE)

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
        # print f, h
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


def getRootDirs(bMulti_="", jMulti_="", htbins = [], sitv_=False):

    dirs = []


    for ht in htbins:
        thisHT = ht.split("_")
        if len(thisHT)<2: thisHT.append(None)
        aT = pdeets.alphaTDict["%s,%s" % (thisHT[0], thisHT[1])][0]
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
                    # nbins, int(minVal), int(maxVal)+10.)
                    nbins+2, 0., int(maxVal)+10.)
    print "delta mass bin hack!!!!"

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

###-------------------------------------------------------------------###

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

###-------------------------------------------------------------------###

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
def reject_outliers(data, m = 2.):
    """return listen of indices to reject"""

    reject = []
    new = [np.abs(d - np.median(data)) for d in data]
    mdev = np.median(new)
    snew = []
    for d in new:
        if mdev:
            snew.append(d/mdev)
        else:
            snew.append(0.)
    for n in range(len(snew)):
        if snew[n]>m:
            reject.append(n)
        if snew[n] == 0.:
            reject.append(n)

    return reject

def syst_smooth(eff = None, err = None, iterations = 1):
    """smooth the eff"""

    print "Smoothing%s, %d iterations" % (" with err" if err else "", iterations)

    new_hist = eff.Clone()
    itera = 0
    # note: ONLY WORKS FOR T2CC AND T24BODY because of binning assumptions
    while itera < iterations:
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
                            # get relative err
                            tmp_err = float(err.GetBinContent( xbin+xtmp, ybin+xtmp*5)/tmp_val)
                        else:
                            tmp_err = 0.
                        vals.append(tmp_val)
                        if tmp_err > 0.:
                            errs.append(float(1./math.pow(tmp_err, 1)))
                        else:
                            errs.append(0.)

                for ytmp in range(-1,2):
                    tmp_val = eff.GetBinContent(xbin, ybin+ytmp*2)
                    if tmp_val > 0.:
                        if err:
                            tmp_err = float(err.GetBinContent(xbin, ybin+ytmp*2)/tmp_val)
                        else:
                            tmp_err = 0.
                        vals.append(tmp_val)
                        if tmp_err > 0.:
                            errs.append(float(1./math.pow(tmp_err, 1)))
                        else:
                            errs.append(0.)
                # get weighted average considering errs
                err_sum = np.sum(errs)
                vals_array = np.array(vals)
                to_rej = reject_outliers(vals_array)

                new_vals = []
                new_errs = []
                for n in range(len(vals)):
                    if n not in to_rej:
                        new_vals.append(vals[n])
                        new_errs.append(errs[n])
                err_sum = np.sum(new_errs)
                if err_sum>0.:
                    for j in range(len(new_errs)):
                        new_errs[j]/=err_sum
                    ave_val = np.average(new_vals, weights=new_errs)
                else:
                    # should only be in here in err == None
                    ave_val = np.average(new_vals)

                if ave_val == np.nan:
                    ave_val = 1.
                new_hist.SetBinContent(xbin, ybin, ave_val)
        itera+=1

        eff = new_hist.Clone()
    return new_hist

###-------------------------------------------------------------------###

def dict_printer(dicto = {}, indent = 1):
  
  print "{ (%d keys)\n" % len(dicto)
  for key in dicto:
    print "\t"*indent, "'%s': " % key,
    if dict == type(dicto[key]):
      dict_printer(dicto[key], indent+1)
    else:
      print dicto[key]
  print "\t"*indent, "}\n"

###-------------------------------------------------------------------###
###-------------------------------------------------------------------###