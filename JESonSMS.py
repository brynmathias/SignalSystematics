#!/usr/bin/env python
import sys
import os
import ROOT as r
from plottingstuff import *
from plottingUtils import Print
import math
# Make central efficiency map

# def threeToTwo(h3) :
#     name = h3.GetName()
#     binsz = h3.GetNbinsZ()
#     print binsz
#     xBins = int(h3.GetXaxis().GetXmax()/h3.GetXaxis().GetBinWidth(1))
#     yBins = int(h3.GetYaxis().GetXmax()/h3.GetYaxis().GetBinWidth(1))
#     h2 = r.TH2D(name+"_2D",h3.GetTitle(),
#                 xBins, 0, h3.GetXaxis().GetXmax(),
#                 yBins, 0, h3.GetYaxis().GetXmax(),
#                 )
#                 
#     for iX in range(1, xBins) :
#         for iY in range(1, yBins) :
#             content = h3.GetBinContent(iX, iY, 1) + h3.GetBinContent(iX, iY, 2)+ h3.GetBinContent(iX, iY, 0)
#             h2.SetBinContent(iX, iY, content)
#     h2.GetZaxis().SetTitle(h3.GetZaxis().GetTitle())
#     return h2




def threeToTwo(h3) :
    name = h3.GetName()
    binsz = h3.GetNbinsZ()
    print binsz
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


settings = {
    "HTBins":["a"],
    "SubProcesses":["nn","ns","ng","ss","ll","sb","tb","gg","bb","sg"]
}
def GetHist(File = None, folder = None, hist = None, Norm = None, rebin = None):
    h = None
    for f in folder:
        print f
        print hist
        directory = File.Get(f)
        # print directory.ls()
        a = directory.Get(hist)
        if h is None:
            h = a.Clone()
        else: h.Add(a)
    return h
    
    
def nloTotalXsecMaker(individualXSecs = None):
    out = None
    for h in individualXSecs:
        if out is None: out = h.Clone()
        else: out.Add(h)
    return out

c1 = Print("JESSignalSMST1_250.pdf")
c1.open()
r.gPad.SetRightMargin(0.15)
centalRootFile73 = r.TFile.Open("./rootFiles/results_had_T1_73.7_.root")
centalRootFile86 = r.TFile.Open("./rootFiles/results_had_T1_86.7_.root")
centalRootFile100 = r.TFile.Open("./rootFiles/results_had_T1_100_.root")
jesPlusRootFile73 = r.TFile.Open("./rootFiles/results_had_T1_73.7_+ve.root")
jesPlusRootFile86 = r.TFile.Open("./rootFiles/results_had_T1_86.7_+ve.root")
jesPlusRootFile100 = r.TFile.Open("./rootFiles/results_had_T1_100_+ve.root")
jesNegRootFile73 = r.TFile.Open("./rootFiles/results_had_T1_73.7_-ve.root")
jesNegRootFile86 = r.TFile.Open("./rootFiles/results_had_T1_86.7_-ve.root")
jesNegRootFile100 = r.TFile.Open("./rootFiles/results_had_T1_100_-ve.root")
# jesRanRootFile73 = r.TFile.Open("./results_had_T1_73.7_ran.root")
# jesRanRootFile86 = r.TFile.Open("./results_had_T1_86.7_ran.root")
# jesRanRootFile100 = r.TFile.Open("./results_had_T1_100_ran.root")



# centalRootFile73.ls()
# Make cross sections/ efficiencies
for bin in settings["HTBins"]:        
    processCrossSections = []
    cuts = []
    cutsJESNeg = []
    cutsJESPlus = []
    # cutsJESRan = []
    nocuts = []

    nocuts = GetHist(File = centalRootFile73,folder = ["smsScan_before",], hist = "m0_m12_mChi_noweight", Norm = None ,rebin= 2)
    nocuts = threeToTwo(nocuts)
    # for bin in range(nocuts.GetNbinsX()*nocuts.GetNbinsY()):
    #   if nocuts.GetBinContent(bin) < 40000:
    #     nocuts.SetBinContent(bin,0.)
    cutsHist = GetHist(File = centalRootFile73,folder = ["smsScan_AlphaT55_275_325"],hist = "m0_m12_mChi_noweight", Norm = None ,rebin= 2).Clone()
    cutsHist.Add(GetHist(File = centalRootFile86,folder = ["smsScan_AlphaT55_325_375"],hist = "m0_m12_mChi_noweight", Norm = None ,rebin= 2))
    cutsHist.Add(GetHist(File = centalRootFile100,folder = ["smsScan_AlphaT55_375_475","smsScan_AlphaT55_475_575","smsScan_AlphaT55_575_675",
                                                         "smsScan_AlphaT55_675_775","smsScan_AlphaT55_775_875","smsScan_AlphaT55_875"],hist = "m0_m12_mChi_noweight", Norm = None ,rebin= 2))
    cutsHist = threeToTwo(cutsHist)                                                    
    cutsJESPlusHist = GetHist(File =   jesPlusRootFile73,folder = ["smsScan_AlphaT55_275_325"],hist = "m0_m12_mChi_noweight", Norm = None ,rebin= 2).Clone()
    cutsJESPlusHist.Add(GetHist(File = jesPlusRootFile86,folder = ["smsScan_AlphaT55_325_375"],hist = "m0_m12_mChi_noweight", Norm = None ,rebin= 2))
    cutsJESPlusHist.Add(GetHist(File = jesPlusRootFile100,folder = ["smsScan_AlphaT55_375_475","smsScan_AlphaT55_475_575","smsScan_AlphaT55_575_675",
                                                         "smsScan_AlphaT55_675_775","smsScan_AlphaT55_775_875","smsScan_AlphaT55_875"],hist = "m0_m12_mChi_noweight", Norm = None ,rebin= 2))

    cutsJESPlusHist = threeToTwo(cutsJESPlusHist)
    cutsJESNegHist = GetHist(File =   jesNegRootFile73,folder = ["smsScan_AlphaT55_275_325"],hist = "m0_m12_mChi_noweight", Norm = None ,rebin= 2).Clone()
    cutsJESNegHist.Add(GetHist(File = jesNegRootFile86,folder = ["smsScan_AlphaT55_325_375"],hist = "m0_m12_mChi_noweight", Norm = None ,rebin= 2))
    cutsJESNegHist.Add(GetHist(File = jesNegRootFile100,folder = ["smsScan_AlphaT55_375_475","smsScan_AlphaT55_475_575","smsScan_AlphaT55_575_675",
                                                         "smsScan_AlphaT55_675_775","smsScan_AlphaT55_775_875","smsScan_AlphaT55_875"],hist = "m0_m12_mChi_noweight", Norm = None ,rebin= 2))
    cutsJESNegHist = threeToTwo(cutsJESNegHist)                                       
        
        
    # cutsJESRanHist = GetHist(File =   jesRanRootFile73,folder = ["smsScan_AlphaT55_275_325"],hist = "m0_m12_mChi_noweight", Norm = None ,rebin= 2).Clone()
    #   cutsJESRanHist.Add(GetHist(File = jesRanRootFile86,folder = ["smsScan_AlphaT55_325_375"],hist = "m0_m12_mChi_noweight", Norm = None ,rebin= 2))
    #   cutsJESRanHist.Add(GetHist(File = jesRanRootFile100,folder = ["smsScan_AlphaT55_375_475","smsScan_AlphaT55_475_575","smsScan_AlphaT55_575_675",
    #                                                        "smsScan_AlphaT55_675_775","smsScan_AlphaT55_775_875","smsScan_AlphaT55_875"],hist = "m0_m12_mChi_noweight", Norm = None ,rebin= 2))
        
    # cutsJESRanHist = threeToTwo(cutsJESRanHist)
    cutsHist.GetXaxis().SetRangeUser(0.,1200.)
    cutsJESPlusHist.GetXaxis().SetRangeUser(0.,1200.)
    cutsJESNegHist.GetXaxis().SetRangeUser(0.,1200.)
    # cutsJESRanHist.GetXaxis().SetRangeUser(0.,1200.)
    cutsHist.GetYaxis().SetRangeUser(0.,1200.)
    cutsJESPlusHist.GetYaxis().SetRangeUser(0.,1200.)
    cutsJESNegHist.GetYaxis().SetRangeUser(0.,1200.)
    # cutsJESRanHist.GetYaxis().SetRangeUser(0.,1200.)
    
        # cutsHist.Draw("COLZ")
        # c1.Print()
        # cutsJESPlusHist.Draw("COLZ")
        # c1.Print()
        # cutsJESNegHist.Draw("COLZ")
        # c1.Print()
        # cutsJESRanHist.Draw("COLZ")
        # c1.Print()
    l =  [i * 5 for i in range(301)]
    print l
    for a in l:
      for b in l:
        if a - b < 250 - 75 :
          bin = cutsHist.FindBin(float(a),float(b))
          cutsHist.SetBinContent(bin,0.)
          cutsJESPlusHist.SetBinContent(bin,0.)
          cutsJESNegHist.SetBinContent(bin,0.)
          # cutsJESRanHist.SetBinContent(bin,0.)
          

    x_width = cutsHist.GetXaxis().GetBinWidth(10)
    y_width = cutsHist.GetYaxis().GetBinWidth(10)
    cutValue = 250.
    print x_width,y_width
    #  use Y = mX + c
    mini = 0.9
    maxi =1.1
        

    c1.canvas.SetLogz()

    # c1.Print()    
    c1.canvas.SetLogz(False)
    TotalEff = cutsHist.Clone()
    TotalEff.SetTitle("Efficiency (No JES)")
    TotalEff.Divide(nocuts)
    TotalEff.Draw("COLZ")
    c1.Print()
    TotalEffPlus = cutsJESPlusHist.Clone()
    TotalEffPlus.Divide(nocuts)
    TotalEffNeg = cutsJESNegHist.Clone()
    TotalEffNeg.Divide(nocuts)
    # TotalEffRan = cutsJESRanHist.Clone()
    # TotalEffRan.Divide(nocuts)


    EffOverJESNeg = TotalEffNeg.Clone()
    EffOverJESNeg.Divide(TotalEff)
    EffOverJESNeg.SetTitle("EffOverJESNeg")
    EffOverJESNeg.Draw("COLZ")
    for bin in range(EffOverJESNeg.GetNbinsX()*EffOverJESNeg.GetNbinsY()):
        if EffOverJESNeg.GetBinContent(bin) > maxi or EffOverJESNeg.GetBinContent(bin) < mini: EffOverJESNeg.SetBinContent(bin,0.)
    
    EffOverJESNeg.SetMinimum(mini)
    EffOverJESNeg.SetMaximum(maxi)
    c1.Print()

    EffOverJESPlus = TotalEffPlus.Clone()
    EffOverJESPlus.Divide(TotalEff)
    EffOverJESPlus.SetTitle("EffOverJESPlus")
    EffOverJESPlus.Draw("COLZ")
    for bin in range(EffOverJESPlus.GetNbinsX()*EffOverJESPlus.GetNbinsY()):
        if EffOverJESPlus.GetBinContent(bin) > maxi or EffOverJESPlus.GetBinContent(bin) < mini: EffOverJESPlus.SetBinContent(bin,0.)
    
    EffOverJESPlus.SetMinimum(mini)
    EffOverJESPlus.SetMaximum(maxi)
    c1.Print()

    # EffOverJESRan = TotalEffRan.Clone()
    # EffOverJESRan.Divide(TotalEff)
    # EffOverJESRan.SetTitle("JES variation for signal efficiency")
    # EffOverJESRan.Draw("COLZ")
    # for bin in range(EffOverJESRan.GetNbinsX()*EffOverJESRan.GetNbinsY()):
        # if EffOverJESRan.GetBinContent(bin) > maxi or EffOverJESRan.GetBinContent(bin) < mini: EffOverJESRan.SetBinContent(bin,0.)
    # EffOverJESRan.SetMinimum(mini)
    # EffOverJESRan.SetMaximum(maxi)
    # c1.Print()


    





    oneDJesMinus = r.TH1D("oneDJesMinus","oneDJesMinus",1000,-5.,5.)
    oneDJesPlus = r.TH1D("oneDJesPlus","oneDJesPlus",1000,-5.,5.)
    # oneDJesRan = r.TH1D("oneDJesRan","JES variation for signal efficiency 1D Projection",1000,-5.,5.)
    nEvents = r.TH1D("totEv","totEv",2500,0,2500)
    minf=0.9
    maxf=1.1
    fit =r.TF1("Gaussian","gaus",minf,maxf)
    
    totalBins = TotalEff.GetNbinsX()*TotalEff.GetNbinsY()
    for bin in range(totalBins):
        contentMinus = EffOverJESNeg.GetBinContent(bin)
        contentPlus =  EffOverJESPlus.GetBinContent(bin)
        # contentRan =   EffOverJESRan.GetBinContent(bin)
        content = nocuts.GetBinContent(bin)
        if content == 0: continue
        if EffOverJESPlus.GetBinContent(bin) < 0.01: continue
        if EffOverJESPlus.GetBinContent(bin) < 0.95  or EffOverJESPlus.GetBinContent(bin) > 1.15:
           print "Efficiency is %f, noscaling events = %f, scaling events = %f, difference"%(EffOverJESPlus.GetBinContent(bin),cutsHist.GetBinContent(bin),cutsJESPlusHist.GetBinContent(bin) )
           nEvents.Fill(math.fabs(cutsJESPlusHist.GetBinContent(bin)-cutsHist.GetBinContent(bin)))
           
        if contentMinus > 0.:oneDJesMinus.Fill(abs(contentMinus))
        if contentPlus > 0.: oneDJesPlus.Fill(abs(contentPlus))
        # if contentRan > 0.: oneDJesRan.Fill(contentRan)
    # nEvents.Draw()
    # c1.Print()
    # c1.canvas.SetLogy()
    r.gStyle.SetOptStat(1111)
    # oneDJesMinus.GetXaxis().SetRangeUser(0.8,1.2)
    oneDJesMinus.Draw("hist")
    oneDJesMinus.Fit("Gaussian","RQ")
    fit.Draw("lsame")

    c1.Print()
    # oneDJesPlus.GetXaxis().SetRangeUser(0.8,1.2)
    oneDJesPlus.Draw("hist")
    oneDJesPlus.Fit("Gaussian","RQ")
    fit.Draw("lsame")

    c1.Print()
    # minf=0.97
    # maxf=1.03
    # fit =r.TF1("Gaussian","gaus",minf,maxf)

    # oneDJesRan.GetXaxis().SetRangeUser(0.8,1.2)
    # oneDJesRan.Draw("hist")
    # oneDJesRan.Fit("Gaussian","RQ")
    # fit.Draw("lsame")
    # print "RMS is =", oneDJesRan.GetRMS()
    c1.Print()

    JesTotal = oneDJesPlus.Clone()
    JesTotal.Add(oneDJesMinus)
    JesTotal.SetTitle("Total (Sum of up and down)")
    JesTotal.Draw("hist")
    JesTotal.Fit("Gaussian","RQ")
    fit.Draw("lsame")
    # print "RMS is =", oneDJesRan.GetRMS()
    c1.Print()
    


c1.close()
