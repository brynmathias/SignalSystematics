#!/usr/bin/env python
import sys
import os
import ROOT as r
from plottingstuff import *
from plottingUtils import Print, MakeCumu
import math

r.gStyle.SetOptStat(0)
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


settings = {
    "HTBins":["a"],
    "SubProcesses":["nn","ns","ng","ss","ll","sb","tb","gg","bb","sg"]
}
def GetHist(File = None, folder = None, hist = None, Norm = None, rebin = None):
    h = None
    for f in folder:
        # print f
        # print hist
        directory = File.Get(f)
        # print directory.ls()
        a = directory.Get(hist)
        if h is None:
            h = a.Clone()
        else: h.Add(a)
    return h
    

subDirList = [
    "smsScan_geq1b_geq2j_AlphaT55_375_475",
    "smsScan_geq1b_geq2j_AlphaT55_475_575",
    "smsScan_geq1b_geq2j_AlphaT55_575_675",
    "smsScan_geq1b_geq2j_AlphaT55_675_775",
    "smsScan_geq1b_geq2j_AlphaT55_775_875",
    "smsScan_geq1b_geq2j_AlphaT55_875",
    "smsScan_eq0b_geq2j_AlphaT55_375_475",
    "smsScan_eq0b_geq2j_AlphaT55_475_575",
    "smsScan_eq0b_geq2j_AlphaT55_575_675",
    "smsScan_eq0b_geq2j_AlphaT55_675_775",
    "smsScan_eq0b_geq2j_AlphaT55_775_875",
    "smsScan_eq0b_geq2j_AlphaT55_875",
    
]

    
def nloTotalXsecMaker(individualXSecs = None):
    out = None
    for h in individualXSecs:
        if out is None: out = h.Clone()
        else: out.Add(h)
    return out
models = ["T1","T2"]
for model in models:

    xTitle = None
    yTitle = None
    if model == "T2bb":
      xTitle = "m_{sbottom} (GeV)"
      yTitle = "m_{LSP} (GeV)"
    if model == "T2tt":
      xTitle = "m_{stop} (GeV)"
      yTitle = "m_{LSP} (GeV)"
    if model == "T1":
      xTitle = "m_{gluino} (GeV)"
      yTitle = "m_{LSP} (GeV)"
    if model == "T2":
      xTitle = "m_{squark} (GeV)"
      yTitle = "m_{LSP} (GeV)"
    if model == "T1tttt":
      xTitle = "m_{gluino} (GeV)"
      yTitle = "m_{LSP} (GeV)"
    if model == "T1bbbb":
      xTitle = "m_{gluino} (GeV)"
      yTitle = "m_{LSP} (GeV)"
    if model == "T1ttttProto":
      xTitle = "m_{gluino} (GeV)"
      yTitle = "m_{LSP} (GeV)"
    
    
    
    c1 = Print("./out/BTagScale_SMS%s.pdf"%(model))
    c1.DoPageNum = False
    
    r.gPad.SetRightMargin(0.175)
    r.gPad.SetLeftMargin(0.15)
    r.gPad.SetTopMargin(0.05)
    r.gPad.SetBottomMargin(0.15)
    centalRootFile73 = r.TFile.Open("./rootFiles/sigScan_%s_had_2012_73.7_bt0.0_MChi-1.root"%(model))
    centalRootFile86 = r.TFile.Open("./rootFiles/sigScan_%s_had_2012_86.7_bt0.0_MChi-1.root"%(model))
    centalRootFile100 = r.TFile.Open("./rootFiles/sigScan_%s_had_2012_100.0_bt0.0_MChi-1.root"%(model))
    

    jesPlusRootFile73 = r.TFile.Open("./rootFiles/sigScan_%s_had_2012_73.7_bt1.0_MChi-1.root"%(model))
    jesPlusRootFile86 = r.TFile.Open("./rootFiles/sigScan_%s_had_2012_86.7_bt1.0_MChi-1.root"%(model))
    jesPlusRootFile100 = r.TFile.Open("./rootFiles/sigScan_%s_had_2012_100.0_bt1.0_MChi-1.root"%(model))
    

    jesNegRootFile73 = r.TFile.Open("./rootFiles/sigScan_%s_had_2012_73.7_bt-1.0_MChi-1.root"%(model))
    jesNegRootFile86 = r.TFile.Open("./rootFiles/sigScan_%s_had_2012_86.7_bt-1.0_MChi-1.root"%(model))
    jesNegRootFile100 = r.TFile.Open("./rootFiles/sigScan_%s_had_2012_100.0_bt-1.0_MChi-1.root"%(model))
      
    # jesRanRootFile73 = r.TFile.Open("./results_had_T2bb_73.7_ran.root")
    # jesRanRootFile86 = r.TFile.Open("./results_had_T2bb_86.7_ran.root")
    # jesRanRootFile100 = r.TFile.Open("./results_had_T2bb_100_ran.root")


    processCrossSections = []
    cuts = []
    cutsJESNeg = []
    cutsJESPlus = []
    # cutsJESRan = []
    nocuts = []
    # centalRootFile73.ls()
    # Make cross sections/ efficiencies
    for bin in settings["HTBins"]:        
        

        nocuts = GetHist(File = centalRootFile100,folder = ["smsScan_before",],
                         hist = "m0_m12_mChi_noweight", Norm = None ,rebin= 2)
        nocuts = threeToTwo(nocuts)
   
        cutsHist = GetHist(File = centalRootFile73,folder = ["smsScan_geq1b_geq2j_AlphaT55_275_325","smsScan_eq0b_geq2j_AlphaT55_275_325",],hist = "m0_m12_mChi_noweight", Norm = None ,rebin= 2).Clone()
        cutsHist.Add(GetHist(File = centalRootFile86,folder = ["smsScan_geq1b_geq2j_AlphaT55_325_375","smsScan_eq0b_geq2j_AlphaT55_325_375",],hist = "m0_m12_mChi_noweight", Norm = None ,rebin= 2))
        cutsHist.Add(GetHist(File = centalRootFile100,folder = subDirList,hist 
        = "m0_m12_mChi_noweight", Norm = None ,rebin= 2))
        cutsHist = threeToTwo(cutsHist)                                                    
        cutsJESPlusHist = GetHist(File =   jesPlusRootFile73,folder = ["smsScan_geq1b_geq2j_AlphaT55_275_325","smsScan_eq0b_geq2j_AlphaT55_275_325",],hist = "m0_m12_mChi_noweight", Norm = None ,rebin= 2).Clone()
        cutsJESPlusHist.Add(GetHist(File = jesPlusRootFile86,folder = ["smsScan_geq1b_geq2j_AlphaT55_325_375","smsScan_eq0b_geq2j_AlphaT55_325_375",],hist = "m0_m12_mChi_noweight", Norm = None ,rebin= 2))
        cutsJESPlusHist.Add(GetHist(File = jesPlusRootFile100,folder = 
        subDirList,hist = "m0_m12_mChi_noweight", Norm = None ,rebin= 2))

        cutsJESPlusHist = threeToTwo(cutsJESPlusHist)
        cutsJESNegHist = GetHist(File =   jesNegRootFile73,folder = ["smsScan_geq1b_geq2j_AlphaT55_275_325","smsScan_eq0b_geq2j_AlphaT55_275_325",],hist = "m0_m12_mChi_noweight", Norm = None ,rebin= 2).Clone()
        cutsJESNegHist.Add(GetHist(File = jesNegRootFile86,folder = [    
"smsScan_geq1b_geq2j_AlphaT55_325_375","smsScan_eq0b_geq2j_AlphaT55_325_375",]
,hist = "m0_m12_mChi_noweight", Norm = None ,rebin= 2))
        cutsJESNegHist.Add(GetHist(File = jesNegRootFile100,folder = 
        subDirList,hist = "m0_m12_mChi_noweight", Norm = None ,rebin= 2))
        cutsJESNegHist = threeToTwo(cutsJESNegHist)                                       
        
        

        
        # cutsJESRanHist = threeToTwo(cutsJESRanHist)
        cutsHist.GetXaxis().SetRangeUser(0.,1200.)
        cutsJESPlusHist.GetXaxis().SetRangeUser(0.,1200.)
        cutsJESNegHist.GetXaxis().SetRangeUser(0.,1200.)
        # cutsJESRanHist.GetXaxis().SetRangeUser(0.,1200.)
        cutsHist.GetYaxis().SetRangeUser(0.,1200.)
        cutsJESPlusHist.GetYaxis().SetRangeUser(0.,1200.)
        cutsJESNegHist.GetYaxis().SetRangeUser(0.,1200.)

        l =  [i * 5 for i in range(301)]
        # print l
        for a in l:
          for b in l:
            # m_sq (m_gl) - m_LSP >= 175 && m_sq (m_gl) >= 300
            if  a - b <= 175 or a < 300. :
              bin = cutsHist.FindBin(float(a),float(b))
              cutsHist.SetBinContent(bin,0.)
              cutsJESPlusHist.SetBinContent(bin,0.)
              cutsJESNegHist.SetBinContent(bin,0.)
              # cutsJESRanHist.SetBinContent(bin,0.)
          

        x_width = cutsHist.GetXaxis().GetBinWidth(10)
        y_width = cutsHist.GetYaxis().GetBinWidth(10)
        cutValue = 250.
        # print x_width,y_width
        #  use Y = mX + c
        mini = 0.8
        maxi =1.2
        

        c1.canvas.SetLogz()
        offset = 1.1
        # c1.Print()    
        c1.canvas.SetLogz(False)
        TotalEff = cutsHist.Clone()
        TotalEff.GetZaxis().SetTitle("Fraction of expected signal yield")
        TotalEff.GetZaxis().SetTitleOffset(offset)
        TotalEff.GetZaxis().SetTitleSize(0.05)
        TotalEff.GetXaxis().SetTitle(xTitle)

        # TotalEff.GetYaxis().SetLabelSize(0.04)
        TotalEff.GetYaxis().SetTitleOffset(1.3)
        TotalEff.GetYaxis().SetTitleSize(0.05)        
        TotalEff.GetYaxis().SetTitle(yTitle)
        # TotalEff.SetTitle("Efficiency (No JES)")
        TotalEff.Divide(nocuts)
        TotalEff.Draw("COLZ")
        # TotalEff.Scale(100.)
        c1.Print()
        # TotalEff.Scale(1./100.)
        TotalEffPlus = cutsJESPlusHist.Clone()
        TotalEffPlus.GetZaxis().SetTitle("Relative change in efficiency")
        TotalEffPlus.GetZaxis().SetTitleOffset(offset)
        TotalEffPlus.Divide(nocuts)
        TotalEffPlus.GetXaxis().SetTitle(xTitle)
        TotalEffPlus.GetYaxis().SetTitle(yTitle)
        TotalEffPlus.GetZaxis().SetTitleSize(0.05)
        
        # TotalEffPlus.GetYaxis().SetLabelSize(0.04)
        TotalEffPlus.GetYaxis().SetTitleOffset(1.3)
        TotalEffPlus.GetYaxis().SetTitleSize(0.05)    
        TotalEffPlus.Draw("COLZ")
        c1.Print()    
        TotalEffNeg = cutsJESNegHist.Clone()
        TotalEffNeg.GetZaxis().SetTitle("Relative change in efficiency")
        TotalEffNeg.GetZaxis().SetTitleOffset(offset)
        TotalEffNeg.Divide(nocuts)
        TotalEffNeg.GetZaxis().SetTitleSize(0.05)
        
        TotalEffNeg.GetXaxis().SetTitle(xTitle)
        TotalEffNeg.GetYaxis().SetTitle(yTitle)
        # TotalEffNeg.GetYaxis().SetLabelSize(0.04)
        TotalEffNeg.GetYaxis().SetTitleOffset(1.3)
        TotalEffNeg.GetYaxis().SetTitleSize(0.05)
        TotalEffNeg.Draw("COLZ")
        c1.Print()

        EffOverJESNeg = TotalEffNeg.Clone()
        EffOverJESNeg.Divide(TotalEff)
        TotalEffNeg.Draw("COLZ")


        EffOverJESPlus = TotalEffPlus.Clone()
        EffOverJESPlus.Divide(TotalEff)
        TotalEffPlus.Draw("COLZ")

        oneDJesMinus = r.TH1D("oneDJesMinus","",400,0.,0.4)
        oneDJesPlus = r.TH1D("oneDJesPlus","",400,0.,0.4)
        # oneDJesRan = r.TH1D("oneDJesRan","JES variation for signal efficiency 1D Projection",1000,-5.,5.)
        nEvents = r.TH1D("totEv","totEv",2500,0,2500)
        minf=0.9
        maxf=1.1
        fit =r.TF1("Gaussian","gaus",minf,maxf)
    
        totalBins = TotalEff.GetXaxis().GetNbins()*TotalEff.GetYaxis().GetNbins()
        for bin in range(totalBins):
            contentMinus = math.fabs(EffOverJESNeg.GetBinContent(bin)-1.)
            contentPlus =  math.fabs(EffOverJESPlus.GetBinContent(bin)-1.)
            # contentRan =   EffOverJESRan.GetBinContent(bin)
            content = nocuts.GetBinContent(bin)
            if content == 0: continue
            if EffOverJESPlus.GetBinContent(bin) < 0.01: continue
            if EffOverJESPlus.GetBinContent(bin) < 0.  or EffOverJESPlus.GetBinContent(bin) > 100.:
               nEvents.Fill(math.fabs(cutsJESPlusHist.GetBinContent(bin)-cutsHist.GetBinContent(bin)))
           
            if contentMinus > 0.:oneDJesMinus.Fill(math.fabs(contentMinus))
            if contentPlus > 0.: oneDJesPlus.Fill(math.fabs(contentPlus))


        l =  [i * 25 for i in range(60)]


        binlist = []
        closeBins = []
        farBins = []
        closeToLine = r.TH1D("OneD_Projection_closeToLine","",400,0,0.4)
        farToLine = r.TH1D("OneD_Projection_farToLine","",400,0,0.4)
        closeToLine.GetXaxis().SetTitle("Relative change in efficiency")
        farToLine.GetXaxis().SetTitle("Relative change in efficiency")
        line = 350.
        cutoff = 450.
        for a in l:
          for b in l:
            # m_sq (m_gl) - m_LSP >= 175 && m_sq (m_gl) >= 300
            bin = EffOverJESPlus.FindBin(float(a),float(b))
            if  a - b >= line and a > cutoff :                
              # print bin
              farBins.append(bin)
            else:
              closeBins.append(bin)


        # 
        closeBins = set(closeBins)
        farBins = set(farBins)
        for bin in farBins:
          if bin in closeBins: 
            closeBins.remove(bin)
            print "Bin %d is in both near and Far!!!!!!!"%bin
        TestMe = r.TH2D(EffOverJESPlus)
        TestMe.SetMinimum(0.)
        TestMe.SetMaximum(1.)
        for b in farBins: 
          if EffOverJESPlus.GetBinContent(b) >0.:farToLine.Fill(math.fabs(EffOverJESPlus.GetBinContent(b)-1.))
          if EffOverJESNeg.GetBinContent(b) >0.: farToLine.Fill(math.fabs(EffOverJESNeg.GetBinContent(b)-1.))
        for b in closeBins: 
          if EffOverJESPlus.GetBinContent(b) >0.:closeToLine.Fill(math.fabs(EffOverJESPlus.GetBinContent(b)-1.))
          if EffOverJESNeg.GetBinContent(b) >0.: closeToLine.Fill(math.fabs(EffOverJESNeg.GetBinContent(b)-1.))

        closeBinNom = closeToLine.Integral()
        farBinNorm = farToLine.Integral()
        closeToLineClone = MakeCumu(closeToLine)
        farToLineClone = MakeCumu(farToLine)
        if closeBinNom > 0.:
            closeToLineClone.Scale(1./closeBinNom)
        if farBinNorm > 0.:
            farToLineClone.Scale(1./farBinNorm)
        bin68Close = 0
        for bin in range(closeToLine.GetNbinsX()):
          # print model, closeToLine.GetBinContent(bin) , bin
          if closeToLineClone.GetBinContent(bin)  <= 0.68:
            bin68Close = bin
        bin68far = 0
        for bin in range(farToLine.GetNbinsX()):
          # print model, farToLine.GetBinContent(bin) , bin
          if farToLineClone.GetBinContent(bin)  <= 0.68:
            bin68far = bin

        # r.gStyle.SetOptStat(0)
        
        closeToLine68 = r.TH1D(closeToLine)
        for bin in range(closeToLine.GetNbinsX()):
           if bin > bin68Close:
             closeToLine68.SetBinContent(bin,0.)
             closeToLine68.SetBinError(bin,0.)
        closeToLine.Draw("hist")
        closeToLine68.SetFillColor(r.kRed)
        closeToLine68.SetLineColor(r.kRed)
        closeToLine68.Draw("samehist")
        num = r.TLatex(0.4,0.8,"68%% of events below %.3f"%(closeToLineClone.GetBinLowEdge(bin68Close)))
        num.SetNDC()
        num.Draw("same")
        
        c1.Print()
        farToLine68 = r.TH1D(farToLine)
        for bin in range(farToLine.GetNbinsX()):
           if bin > bin68far:
             farToLine68.SetBinContent(bin,0.)
             farToLine68.SetBinError(bin,0.)
        
        farToLine.Draw("hist")
        farToLine68.SetFillColor(r.kRed)
        farToLine68.SetLineColor(r.kRed)
        farToLine68.Draw("samehist")
        num = r.TLatex(0.4,0.8,"68%% of events below %.3f"%(farToLineClone.GetBinLowEdge(bin68far)))
        num.SetNDC()
        num.Draw("same")
        
        c1.Print()

        r.gStyle.SetOptStat(0)
        

        for b in range(TestMe.GetXaxis().GetNbins()*TestMe.GetYaxis().GetNbins()):
          # if b in closeBins: TestMe.SetBinContent(b,1.)
          if b in farBins: 
            TestMe.SetBinContent(b,0.5)
          elif b in closeBins: TestMe.SetBinContent(b,1.0)
          else:
            TestMe.SetBinContent(b,0.)
          if EffOverJESNeg.GetBinContent(b) < 0.01: TestMe.SetBinContent(b,0.)
        TestMe.Draw("COLZ")
        r.gStyle.SetOptStat(0)
        c1.Print()





        r.gStyle.SetOptStat(0)
        scalValMinus = oneDJesMinus.Integral()
        # oneDJesMinus.GetXaxis().SetRangeUser(0.8,1.2)
        oneDJesMinus = MakeCumu(oneDJesMinus)
        oneDJesMinus.Scale(1./scalValMinus)
        oneDJesMinus.Draw("hist")
        oneDJesMinus.GetXaxis().SetTitle("Relative change in efficiency")
        oneDJesMinus.GetXaxis().SetTitleSize(0.05)
        oneDJesMinus.Fit("Gaussian","RQ")
        fit.Draw("lsame")

        # c1.Print()
        # oneDJesPlus.GetXaxis().SetRangeUser(0.8,1.2)
        scalValPlus = oneDJesPlus.Integral()
        oneDJesPlus = MakeCumu(oneDJesPlus)
        oneDJesPlus.Scale(1./scalValPlus)
        oneDJesPlus.Draw("hist")
        oneDJesPlus.GetXaxis().SetTitle("Relative change in efficiency")
        oneDJesPlus.GetXaxis().SetTitleSize(0.05)
        oneDJesPlus.Fit("Gaussian","RQ")
        fit.Draw("lsame")

        c1.Print()
 
        
        
        JesTotal = oneDJesPlus.Clone()
        JesTotal.Add(oneDJesMinus)
        JesTotal.Scale(1./2.)
        JesTotal.SetName("JESTotal")
        bin68 = 0
        for bin in range(JesTotal.GetNbinsX()):
          # print model, JesTotal.GetBinContent(bin) , bin
          if JesTotal.GetBinContent(bin) <= 0.68:
            bin68 = bin
        JesTotClone = r.TH1D(JesTotal)
        for bin in range(JesTotal.GetNbinsX()):
          if bin > bin68: 
            JesTotClone.SetBinContent(bin,0.)
            JesTotClone.SetBinError(bin,0.)
        # JesTotal.SetTitle("Total (Sum of up and down)")
        JesTotal.Draw("h")
        JesTotClone.SetFillColor(r.kRed)
        JesTotClone.Draw("sameh")
        JesTotal.Fit("Gaussian","RQ")
        fit.Draw("lsame")
        # print "RMS is =", oneDJesRan.GetRMS()
        c1.Print()
        r.gStyle.SetOptStat(0)
        for bin in range(EffOverJESNeg.GetNbinsX()*EffOverJESNeg.GetNbinsY()):
          if EffOverJESNeg.GetBinContent(bin) > 0:
            if EffOverJESNeg.GetBinContent(bin) < mini: EffOverJESNeg.SetBinContent(bin,mini)
            if EffOverJESNeg.GetBinContent(bin) > maxi: EffOverJESNeg.SetBinContent(bin,maxi)
        for bin in range(EffOverJESPlus.GetNbinsX()*EffOverJESPlus.GetNbinsY()):
          if EffOverJESNeg.GetBinContent(bin) > 0:
            if EffOverJESPlus.GetBinContent(bin) < mini: EffOverJESPlus.SetBinContent(bin,mini)
            if EffOverJESPlus.GetBinContent(bin) > maxi: EffOverJESPlus.SetBinContent(bin,maxi)
            
        for bin in range(EffOverJESPlus.GetXaxis().GetNbins()*EffOverJESPlus.GetYaxis().GetNbins()+10000):
          if EffOverJESNeg.GetBinContent(bin) > 0:
            EffOverJESNeg.SetBinContent(bin,EffOverJESNeg.GetBinContent(bin)-1.)
          else:EffOverJESNeg.SetBinContent(bin,-1000)
          if EffOverJESPlus.GetBinContent(bin) > 0:
            EffOverJESPlus.SetBinContent(bin,EffOverJESPlus.GetBinContent(bin)-1.)
          else:EffOverJESPlus.SetBinContent(bin,-1000)
        EffOverJESNeg.SetMinimum(mini-1.)
        EffOverJESNeg.SetMaximum(maxi-1.)
        EffOverJESPlus.SetMinimum(mini-1.)
        EffOverJESPlus.SetMaximum(maxi-1.)

        EffOverJESNeg.Draw("COLZ")
        c1.Print()
        EffOverJESPlus.Draw("COLZ")
        c1.Print()
        
        


    c1.close()
