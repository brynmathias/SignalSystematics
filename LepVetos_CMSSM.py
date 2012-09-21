#!/usr/bin/env python
import sys
import os
import ROOT as r
from plottingstuff import *
from plottingUtils import *
from pointsAlongLimit import *
import math


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
maxi = .1
mini = 1E-3
    
def nloTotalXsecMaker(individualXSecs = None):
    out = None
    for h in individualXSecs:
        if out is None: out = h.Clone()
        else: out.Add(h)
    return out
models = ["tanB10",]
Cuts = ["MHTovMET","DECAL"]
Dirs = ["mSuGraScan_AlphaT55_375_475_scale1","mSuGraScan_AlphaT55_475_575_scale1","mSuGraScan_AlphaT55_575_675_scale1","mSuGraScan_AlphaT55_675_775_scale1","mSuGraScan_AlphaT55_775_875_scale1","mSuGraScan_AlphaT55_875_scale1"]
Dirs2 = ["mSuGraScan_noOddLep_AlphaT55_375_475_scale1","mSuGraScan_noOddLep_AlphaT55_475_575_scale1","mSuGraScan_noOddLep_AlphaT55_575_675_scale1","mSuGraScan_noOddLep_AlphaT55_675_775_scale1","mSuGraScan_noOddLep_AlphaT55_775_875_scale1","mSuGraScan_noOddLep_AlphaT55_875_scale1"]
xTitle = "m_{0}"
yTitle = "m_{1/2}"
for cut in Cuts:
  for model in models:  
      c1 = Print("./out/LEPVETO_CMSSM.pdf")
      c1.DoPageNum = False
      c1.open()
      r.gPad.SetRightMargin(0.175)
      r.gPad.SetLeftMargin(0.15)
      r.gPad.SetTopMargin(0.05)
      r.gPad.SetBottomMargin(0.15)
      print model, cut
      centalRootFile100 = "./LepVetoOnSignalFiles/results_had_%s_100.root"%(model)
      jesPlusRootFile100 = "./LepVetoOnSignalFiles/results_had_%s_100.root"%(model)

      processCrossSections = []
      cuts = []
      cutsNoWeight = []
      cutsJESNeg = []
      cutsJESPlus = []
      # cutsJESRan = []
      nocuts = []
      

      for process in settings["SubProcesses"]:                
          p_xsec = GetSumHist(File = [centalRootFile100],Directories = ["mSuGraScan_before_scale1",], Hist = "m0_m12_%s"%(process),Norm = None ,LegendText = None).hObj
          p_xsec.Divide(GetSumHist(File = [centalRootFile100],Directories = ["mSuGraScan_before_scale1",], Hist = "m0_m12_%s_noweight"%(process),Norm = None ,LegendText = None).hObj)
          p_xsec.SetTitle(process)
          processCrossSections.append(p_xsec)
          nocutsHist = GetSumHist(File = [centalRootFile100],Directories = ["mSuGraScan_before_scale1",], Hist = "m0_m12_%s"%(process), Norm = None ,LegendText = None).hObj

          cutsHist=GetSumHist(File = [centalRootFile100],Directories =Dirs,Hist = "m0_m12_%s"%(process), Norm = None ,LegendText = None).hObj.Clone()

          cutsHistNoWeight=GetSumHist(File = [centalRootFile100],Directories = Dirs, Hist = "m0_m12_%s_noweight"%(process), Norm = None ,LegendText = None).hObj.Clone()

          cutsJESPlusHist=GetSumHist(File = [jesPlusRootFile100],Directories = Dirs2,Hist = "m0_m12_%s"%(process), Norm = None ,LegendText = None).hObj.Clone()



          b = None
          nocuts.append(nocutsHist)
          cuts.append(cutsHist)
          cutsJESPlus.append(cutsJESPlusHist)   
          cutsNoWeight.append(cutsHistNoWeight)                                                        
      


      totalXsec =  nloTotalXsecMaker(processCrossSections)
      FullCutFlow =  r.TH2D(NloEffHisto(cuts,nocuts,processCrossSections,totalXsec))           
      for bin in range(FullCutFlow.GetXaxis().GetNbins()*FullCutFlow.GetYaxis().GetNbins() + 1000):
        if FullCutFlow.GetBinContent(bin) > 0.9: FullCutFlow.SetBinContent(bin,0.)
      NMinus1 = r.TH2D(NloEffHisto(cutsJESPlus,nocuts,processCrossSections,totalXsec))
      for bin in range(NMinus1.GetXaxis().GetNbins()*NMinus1.GetYaxis().GetNbins() + 1000):
        if NMinus1.GetBinContent(bin) > 0.9: NMinus1.SetBinContent(bin,0.)

                              
        

        

      c1.canvas.SetLogz()
      offset = 1.1
      # c1.Print()    
      c1.canvas.SetLogz(False)
      FullCutFlowEff = FullCutFlow.Clone()
      FullCutFlowEff.GetZaxis().SetTitle("Fraction of expected signal yield")
      FullCutFlowEff.GetZaxis().SetTitleOffset(offset)
      FullCutFlowEff.GetZaxis().SetTitleSize(0.05)
      FullCutFlowEff.GetXaxis().SetTitle(xTitle)

      # FullCutFlowEff.GetYaxis().SetLabelSize(0.04)
      FullCutFlowEff.GetYaxis().SetTitleOffset(1.3)
      FullCutFlowEff.GetYaxis().SetTitleSize(0.05)        
      FullCutFlowEff.GetYaxis().SetTitle(yTitle)
      # FullCutFlowEff.SetTitle("Efficiency (No JES)")
      # FullCutFlowEff.Divide(nocuts)
      FullCutFlowEff.Draw("COLZ")
      c1.Print()



      NMinus1Eff = NMinus1.Clone()
      NMinus1Eff.GetZaxis().SetTitle("Fraction of expected signal yield")
      NMinus1Eff.GetZaxis().SetTitleOffset(offset)
      NMinus1Eff.GetZaxis().SetTitleSize(0.05)
      NMinus1Eff.GetXaxis().SetTitle(xTitle)
      NMinus1Eff.GetYaxis().SetTitleOffset(1.3)
      NMinus1Eff.GetYaxis().SetTitleSize(0.05)        
      NMinus1Eff.GetYaxis().SetTitle(yTitle)
      # NMinus1Eff.Divide(nocuts)
      NMinus1Eff.Draw("COLZ")
      c1.Print()
      
      EffChange = r.TH2D(FullCutFlowEff)
      EffChange.Divide(NMinus1Eff)
      EffChange.GetZaxis().SetTitle("Inefficiency from %s cut"%cut)
      for bin in range(EffChange.GetNbinsX()*EffChange.GetNbinsY()):
        if EffChange.GetBinContent(bin) > 0.:
          EffChange.SetBinContent(bin, 1.- EffChange.GetBinContent(bin))
      
          
          

      # c1.Print()
         
      oneDaroundLim = r.TH1D("OneD_Projection_AroundLimit","",500,0,0.5)

      oneD = r.TH1D("OneD_Projection","",500,0,0.5)
      for bin in range(EffChange.GetNbinsX()*EffChange.GetNbinsY()):
        events = 0
        for h in cutsNoWeight:
          events+=h.GetBinContent(bin)
        if events >0.: oneD.Fill(EffChange.GetBinContent(bin))
  
      c1.canvas.SetLogy(True)
      oneD.GetXaxis().SetTitle("Fraction of expected signal yield rejected")
      r.gStyle.SetOptStat(1111)
      r.gPad.SetRightMargin(0.05)
      oneD.SetName("")
      oneD.Draw("h")
      c1.Print()

      # Now look at the error around the limit
      oneDaroundLim = r.TH1D("OneD_Projection_AroundLimit","",200,0,0.2)
      limx=[100,300,500,700,900,1100,1300,1500,1700,1900,2000,3000]
      limy=[640,625,600,560,480, 400, 340, 330, 320, 320, 320, 320]
      bins = returnBinsAroundLimit(histogram = EffChange, limitX = limx, limitY =limy, nBinsAroundLimit = 5)
      bins = set(bins)
      limit = ExclusionCurve(xVals = limx, yVals = limy)
      limit.SetLineColor(r.kBlack)

      for bin in bins:
        events = 0
        for h in cutsNoWeight:
          events+=h.GetBinContent(bin)
        if events >25: oneDaroundLim.Fill(EffChange.GetBinContent(bin))
      oneDaroundLim.GetXaxis().SetTitle("Fraction of expected signal yield rejected")

      oneDaroundLim.Draw("h")
      c1.Print()


      
      # c1.Print()
      r.gStyle.SetOptStat(0)
      
      r.gPad.SetRightMargin(0.175)

      for bin in range(EffChange.GetXaxis().GetNbins()*EffChange.GetYaxis().GetNbins() + 10000):
        events = 0
        for h in cutsNoWeight:
          events+=h.GetBinContent(bin)
        if events > 25:
          if bin > EffChange.GetXaxis().GetNbins()*EffChange.GetYaxis().GetNbins(): EffChange.SetBinContent(bin,mini)
          if EffChange.GetBinContent(bin) > maxi: EffChange.SetBinContent(bin,maxi)
          if EffChange.GetBinContent(bin) < mini: EffChange.SetBinContent(bin,mini)          
      c1.canvas.SetLogy(False)

      c1.canvas.Clear()
      EffChange.Draw("COLZ")
      EffChange.SetMinimum(mini)
      EffChange.SetMaximum(maxi)
      EffChange.GetZaxis().SetTitle("Fraction of expected signal yield rejected")
      limit.Draw("same")
      c1.canvas.SetLogz()
      c1.Print()

      for b in range(EffChange.GetXaxis().GetNbins()*EffChange.GetXaxis().GetNbins()):
            
        if b not in bins: 
          EffChange.SetBinContent(b,0.)
      EffChange.Draw("COLZ")
      limit.Draw("same")
      c1.Print()



      c1.close()
