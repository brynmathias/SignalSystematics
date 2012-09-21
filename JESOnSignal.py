#!/usr/bin/env python
import sys
import os
import ROOT as r
from plottingstuff import *
from plottingUtils import Print,MakeCumu
import math
from pointsAlongLimit import *

rmbins = []



def AddHits(histList = None):
  """docstring for AddHits"""
  out = None
  for h in histList:
    if out is None: out = h.Clone()
    else: out.Add(h)
  return out
  

# Make central efficiency map
def threeToTwo(h3) :
    name = h3.GetName()
    binsz = h3.GetNbinsZ()
    #print binsz
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
    "HTBins":["a",],#"mSuGraScan_AlphaT55_375_475_scale1","mSuGraScan_AlphaT55_475_575_scale1","mSuGraScan_AlphaT55_575_675_scale1",
                    #         "mSuGraScan_AlphaT55_675_775_scale1","mSuGraScan_AlphaT55_775_875_scale1","mSuGraScan_AlphaT55_875_scale1"],
    "SubProcesses":["nn","ns","ng","ss","ll","sb","tb","gg","bb","sg"]
}
def GetHist(File = None, folder = None, hist = None, Norm = None, rebin = None):
    h = None
    for f in folder:
        #print f
        #print hist
        directory = File.Get(f)    
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

c1 = Print("./out/JES_CMSSM_NEWSCAN.pdf")

c1.DoPageNum = False

c1.open()
r.gPad.SetRightMargin(0.175)

centalRootFile73 = r.TFile.Open("./JESandLepVetoRootFiles/results_NEWSCAN_JESCEN_had_tanB10_73.7.root")
centalRootFile86 = r.TFile.Open("./JESandLepVetoRootFiles/results_NEWSCAN_JESCEN_had_tanB10_86.7.root")
centalRootFile100 = r.TFile.Open("./JESandLepVetoRootFiles/results_NEWSCAN_JESCEN_had_tanB10_100.root")
jesPlusRootFile73 = r.TFile.Open("./JESandLepVetoRootFiles/results_NEWSCAN_had_tanB10_73.7_+ve.root")
jesPlusRootFile86 = r.TFile.Open("./JESandLepVetoRootFiles/results_NEWSCAN_had_tanB10_86.7_+ve.root")
jesPlusRootFile100 = r.TFile.Open("./JESandLepVetoRootFiles/results_NEWSCAN_had_tanB10_100_+ve.root")
jesNegRootFile73 = r.TFile.Open("./JESandLepVetoRootFiles/results_NEWSCAN_had_tanB10_73.7_-ve.root")
jesNegRootFile86 = r.TFile.Open("./JESandLepVetoRootFiles/results_NEWSCAN_had_tanB10_86.7_-ve.root")
jesNegRootFile100 = r.TFile.Open("./JESandLepVetoRootFiles/results_NEWSCAN_had_tanB10_100_-ve.root")
# jesRanRootFile73 = r.TFile.Open("./results_NEWSCAN_had_tanB10_73.7_ran.root")
# jesRanRootFile86 = r.TFile.Open("./results_NEWSCAN_had_tanB10_86.7_ran.root")
# jesRanRootFile100 = r.TFile.Open("./results_NEWSCAN_had_tanB10_100_ran.root")



centalRootFile100.ls()
# Make cross sections/ efficiencies
RawEvNums = GetHist(File = centalRootFile100,folder = ["mSuGraScan_before_scale1",], hist = "m0_m12_mChi_noweight",Norm = 8. ,rebin= 2)
RawEvNums = threeToTwo(RawEvNums)
BinsToExclude = []
for bin in range(RawEvNums.GetNbinsX()*RawEvNums.GetNbinsY()):
    if RawEvNums.GetBinContent(bin) != 10000 : BinsToExclude.append(bin)
Dirs = ["mSuGraScan_AlphaT55_375_475_scale1","mSuGraScan_AlphaT55_475_575_scale1","mSuGraScan_AlphaT55_575_675_scale1","mSuGraScan_AlphaT55_675_775_scale1","mSuGraScan_AlphaT55_775_875_scale1","mSuGraScan_AlphaT55_875_scale1"]
for HTbin in settings["HTBins"]:        
    processCrossSections = []
    cuts = []
    cutsNoWeight = []
    cutsJESNeg = []
    cutsJESPlus = []
    # cutsJESRan = []
    nocuts = []
    for process in settings["SubProcesses"]:                
        p_xsec = GetHist(File = centalRootFile100,folder = ["mSuGraScan_before_scale1",], hist = "m0_m12_%s"%(process),Norm = 8. ,rebin= 2)
        p_xsec.Divide(GetHist(File = centalRootFile100,folder = ["mSuGraScan_before_scale1",], hist = "m0_m12_%s_noweight"%(process),Norm = 8. ,rebin= 2))
        p_xsec.SetTitle(process)
        processCrossSections.append(p_xsec)
        nocutsHist = GetHist(File = centalRootFile100,folder = ["mSuGraScan_before_scale1",], hist = "m0_m12_%s"%(process), Norm = None ,rebin= 2)



        # cutsHist = GetHist(File = centalRootFile73,folder = ["mSuGraScan_AlphaT55_275_325_scale1"],hist = "m0_m12_%s"%(process), Norm = None ,rebin= 2).Clone()
        # cutsHist.Add(GetHist(File = centalRootFile86,folder = ["mSuGraScan_AlphaT55_325_375_scale1"],hist = "m0_m12_%s"%(process), Norm = None ,rebin= 2))
        cutsHist=GetHist(File = centalRootFile100,folder =Dirs,hist = "m0_m12_%s"%(process), Norm = None ,rebin= 2).Clone()
        cutsHistNoWeight=GetHist(File = centalRootFile100,folder =Dirs,hist = "m0_m12_%s_noweight"%(process), Norm = None ,rebin= 2).Clone()

        # cutsHistNoWeight = GetHist(File = centalRootFile73,folder = ["mSuGraScan_AlphaT55_275_325_scale1"],hist = "m0_m12_%s_noweight"%(process), Norm = None ,rebin= 2).Clone()
        # cutsHistNoWeight.Add(GetHist(File = centalRootFile86,folder = ["mSuGraScan_AlphaT55_325_375_scale1"],hist = "m0_m12_%s_noweight"%(process), Norm = None ,rebin= 2))
        cutsHistNoWeight=GetHist(File = centalRootFile100,folder = Dirs, hist = "m0_m12_%s_noweight"%(process), Norm = None ,rebin= 2).Clone()



        
        # cutsJESPlusHist = GetHist(File =   jesPlusRootFile73,folder = ["mSuGraScan_AlphaT55_275_325_scale1"],hist = "m0_m12_%s"%(process), Norm = None ,rebin= 2).Clone()
        # cutsJESPlusHist.Add(GetHist(File = jesPlusRootFile86,folder = ["mSuGraScan_AlphaT55_325_375_scale1"],hist = "m0_m12_%s"%(process), Norm = None ,rebin= 2))
        cutsJESPlusHist=GetHist(File = jesPlusRootFile100,folder = Dirs,hist = "m0_m12_%s"%(process), Norm = None ,rebin= 2).Clone()
  
                                                             


        # cutsJESNegHist = GetHist(File =   jesNegRootFile73,folder = ["mSuGraScan_AlphaT55_275_325_scale1"],hist = "m0_m12_%s"%(process), Norm = None ,rebin= 2).Clone()
        # cutsJESNegHist.Add(GetHist(File = jesNegRootFile86,folder = ["mSuGraScan_AlphaT55_325_375_scale1"],hist = "m0_m12_%s"%(process), Norm = None ,rebin= 2))
        cutsJESNegHist=GetHist(File = jesNegRootFile100,folder = Dirs,hist = "m0_m12_%s"%(process), Norm = None ,rebin= 2).Clone()

        
        


        b = None
        nocuts.append(nocutsHist)
        cuts.append(cutsHist)
        # cutsJESRan.append(cutsJESRanHist)                                                     
        cutsJESNeg.append(cutsJESNegHist)                                                             
        cutsJESPlus.append(cutsJESPlusHist)   
        cutsNoWeight.append(cutsHistNoWeight)                                                        
        cutsHistNoWeight.Draw("COLZ")
        c1.canvas.SetLogz()
        cutsHistNoWeight.SetTitle("%s Central"%(process))
        b = cutsHistNoWeight.Clone()
        for bin in range(b.GetNbinsX()*b.GetNbinsY()):
          if b.GetBinContent(bin)> 0: b.SetBinContent(bin,math.sqrt(b.GetBinContent(bin))/b.GetBinContent(bin))
        c1.Print()
        b.SetTitle("%s Relitive Error"%(process))
        b.Draw("COLZ")
        c1.Print()

    
    
    a = AddHits(cutsNoWeight)
    a.Draw("COLZ")
    c1.canvas.SetLogz()
    c1.Print()
    ents = r.TH1D("ents","ents",200,0,200)
    fracError = a.Clone()
    for bin in range(a.GetNbinsX()*a.GetNbinsY()):
      if a.GetBinContent(bin) > 0:   ents.Fill(a.GetBinContent(bin))
    ents.Draw("hist")
    c1.Print()
    for bin in range(fracError.GetNbinsX()*fracError.GetNbinsY()):
      if fracError.GetBinContent(bin) > 0: fracError.SetBinContent(bin, math.sqrt(fracError.GetBinContent(bin))/fracError.GetBinContent(bin))
    fracError.Draw("COLZ")
    fracError.SetMinimum(1E-2)
    c1.Print()
    c1.canvas.SetLogz(False)
    # for bin in range(a.GetNbinsX()*a.GetNbinsY()):
    #   if a.GetBinContent(bin) < 9: rmbins.append(bin)
    totalXsec =  nloTotalXsecMaker(processCrossSections)
    for bin in BinsToExclude+rmbins:
        totalXsec.SetBinContent(bin,0.)
    # for bin in rmbins:
        # totalXsec.SetBinContent(bin,0.00000000)

    c1.canvas.SetLogz()
    a = None
    TotalEff =  r.TH2D(NloEffHisto(cuts,nocuts,processCrossSections,totalXsec))
    # TotalEff.GetXaxis().SetTitle("m_0 (GeV)")
    # TotalEff.GetXaxis().SetTitle("m_1/2 (GeV)")
    noCutsTot = AddHits(nocuts).Clone()

    for p_xs, c , nc in zip(processCrossSections,cuts,nocuts):
      c1.canvas.SetLogz(False)
      processEfficiency = c.Clone()
      processEfficiency.Divide(noCutsTot)
      processEfficiency.Divide(TotalEff)
      processEfficiency.Draw("COLZ")
      processEfficiency.SetTitle("Fraction of total efficiency contributed by subprocess = "+ p_xs.GetTitle())
      processEfficiency.SetMaximum(1)
      processEfficiency.SetMinimum(0)
      c1.Print()
      a = p_xs.Clone()
      
      a.Divide(totalXsec)
      a.Multiply(processEfficiency)
      a.SetTitle("Fraction of total cross section "+a.GetTitle())
      a.Draw("COLZ")
      a.SetMinimum(0)
      a.SetMaximum(1)
      c1.Print()
    
    totalXsec.Draw("COLZ")
    totalXsec.SetTitle("xs")
    totalXsec.SetMinimum(1E-2)
    totalXsec.SetMaximum(1E5)    
    c1.Print()    
    c1.canvas.SetLogz(False)
    TotalEff =  r.TH2D(NloEffHisto(cuts,nocuts,processCrossSections,totalXsec))
    TotalEff.SetTitle("Efficiency")
    c1.canvas.SetLogz(False)
    
    # for b in rmbins:
    #   TotalEff.SetBinContent(bin,0.)
    # TotalEff.Multiply(totalXsec)
    c1.canvas.SetLogz()
    # TotalEff.Scale(4650.)
    # TotalEff.SetMaximum(1)
    for bin in range(TotalEff.GetNbinsX()*TotalEff.GetNbinsY()):
      if TotalEff.GetBinContent(bin) > 0.9: TotalEff.SetBinContent(bin,0.)
    
    TotalEff.SetMinimum(1E-3)
    TotalEff.Draw("COLZ")
    c1.Print()
    Yeild = TotalEff.Clone()
    c1.Print()
    Yeild.Multiply(totalXsec)
    c1.canvas.SetLogz()
    Yeild.Scale(4650.)
    Yeild.SetTitle("Yield for 4.65fb^{-1}")
    Yeild.Draw("COLZ")
    Yeild.SetMinimum(1)
    Yeild.SetMaximum(1E5)
    c1.Print()
    c1.canvas.SetLogz(False)
    
    TotalEffPlus = NloEffHisto(cutsJESPlus,nocuts,processCrossSections,totalXsec).Clone()
    for bin in range(TotalEffPlus.GetNbinsX()*TotalEffPlus.GetNbinsY()):
      if TotalEffPlus.GetBinContent(bin) > 0.8: TotalEffPlus.SetBinContent(bin,0.)


    TotalEffNeg =  NloEffHisto(cutsJESNeg,nocuts,processCrossSections,totalXsec).Clone()
    # TotalEffRan =  NloEffHisto(cutsJESRan,nocuts,processCrossSections,totalXsec)
    for bin in range(TotalEffNeg.GetNbinsX()*TotalEffNeg.GetNbinsY()):
      if TotalEffNeg.GetBinContent(bin) > 0.8: TotalEffNeg.SetBinContent(bin,0.)


    EffOverJESNeg = r.TH2D(TotalEff)
    EffOverJESNeg.Divide(TotalEffNeg)

    c1.Print()

    EffOverJESPlus = r.TH2D(TotalEff)
    EffOverJESPlus.Divide(TotalEffPlus)

    c1.Print()

    oneDJesMinus = r.TH1D("oneDJesMinus","",250,0.,2.5)
    oneDJesPlus = r.TH1D("oneDJesPlus","",250,0.,2.5)
    # oneDJesRan = r.TH1D("oneDJesRan","JES variation for signal efficiency 1D Projection",250,0.,2.5)

    minf=0.9
    maxf=1.1
    fit =r.TF1("Gaussian","gaus",minf,maxf)
    
    totalBins = TotalEff.GetXaxis().GetNbins()*TotalEff.GetYaxis().GetNbins()
    # #print BinsToExclude
    for bin in range(totalBins):
        if bin in BinsToExclude: continue
        contentMinus = EffOverJESNeg.GetBinContent(bin)
        contentPlus =  EffOverJESPlus.GetBinContent(bin)
        # contentRan =   EffOverJESRan.GetBinContent(bin)
        if contentMinus > 0.:oneDJesMinus.Fill(contentMinus,1./fracError.GetBinContent(bin)**2)
        if contentPlus > 0.: oneDJesPlus.Fill(contentPlus,1./fracError.GetBinContent(bin)**2)


    # c1.Print()
    minf=0.9
    maxf=1.1
    fit =r.TF1("Gaussian","gaus",minf,maxf)



    # c1.canvas.SetLogy()
    r.gStyle.SetOptStat(1111)
    oneDJesMinus.GetXaxis().SetRangeUser(0.2,1.8)
    oneDJesMinus.Draw("hist")
    oneDJesMinus.GetXaxis().SetTitle("Relative change in efficiency")
    oneDJesMinus.GetXaxis().SetTitleSize(0.05)
    
    # oneDJesMinus.Fit("Gaussian","RQ")
    # fit.Draw("lsame")

    c1.Print()
    oneDJesPlus.GetXaxis().SetRangeUser(0.2,1.8)
    oneDJesPlus.Draw("hist")
    # oneDJesPlus.SetName("")
    oneDJesPlus.GetXaxis().SetTitle("Relative change in efficiency")
    oneDJesPlus.GetXaxis().SetTitleSize(0.05)
    
    # oneDJesPlus.Fit("Gaussian","RQ")
    # fit.Draw("lsame")
    c1.Print()


    # oneDJesRan.GetXaxis().SetRangeUser(0.8,1.2)
    # oneDJesRan.Draw("hist")
    # oneDJesRan.Fit("Gaussian","RQ")
    # fit.Draw("lsame")
    # #print "RMS is =", oneDJesRan.GetRMS()


    JesTotal = oneDJesPlus.Clone()
    JesTotal.Add(oneDJesMinus)
    JesTotal.Draw("hist")
    JesTotal.SetName("JESTotal")
    # JesTotal.Fit("Gaussian","RQ")
    # fit.Draw("lsame")
    # #print "RMS is =", oneDJesRan.GetRMS()
    c1.Print()


    # Now look at the error around the limit
    oneDaroundLim = r.TH1D("OneD_Projection_AroundLimit","",400,0.,0.4)
    oneDaroundLim.GetXaxis().SetTitle("Relative change in efficiency")
    oneDaroundLim.GetXaxis().SetTitleSize(0.05)
    
    limx=[100,300,500,700,900,1100,1300,1500,1700,1900,2000,3000]
    limy=[640,625,600,560,480, 400, 340, 330, 320, 320, 320, 320]
    
    limit = ExclusionCurve(xVals = limx, yVals = limy)
    limit.SetLineColor(r.kBlack)


    twoDaroundLim = r.TH2D("2dLook","",200,0.9,1.1,100,0.,500.) 
    
    bins = returnBinsAroundLimit(histogram = EffOverJESPlus, limitX = limx, limitY =limy, nBinsAroundLimit = 3)
    bins = set(bins)
    for bin in bins:      
      cont = 0
      for h in cutsNoWeight:
        cont += h.GetBinContent(bin)
      if cont < 25.: continue
      if EffOverJESPlus.GetBinContent(bin) > 0.:
        oneDaroundLim.Fill(math.fabs(EffOverJESPlus.GetBinContent(bin)-1.))
        twoDaroundLim.Fill(EffOverJESPlus.GetBinContent(bin),cont)        
      if EffOverJESNeg.GetBinContent(bin) > 0. :
        oneDaroundLim.Fill(math.fabs(EffOverJESNeg.GetBinContent(bin)-1.))
        twoDaroundLim.Fill(EffOverJESNeg.GetBinContent(bin),cont)
    
    oneDCumu = MakeCumu(oneDaroundLim)
    oneDCumu.Scale(1./oneDaroundLim.Integral())
    bin68 = 0
    for bin in range(oneDaroundLim.GetNbinsX()):
      if oneDCumu.GetBinContent(bin) <= .68:
        bin68 = bin
    oneDaroundLimClone = r.TH1D(oneDaroundLim)
    for bin in range(oneDaroundLim.GetNbinsX()):
      if bin > bin68: 
        oneDaroundLimClone.SetBinContent(bin,0.)
        oneDaroundLimClone.SetBinError(bin,0.)
    oneDaroundLimClone.SetLineColor(r.kRed)
    oneDaroundLimClone.SetFillColor(r.kRed)
    r.gStyle.SetOptStat(0)
    oneDaroundLim.Draw("h")
    num = r.TLatex(0.4,0.8,"68%% of events below %.3f"%(oneDaroundLim.GetBinLowEdge(bin68)))
    num.SetNDC()
    num.Draw("same")
    oneDaroundLimClone.Draw("samehist")
    c1.Print()
    
    
    
    
    PointsConsidered = r.TH2D(EffOverJESNeg.Clone())
    PointsConsidered.SetTitle("")
    for bin in range(PointsConsidered.GetXaxis().GetNbins()*PointsConsidered.GetYaxis().GetNbins() + 1000):
      PointsConsidered.SetBinContent(bin,0.)
      if bin in bins: PointsConsidered.SetBinContent(bin,1.)
      
    r.gStyle.SetOptStat(0)
    twoDaroundLim.Draw("COLZ")
    c1.Print()
    PointsConsidered.Draw("COLZ")
    PointsConsidered.GetXaxis().SetTitle("m_{0} (GeV)")
    PointsConsidered.GetYaxis().SetTitle("m_{1/2} (GeV)")
    limit.Draw("SAME")
    c1.Print()





    mini = 0.8
    maxi = 1.2

    for bin in range(EffOverJESNeg.GetNbinsX()*EffOverJESNeg.GetNbinsY()):
      if EffOverJESNeg.GetBinContent(bin) > 0:
        if EffOverJESNeg.GetBinContent(bin) < mini: EffOverJESNeg.SetBinContent(bin,mini)
        if EffOverJESNeg.GetBinContent(bin) > maxi: EffOverJESNeg.SetBinContent(bin,maxi)
    for bin in range(EffOverJESPlus.GetNbinsX()*EffOverJESPlus.GetNbinsY()):
      if EffOverJESNeg.GetBinContent(bin) > 0:
        if EffOverJESPlus.GetBinContent(bin) < mini: EffOverJESPlus.SetBinContent(bin,mini)
        if EffOverJESPlus.GetBinContent(bin) > maxi: EffOverJESPlus.SetBinContent(bin,maxi)
    EffOverJESNeg.SetTitle("")
    EffOverJESPlus.SetTitle("")


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


    # 
    # EffOverJESNeg.SetMinimum(mini)
    # EffOverJESNeg.SetMaximum(maxi)
    # EffOverJESPlus.SetMinimum(mini)
    # EffOverJESPlus.SetMaximum(maxi)





    EffOverJESPlus.GetXaxis().SetTitle("m_{0} (GeV)")
    EffOverJESPlus.GetYaxis().SetTitle("m_{1/2} (GeV)")
    EffOverJESNeg.GetXaxis().SetTitle("m_{0} (GeV)")
    EffOverJESNeg.GetYaxis().SetTitle("m_{1/2} (GeV)")
    TotalEff.GetXaxis().SetTitle("m_{0} (GeV)")
    TotalEff.GetYaxis().SetTitle("m_{1/2} (GeV)")
    EffOverJESPlus.GetZaxis().SetTitle("Relative change in efficiency")
    EffOverJESNeg.GetZaxis().SetTitle("Relative change in efficiency")
    TotalEff.GetZaxis().SetTitle("Fraction of expected signal yield")

    
    TotalEff.GetZaxis().SetTitleSize(0.05)
    EffOverJESNeg.GetZaxis().SetTitleSize(0.05)
    EffOverJESPlus.GetZaxis().SetTitleSize(0.05)
    TotalEff.GetZaxis().CenterTitle()
    EffOverJESNeg.GetZaxis().CenterTitle()
    EffOverJESPlus.GetZaxis().CenterTitle()
    #print TotalEff.GetZaxis().GetTitleColor()
    #print bins    
    #print "Total bins =", EffOverJESPlus.GetXaxis().GetNbins()*EffOverJESPlus.GetXaxis().GetNbins()
    # for b in range(EffOverJESPlus.GetXaxis().GetNbins()*EffOverJESPlus.GetXaxis().GetNbins()):
      # if EffOverJESPlus.GetBinContent(b) > 1.01 or EffOverJESPlus.GetBinContent(b) < .99:
      #   EffOverJESPlus.SetBinContent(b,0.)
      #   # TotalEff.SetBinContent(b,0.)
      # 
      # if EffOverJESNeg.GetBinContent(b) > 1.01 or EffOverJESNeg.GetBinContent(b) < .99:
      #   EffOverJESNeg.SetBinContent(b,0.)
      #   # TotalEff.SetBinContent(b,0.)

      # if b not in bins: 
      #   EffOverJESPlus.SetBinContent(b,0.)
      #   EffOverJESNeg.SetBinContent(b,0.)
        # TotalEff.SetBinContent(b,0.)
    TotalEff.SetTitle("")
    a =  GetHist(File = centalRootFile100,folder = ["mSuGraScan_before_scale1",], hist = "m0_m12_mChi_noweight",Norm = 8. ,rebin= 2)
    a = threeToTwo(a)
    a.Draw("COLZ")
    a.SetMinimum(9990)
    c1.Print()
    TotalEff.Draw("COLZ")
    limit.Draw("same")
    c1.Print()
    # raw_input()
    EffOverJESNeg.Draw("COLZ")
    limit.Draw("same")
    c1.Print()
    EffOverJESPlus.Draw("COLZ")
    limit.Draw("same")
    c1.Print()
    
    




c1.close()
