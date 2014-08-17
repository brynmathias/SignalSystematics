#!/usr/bin/env python
import sys
import os
import ROOT as r
from plottingstuff import *
from plottingUtils import Print, MakeCumu
import math

r.gROOT.SetBatch(r.kTRUE)

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

settings = {
    "syst": ["leptVeto", "deadEcal", "mhtmet"][0:3], # must be an array!
    "HTBins":["275_325", "325_375", "375_475", "475_575", "575_675", "675_775", "775_875", "875"],
    "SubProcesses":["nn","ns","ng","ss","ll","sb","tb","gg","bb","sg"],
    "deltaM":[False, True][0],
    "jMulti": ["le3j", "ge4j", "eq2j", "eq3j", "ge2j"][0],
    "bMulti": ["eq0b", "eq1b", "ge0b"][0]
}

def GetHist(File = None, folder = None, hist = None, Norm = None, rebin = None):
    h = None
    for f in folder:
      directory = File.Get(f)
      a = directory.Get(hist)
      if h is None:
          h = a.Clone()
      else:
        h.Add(a)
    return h
    
    
def nloTotalXsecMaker(individualXSecs = None):
    out = None
    for h in individualXSecs:
        if out is None: out = h.Clone()
        else: out.Add(h)
    return out

def getRootDirs():

    dirs = []
    bMulti = settings["bMulti"]
    jMulti = settings["jMulti"]

    ## convert individual string selections to lists for iteration
    if "str" in str(type(bMulti)):
        bMulti = [bMulti]
    if "str" in str(type(jMulti)):
        jMulti = [jMulti]

    for jM in jMulti:
        for bM in bMulti:
            for ht in settings["HTBins"]:
                dirs.append("smsScan_%s_%s_AlphaT55_%s" % (bM, jM, ht))

    return dirs

models = ["T2cc",]#"T1ttttProto",]#"T2bb","T2tt","T1tttt","T1bbbb"]

for cut in settings["syst"]:
  for model in models:

      xTitle = None
      yTitle = None
      if model == "T2cc":
        xTitle = "m_{scharm} (GeV)"
        if settings["deltaM"]:
          yTitle = "deltaM (GeV)"
        else:
          yTitle = "m_{LSP} (GeV)"
        maxi = .1
        mini = 1E-3
        m_gl_m_m_lsp = 0. # necessary to NOT remove the points for this scan, near the diagonal
        m_sq = 0.

      if model == "T2bb":
        xTitle = "m_{sbottom} (GeV)"
        yTitle = "m_{LSP} (GeV)"
        maxi = .1
        mini = 1E-3
        m_gl_m_m_lsp = 175.
        m_sq = 300.
      
      if model == "T2tt":
        xTitle = "m_{stop} (GeV)"
        yTitle = "m_{LSP} (GeV)"
        maxi = .1
        mini = 1E-3
        m_gl_m_m_lsp = 175.
        m_sq = 300.

      if model == "T1":
        xTitle = "m_{gluino} (GeV)"
        yTitle = "m_{LSP} (GeV)"
        maxi = .1
        mini = 1E-3
        m_gl_m_m_lsp = 175.
        m_sq = 300.
        
      if model == "T2":
        xTitle = "m_{squark} (GeV)"
        yTitle = "m_{LSP} (GeV)"
        maxi = .1
        mini = 1E-3
        m_gl_m_m_lsp = 175.
        m_sq = 300.

        
      if model == "T1tttt":
        xTitle = "m_{gluino} (GeV)"
        yTitle = "m_{LSP} (GeV)"
        maxi = .1
        mini = 1E-3
        m_gl_m_m_lsp = 375.
        m_sq = 450.

      if model == "T1ttttProto":
        xTitle = "m_{gluino} (GeV)"
        yTitle = "m_{LSP} (GeV)"
        maxi = .1
        mini = 1E-3
        m_gl_m_m_lsp = 375.
        m_sq = 450.


        
      if model == "T1bbbb":
        xTitle = "m_{gluino} (GeV)"
        yTitle = "m_{LSP} (GeV)"
        maxi = .1
        mini = 1E-3
        m_gl_m_m_lsp = 175.
        m_sq = 300.

      # set some limits for plotting ranges
      if "mhtmet" in cut:
        jmax = 0.4
      elif "leptVeto" in cut:
        jmax = 0.03 
      elif "deadEcal" in cut:
        jmax = 0.3
      else:
        jmax = 0.1

      outName = ""
      if settings["deltaM"]:
        outName = "./out/LepVeto_SMS%s_%s_dM.pdf"%(model, cut)
      else:
        outName = "./out/LepVeto_SMS%s_%s.pdf"%(model, cut)
    
      c1 = Print(outName)
      c1.DoPageNum = False
      #c1.open()
      r.gPad.SetRightMargin(0.175)
      r.gPad.SetLeftMargin(0.15)
      r.gPad.SetTopMargin(0.05)
      r.gPad.SetBottomMargin(0.15)

      FullCutFlowRootFile100 = r.TFile.Open("./rootfiles/T2cc_v8/sigScan_%s_had_2012_100.0_bt0.0_MChi-1.0.root"%(model))
      NMinus1RootFile100 = r.TFile.Open("./rootfiles/T2cc_v8/%s/sigScan_%s_had_2012_100.0_bt0.0_MChi-1.0.root"%(cut, model))
      FullCutFlowRootFile87 = r.TFile.Open("./rootfiles/T2cc_v8/sigScan_%s_had_2012_86.7_bt0.0_MChi-1.0.root"%(model))
      NMinus1RootFile87 = r.TFile.Open("./rootfiles/T2cc_v8/%s/sigScan_%s_had_2012_86.7_bt0.0_MChi-1.0.root"%(cut, model))
      FullCutFlowRootFile74 = r.TFile.Open("./rootfiles/T2cc_v8/sigScan_%s_had_2012_73.7_bt0.0_MChi-1.0.root"%(model))
      NMinus1RootFile74 = r.TFile.Open("./rootfiles/T2cc_v8/%s/sigScan_%s_had_2012_73.7_bt0.0_MChi-1.0.root"%(cut, model))


      # Make cross sections/ efficiencies
      
      processCrossSections = []
      cuts = []
      cutsJESNeg = []
      cutsJESPlus = []
      # cutsJESRan = []
      nocuts = []

      nocuts = GetHist(File = FullCutFlowRootFile100,folder = ["smsScan_before",], hist = "m0_m12_mChi_noweight", Norm = None ,rebin= 2)
      nocuts = threeToTwo(nocuts)

      FullCutFlow = GetHist(File = FullCutFlowRootFile100,folder = getRootDirs()[2:],hist = "m0_m12_mChi_noweight", Norm = None ,rebin= 2).Clone()
      FullCutFlow.Add( GetHist(File = FullCutFlowRootFile87,folder = getRootDirs()[1:2],hist = "m0_m12_mChi_noweight", Norm = None ,rebin= 2) )
      FullCutFlow.Add( GetHist(File = FullCutFlowRootFile74,folder = getRootDirs()[0:1],hist = "m0_m12_mChi_noweight", Norm = None ,rebin= 2) )                                                  

      NMinus1 = GetHist(File = NMinus1RootFile100,folder = getRootDirs()[2:],hist = "m0_m12_mChi_noweight", Norm = None ,rebin= 2).Clone()
      NMinus1.Add( GetHist(File = NMinus1RootFile87,folder = getRootDirs()[1:2],hist = "m0_m12_mChi_noweight", Norm = None ,rebin= 2) )
      NMinus1.Add( GetHist(File = NMinus1RootFile74,folder = getRootDirs()[0:1],hist = "m0_m12_mChi_noweight", Norm = None ,rebin= 2) )

      FullCutFlow = threeToTwo(FullCutFlow)
      NMinus1 = threeToTwo(NMinus1)

      if settings["deltaM"]:
        nocuts = deltaM(nocuts)
        FullCutFlow = deltaM(FullCutFlow)
        NMinus1 = deltaM(NMinus1)

      Range = NMinus1.GetYaxis().GetNbins()*NMinus1.GetXaxis().GetNbins()

      # zero values outside of the range of interest
      l =  [i * 25 for i in range(86)]
      for a in l:
        for b in l:
          # m_sq (m_gl) - m_LSP >= 175 && m_sq (m_gl) >= 300
          if  a - b <= m_gl_m_m_lsp or a < m_sq :
            bin = FullCutFlow.FindBin(float(a),float(b))
            if bin > 2114: continue
            # print bin
            FullCutFlow.SetBinContent(bin,0.)
            NMinus1.SetBinContent(bin,0.)          
    
      r.gStyle.SetOptStat(0)

      c1.canvas.SetLogz()
      offset = 1.1
      # c1.Print()    
      c1.canvas.SetLogz(False)

      nocuts.Draw("COLZ")
      #c1.Print()


      FullCutFlowEff = FullCutFlow.Clone()
      FullCutFlowEff.GetZaxis().SetTitle("Fraction of expected signal yield")
      FullCutFlowEff.GetZaxis().SetTitleOffset(offset)
      FullCutFlowEff.GetZaxis().SetTitleSize(0.05)
      FullCutFlowEff.GetXaxis().SetTitle(xTitle)
      FullCutFlowEff.SetTitle("Total Eff. (w/ cut)")

      # FullCutFlowEff.GetYaxis().SetLabelSize(0.04)
      FullCutFlowEff.GetYaxis().SetTitleOffset(1.3)
      FullCutFlowEff.GetYaxis().SetTitleSize(0.05)        
      FullCutFlowEff.GetYaxis().SetTitle(yTitle)
      # FullCutFlowEff.SetTitle("Efficiency (No JES)")
      FullCutFlowEff.Divide(nocuts)
      imax = 0.01
      FullCutFlowEff.SetMaximum(imax)
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
      NMinus1Eff.SetTitle("Total Eff. (wo/ cut)")
      NMinus1Eff.Divide(nocuts)
      NMinus1Eff.SetMaximum(imax)
      NMinus1Eff.Draw("COLZ")
      c1.Print()
      
      EffChange = r.TH2D(FullCutFlowEff)
      EffChange.Divide(NMinus1Eff)
      EffChange.GetZaxis().SetTitle("Inefficiency from %s cut"%cut)
      EffChange.SetTitle("deltaEff. w/wo")
      for bin in range(EffChange.GetNbinsX()*EffChange.GetNbinsY()+1000):
        if EffChange.GetBinContent(bin) > 0.:
          EffChange.SetBinContent(bin, EffChange.GetBinContent(bin)-1.)
        else:
          EffChange.SetBinContent(bin, -1000.)
      EffChange.SetMinimum(-1.*jmax)
      EffChange.SetMaximum(jmax)
      EffChange.Draw("COLZ")  
      c1.Print()

     
      oneDMap = r.TH2D(FullCutFlowEff)
      oneDMap.SetTitle("OneDMap")
      oneD = r.TH1D("OneD_Projection","",500,0,jmax)

      irange = EffChange.GetXaxis().GetNbins()*EffChange.GetYaxis().GetNbins()

      for bin in range(irange+1000):
        if EffChange.GetBinContent(bin) > (-1.*jmax) :
          oneD.Fill(math.fabs(EffChange.GetBinContent(bin)))
          print EffChange.GetBinContent(bin)

      scaloneD = oneD.Integral()
      oneD = MakeCumu(oneD)
      if scaloneD>0:
        oneD.Scale(1./scaloneD)
  
      oneD.Draw("h")

      oneD.GetXaxis().SetTitle("Inefficiency from cut")

      # now find 68%
      bin68 = 0
      for bin in range(oneD.GetNbinsX()):
        if oneD.GetBinContent(bin) <= 0.68:
          bin68 = bin
      oneDClone = r.TH1D(oneD)
      for bin in range(oneD.GetNbinsX()):
        if bin > bin68: 
          oneDClone.SetBinContent(bin,0.)
          oneDClone.SetBinError(bin,0.)
      num = r.TLatex(0.4,0.2,"68%% of events below %.3f"%(oneDClone.GetBinLowEdge(bin68)))
      num.SetNDC()
      num.Draw("same")

      oneDClone.SetFillColor(r.kRed)
      oneDClone.Draw("sameh")      
      c1.Print()

      

      countera = 0
      counterb = 0
      counterc = 0

      binlist = []
      closeBins = []
      farBins = []
      closeToLine = r.TH1D("OneD_Projection_closeToLine","",200,0,1.)
      farToLine = r.TH1D("OneD_Projection_farToLine","",200,0,1.)
      line = 350.
      cutoff = 450.
      # if model == "T1tttt":
      #   line = 550.
      #   cutoff = 600.
      for a in l:
        for b in l:
          # m_sq (m_gl) - m_LSP >= 175 && m_sq (m_gl) >= 300
          bin = EffChange.FindBin(float(a),float(b))
          if a - b <= m_gl_m_m_lsp or a < m_sq : continue
          if bin > 2114: continue
          if a > 1200: continue
          if b > 1200: continue
          oneDMap.SetBinContent(bin,1.)
          countera +=1
          if  a - b >= line and a > cutoff:                
            # print bin
            # if EffChange.GetBinContent(bin) > 0.:
              farBins.append(bin)
          else:
            # if EffChange.GetBinContent(bin) > 0.:
              closeBins.append(bin)


      print "Counter a = %d, counter b = %d, counter c = %d, b+c = %d"%(countera,counterb,counterc,counterc+counterb)
      oneD.GetXaxis().SetTitle("Fraction of expected signal yield rejected")

      # r.gPad.SetRightMargin(0.05)
      oneD.SetName("")

      oneDMap.Draw("COLZ")
      #c1.Print()
      allbins = closeBins+farBins 
      closeBins = set(closeBins)
      farBins = set(farBins)
      for b in farBins:
        if b in closeBins: closeBins.remove(b)
      TestMe = r.TH2D(EffChange)
      TestMe.SetMinimum(0.)
      TestMe.SetMaximum(1.)
      
      
      for b in farBins: farToLine.Fill(EffChange.GetBinContent(b))
      for b in closeBins: closeToLine.Fill(EffChange.GetBinContent(b))
      

      for b in range(EffChange.GetXaxis().GetNbins()*EffChange.GetYaxis().GetNbins()+1000):
        if EffChange.GetBinContent(b) > 0.: 
          print EffChange.GetBinContent(b)
          oneD.Fill(EffChange.GetBinContent(b)+1.)

      closeBinNom = closeToLine.Integral()
      farBinNorm = farToLine.Integral()
      closeToLineClone = MakeCumu(closeToLine)
      farToLineClone = MakeCumu(farToLine)
      print "Close Bin Norm = %f, farbinnom = %f"%(closeBinNom,farBinNorm)
      if closeBinNom > 0.:  closeToLineClone.Scale(1./closeBinNom)
      if farBinNorm > 0.:farToLineClone.Scale(1./farBinNorm)
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
      c1.canvas.Clear()
      closeToLine.GetXaxis().SetLabelSize(0.03)
      closeToLine.Draw("hist")
      closeToLine68.SetFillColor(r.kRed)
      closeToLine68.SetLineColor(r.kRed)
      closeToLine68.Draw("samehist")
      num = r.TLatex(0.4,0.8,"68%% of events below %.3f"%(closeToLineClone.GetBinLowEdge(bin68Close)))
      num.SetNDC()

      num.Draw("same")
    
      #c1.Print()
      farToLine68 = r.TH1D(farToLine)
      for bin in range(farToLine.GetNbinsX()):
         if bin > bin68far:
           farToLine68.SetBinContent(bin,0.)
           farToLine68.SetBinError(bin,0.)
      c1.canvas.Clear()
      farToLine.GetXaxis().SetLabelSize(0.03)
      farToLine.Draw("hist")
      farToLine68.SetFillColor(r.kRed)
      farToLine68.SetLineColor(r.kRed)
      farToLine68.Draw("samehist")
      num = r.TLatex(0.4,0.8,"68%% of events below %.3f"%(farToLineClone.GetBinLowEdge(bin68far)))
      num.SetNDC()
      num.Draw("same")
    
      #c1.Print()





      r.gStyle.SetOptStat(0)
      for b in range(TestMe.GetXaxis().GetNbins()*TestMe.GetYaxis().GetNbins()):
        # if b in closeBins: TestMe.SetBinContent(b,1.)
        if b in farBins: 
          TestMe.SetBinContent(b,0.5)
        elif b in closeBins: TestMe.SetBinContent(b,1.0)
        else:
          TestMe.SetBinContent(b,0.)
      
      TestMe.Draw("COLZ")
      #c1.Print()

      
      
      
      r.gPad.SetRightMargin(0.175)

      for bin in range(EffChange.GetNbinsX()*EffChange.GetNbinsY()):
        if bin not in allbins: continue
        # if EffChange.GetBinContent(bin) > 0.:
        if EffChange.GetBinContent(bin) > maxi: EffChange.SetBinContent(bin,maxi)
        if EffChange.GetBinContent(bin) < mini: EffChange.SetBinContent(bin,mini)          

      c1.canvas.Clear()
      EffChange.Draw("COLZ")
      EffChange.SetMinimum(mini)
      EffChange.SetMaximum(maxi)
      EffChange.GetZaxis().SetTitle("Fraction of expected signal yield rejected")
      c1.canvas.SetLogz()
      #c1.Print()
      

      r.gPad.SetLeftMargin(0.1)
      r.gPad.SetRightMargin(0.05)
      c1.canvas.SetLogz(False)
      c1.canvas.SetLogy(True)
      r.gStyle.SetOptStat(1111)
      oneD.Draw("h")
      #c1.Print()

      closeToLine.GetXaxis().SetTitle("Fraction of expected signal yield rejected")
      closeToLine.Draw("h")
      #c1.Print()
      farToLine.GetXaxis().SetTitle("Fraction of expected signal yield rejected")
      farToLine.Draw("h")
      #c1.Print()

      c1.canvas.SetLogy(False)
          



      c1.close()
