#!/usr/bin/env python
import sys
import os
import ROOT as r
from plottingstuff import *
from plottingUtils import Print

# Make central efficiency map

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

c1 = Print("EffectsOfLepVeto_GenLepVeto_T2tt.pdf")
c1.open()
r.gPad.SetRightMargin(0.15)

# centalRootFile100 = r.TFile.Open("./NMinusOneT2tt.root")




# centalRootFile73.ls()
# Make cross sections/ efficiencies

nocuts = threeToTwo(GetHist(File = r.TFile.Open("./results_GenLepVeto_had_T2tt_100.root"),folder = ["smsScan_before",], hist = "m0_m12_mChi_noweight", Norm = None ,rebin= 2))
    
nocuts = threeToTwo(nocuts)
AllCuts_GenLep = threeToTwo(GetHist(File = r.TFile.Open("results_GenLepPhotonVeto_had_T2tt_100.root"),folder = ["smsScan_AlphaT55_375_475","smsScan_AlphaT55_475_575","smsScan_AlphaT55_575_675",
                 "smsScan_AlphaT55_675_775","smsScan_AlphaT55_775_875","smsScan_AlphaT55_875"],hist = "m0_m12_mChi_noweight", Norm = None ,rebin= 2))



NoDeadEcal = threeToTwo(GetHist(File = r.TFile.Open("./results_NoDeadEcal_had_T2tt_100.root"),folder = ["smsScan_AlphaT55_375_475","smsScan_AlphaT55_475_575","smsScan_AlphaT55_575_675",
                 "smsScan_AlphaT55_675_775","smsScan_AlphaT55_775_875","smsScan_AlphaT55_875"],hist = "m0_m12_mChi_noweight", Norm = None ,rebin= 2))

NoLeptonVetos = threeToTwo(GetHist(File = r.TFile.Open("./results_GenLepPhotonVeto_NoRecoVetos_had_T2tt_100.root"),folder = ["smsScan_AlphaT55_375_475","smsScan_AlphaT55_475_575","smsScan_AlphaT55_575_675",
                 "smsScan_AlphaT55_675_775","smsScan_AlphaT55_775_875","smsScan_AlphaT55_875"],hist = "m0_m12_mChi_noweight", Norm = None ,rebin= 2))



AllCuts_No = threeToTwo(GetHist(File = r.TFile.Open("results_FullCutFLow_had_T2tt_100.root"),folder = ["smsScan_AlphaT55_375_475","smsScan_AlphaT55_475_575","smsScan_AlphaT55_575_675",
                 "smsScan_AlphaT55_675_775","smsScan_AlphaT55_775_875","smsScan_AlphaT55_875"],hist = "m0_m12_mChi_noweight", Norm = None ,rebin= 2))


NoMHTovMET = threeToTwo(GetHist(File = r.TFile.Open("./results_NoMHTovMET_had_T2tt_100.root"),folder = ["smsScan_AlphaT55_375_475","smsScan_AlphaT55_475_575","smsScan_AlphaT55_575_675",
                 "smsScan_AlphaT55_675_775","smsScan_AlphaT55_775_875","smsScan_AlphaT55_875"],hist = "m0_m12_mChi_noweight", Norm = None ,rebin= 2))




# Make the 4 efficiencies:
# AllCuts.Divide(nocuts)
AllCuts_GenLep.Draw("COLZ")
AllCuts_GenLep.SetTitle("AllCuts + GenLepVeto Efficiency")
c1.Print()
NoLeptonVetos.Draw("COLZ")
NoLeptonVetos.SetTitle("AllCuts - LeptonVetos + GenLepVeto Efficiency")
c1.Print()

# NoDeadEcal.Divide(nocuts)
NoDeadEcal.Draw("COLZ")
NoDeadEcal.SetTitle("AllCuts - DeadEcal")
c1.Print()
# NoLeptonVetos.Divide(nocuts)

# NoMHTovMET.Divide(nocuts)
NoMHTovMET.Draw("COLZ")
NoMHTovMET.SetTitle("AllCuts - MHTovMET")
c1.Print()


# ChangeDueToDeadEcal = AllCuts.Clone()
# ChangeDueToDeadEcal.Divide(NoDeadEcal)
# ChangeDueToDeadEcal.SetTitle("ChangeDueToDeadEcal")
# ChangeDueToDeadEcal.SetMinimum(0.8)
# ChangeDueToDeadEcal.Draw("COLZ")
# c1.Print()

ChangeDueToLeptonVetos = AllCuts_GenLep.Clone()
ChangeDueToLeptonVetos.Divide(NoLeptonVetos)
for bin in range(ChangeDueToLeptonVetos.GetNbinsX()*ChangeDueToLeptonVetos.GetNbinsY()):
  if ChangeDueToLeptonVetos.GetBinContent(bin) > 0.:  ChangeDueToLeptonVetos.SetBinContent(bin,1.0-ChangeDueToLeptonVetos.GetBinContent(bin))
ChangeDueToLeptonVetos.SetTitle("ChangeDueToLeptonVetos")
# ChangeDueToLeptonVetos.SetMinimum(5E-3)
# ChangeDueToLeptonVetos.SetMaximum(0.05)
ChangeDueToLeptonVetos.SetMaximum(0.2)
ChangeDueToLeptonVetos.SetMinimum(1E-2)
c1.canvas.SetLogz()
ChangeDueToLeptonVetos.Draw("COLZ")
c1.Print()



ChangeDueToMHTovMET = AllCuts_No.Clone()
ChangeDueToMHTovMET.Divide(NoMHTovMET)
for bin in range(ChangeDueToMHTovMET.GetNbinsX()*ChangeDueToMHTovMET.GetNbinsY()):
  if ChangeDueToMHTovMET.GetBinContent(bin) > 0.:  ChangeDueToMHTovMET.SetBinContent(bin,1.0-ChangeDueToMHTovMET.GetBinContent(bin))
  if ChangeDueToMHTovMET.GetBinContent(bin) < -1.:  ChangeDueToMHTovMET.SetBinContent(bin,0.)
ChangeDueToMHTovMET.SetTitle("ChangeDueToMHTovMET")
ChangeDueToMHTovMET.SetMaximum(0.2)
ChangeDueToMHTovMET.SetMinimum(1E-2)

# c1.canvas.SetLogz(False)
ChangeDueToMHTovMET.Draw("COLZ")
c1.Print()

ChangeDueToDeadEcal = AllCuts_No.Clone()
ChangeDueToDeadEcal.Divide(NoDeadEcal)
for bin in range(ChangeDueToDeadEcal.GetNbinsX()*ChangeDueToDeadEcal.GetNbinsY()):
  if ChangeDueToDeadEcal.GetBinContent(bin) > 0.:  ChangeDueToDeadEcal.SetBinContent(bin,1.0-ChangeDueToDeadEcal.GetBinContent(bin))
  if ChangeDueToDeadEcal.GetBinContent(bin) < -1.:  ChangeDueToDeadEcal.SetBinContent(bin,0.)
ChangeDueToDeadEcal.SetTitle("ChangeDueToDeadEcal")
ChangeDueToDeadEcal.SetMaximum(0.2)
ChangeDueToDeadEcal.SetMinimum(1E-2)
# c1.canvas.SetLogz(False)
ChangeDueToDeadEcal.Draw("COLZ")
c1.Print()




c1.close()
