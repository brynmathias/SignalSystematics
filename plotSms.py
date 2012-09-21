#!/usr/bin/env python

from plottingUtils import *
import ROOT as r

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

c1 = Print("T2bb.pdf")
c1.open()

r.gPad.SetRightMargin(0.15)


dirs =[ "smsScan_before","smsScan_AlphaT55_275_325","smsScan_AlphaT55_325_375","smsScan_AlphaT55_375_475","smsScan_AlphaT55_475_575","smsScan_AlphaT55_575_675","smsScan_AlphaT55_675_775","smsScan_AlphaT55_775_875","smsScan_AlphaT55_875","smsScan_btag_>_0_AlphaT55_275_325","smsScan_btag_>_0_AlphaT55_325_375","smsScan_btag_>_0_AlphaT55_375_475","smsScan_btag_>_0_AlphaT55_475_575","smsScan_btag_>_0_AlphaT55_575_675","smsScan_btag_>_0_AlphaT55_675_775","smsScan_btag_>_0_AlphaT55_775_875","smsScan_btag_>_0_AlphaT55_875","smsScan_btag_==_0_AlphaT55_275_325","smsScan_btag_==_0_AlphaT55_325_375","smsScan_btag_==_0_AlphaT55_375_475","smsScan_btag_==_0_AlphaT55_475_575","smsScan_btag_==_0_AlphaT55_575_675","smsScan_btag_==_0_AlphaT55_675_775","smsScan_btag_==_0_AlphaT55_775_875","smsScan_btag_==_0_AlphaT55_875","smsScan_btag_==_1_AlphaT55_275_325","smsScan_btag_==_1_AlphaT55_325_375","smsScan_btag_==_1_AlphaT55_375_475","smsScan_btag_==_1_AlphaT55_475_575","smsScan_btag_==_1_AlphaT55_575_675","smsScan_btag_==_1_AlphaT55_675_775","smsScan_btag_==_1_AlphaT55_775_875","smsScan_btag_==_1_AlphaT55_875","smsScan_btag_==_2_AlphaT55_275_325","smsScan_btag_==_2_AlphaT55_325_375","smsScan_btag_==_2_AlphaT55_375_475","smsScan_btag_==_2_AlphaT55_475_575","smsScan_btag_==_2_AlphaT55_575_675","smsScan_btag_==_2_AlphaT55_675_775","smsScan_btag_==_2_AlphaT55_775_875","smsScan_btag_==_2_AlphaT55_875","smsScan_btag_>_2_AlphaT55_275_325","smsScan_btag_>_2_AlphaT55_325_375","smsScan_btag_>_2_AlphaT55_375_475","smsScan_btag_>_2_AlphaT55_475_575","smsScan_btag_>_2_AlphaT55_575_675","smsScan_btag_>_2_AlphaT55_675_775","smsScan_btag_>_2_AlphaT55_775_875","smsScan_btag_>_2_AlphaT55_875","smsScan_AlphaT52_53_275_325","smsScan_AlphaT52_53_325_375","smsScan_AlphaT52_53_375_475","smsScan_AlphaT52_53_475_575","smsScan_AlphaT52_53_575_675","smsScan_AlphaT52_53_675_775","smsScan_AlphaT52_53_775_875","smsScan_AlphaT52_53_875","smsScan_AlphaT53_55_275_325","smsScan_AlphaT53_55_325_375","smsScan_AlphaT53_55_375_475","smsScan_AlphaT53_55_475_575","smsScan_AlphaT53_55_575_675","smsScan_AlphaT53_55_675_775","smsScan_AlphaT53_55_775_875","smsScan_AlphaT53_55_875",]


for f in ["results_FIX_NEWSCAN_had_T2bb_100.root"]:
  for d in dirs:
    histo = GetSumHist(File = [f], Directories = [d], Hist = "m0_m12_mChi_noweight", Col = r.kBlack, Norm = None, LegendText = "No btag requirement")
    a = threeToTwo(histo.hObj)
    a.SetTitle(d)
    a.Draw("COLZ")
    c1.Print()

c1.close()