#!/usr/bin/env python
from plottingUtils import GetSumHist, Print, Legend
import ROOT as r



def MakeCumu(inHist):
    # eh =  [1.15, 1.36, 1.53, 1.73, 1.98, 2.21, 2.42, 2.61, 2.80, 3.00 ]
    # el =  [0.00, 1.00, 2.00, 2.14, 2.30, 2.49, 2.68, 2.86, 3.03, 3.19 ]
    cumulativeHist = inHist.Clone()
    maxbin = inHist.GetNbinsX()+1
    for bin in range(0,maxbin):
      err = r.Double(0)
      val = inHist.IntegralAndError(0, bin, err)
      # err = math.sqrt(val) if val > 9 else max(el[int(val)],eh[int(val)])
      cumulativeHist.SetBinContent(bin,val)
      cumulativeHist.SetBinError(bin,err)
    return cumulativeHist



dirs = ["OneMuon_375_475","OneMuon_475_575","OneMuon_575_675","OneMuon_675_775","OneMuon_775_875","OneMuon_875",]

def main():
  """docstring for main"""
  c1 = Print("MHToverMET_data_to_mc.pdf")
  c1.DoPageNum = False
  c1.open()
  Data = GetSumHist( File = ["Data.root",], Directories = dirs, Hist = "MHTovMET_all", Col = r.kBlack, Norm = None, LegendText = "Data")
  Data.HideOverFlow()
  # Data.hObj.Rebin(5)
  Mc = GetSumHist( File = ["WJets.root",], Directories = dirs, Hist = "MHTovMET_all", Col = r.kRed, Norm =None, LegendText = "MC")
  DY = GetSumHist( File = ["DY.root",], Directories = dirs, Hist = "MHTovMET_all", Col = r.kRed, Norm =None, LegendText = "MC")
  TTBar = GetSumHist( File = ["AllTop.root",], Directories = dirs, Hist = "MHTovMET_all", Col = r.kRed, Norm =None, LegendText = "MC")
  print "Data integral = %f"%(Data.hObj.Integral())
  print "WJets integral = %f"%(Mc.hObj.Integral())
  print "DY integral = %f"%(DY.hObj.Integral())
  print "TT integral = %f"%(TTBar.hObj.Integral())
  Mc.hObj.Add(DY.hObj)
  Mc.hObj.Add(TTBar.hObj)
  Mc.HideOverFlow()
  print "MC Total integral = %f"%(Mc.hObj.Integral())

  
  Data.hObj.Scale(1./Data.hObj.Integral())
  Mc.hObj.Scale(1./Mc.hObj.Integral())
  
  Data.hObj.Rebin(25)
  Mc.hObj.Rebin(25)


  Data.hObj.GetXaxis().SetRangeUser(0.,5.)

  # Mc.hObj.Rebin(5)

  DataCumu = r.TH1D(Data.hObj)
  McCumu = r.TH1D(Mc.hObj)
  DataCumu = MakeCumu(DataCumu)
  McCumu = MakeCumu(McCumu)
  # Mc.hObj.Scale(4.95)
  Data.Draw("h")
  Mc.Draw("sameh")
  leg = Legend(x1 = 0.6 , x2 = 0.8, y1 = 0.5 , y2 = 0.8 )
  leg.AddEntry(Data.hObj,Data.legendText,"lp")
  leg.AddEntry(Mc.hObj,Mc.legendText,"lp")
  leg.Draw("")
  c1.Print()
  
  DataCumu.Draw("h")
  McCumu.Draw("sameh")
  leg = Legend(x1 = 0.6 , x2 = 0.8, y1 = 0.5 , y2 = 0.8 )
  leg.AddEntry(Data.hObj,Data.legendText,"lp")
  leg.AddEntry(Mc.hObj,Mc.legendText,"lp")
  leg.Draw("")
  
  
  c1.Print()
  

  
  
  ratio = r.TH1D(Data.hObj)
  for bin in range(ratio.GetNbinsX()):
    if Mc.hObj.GetBinContent(bin) > 0 :print "BYHAND Bin lower edge = %f, value = %f"%(ratio.GetBinLowEdge(bin),Data.hObj.GetBinContent(bin)/Mc.hObj.GetBinContent(bin))

  ratio.Divide(Mc.hObj)
  ratio.Draw("h")
  for bin in range(ratio.GetNbinsX()):
    print "Bin lower edge = %f, value = %f pm %f"%(ratio.GetBinLowEdge(bin),ratio.GetBinContent(bin),ratio.GetBinError(bin))
  ratio.GetYaxis().SetRangeUser(0.,3.)
  ratio.GetYaxis().SetRangeUser(0.75,1.25)
  c1.Print()

  
  ratio = r.TH1D(DataCumu)
  ratio.Divide(McCumu)

  ratio.Draw("h")
  for bin in range(ratio.GetNbinsX()):
    print "Cumulative Bin lower edge = %f, value = %f pm %f"%(ratio.GetBinLowEdge(bin),ratio.GetBinContent(bin),ratio.GetBinError(bin))
  ratio.GetYaxis().SetRangeUser(0.,3.)
  ratio.GetYaxis().SetRangeUser(0.75,1.25)
  c1.Print()
  c1.close()
  
  
  
  pass
  







if __name__ == "__main__":
  main()
  
  
  
