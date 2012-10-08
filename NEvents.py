#!/usr/bin/env python
import ROOT as r
someFile = "rootFiles/sigScan_T2bb_muon_2012_100.0_bt0.0_MChi-1.root"
rFile = r.TFile.Open(someFile)
for k in rFile.GetListOfKeys():
  h = rFile.Get(k.GetName()+"/m0_m12_mChi_noweight")
  if h.Integral() > 0:
    print k.GetName() , h.Integral()
