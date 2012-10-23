#!/usr/bin/env python

from plottingUtils import GetSumHist, Print, threeToTwo
import ROOT as r




def main():
    model = "T1bbbb"
    c1 = Print("{model}_EffMap.pdf".format(model = model))
    rFile = r.TFile.Open("rootFiles/sigScan_{model}_had_2012_100.0_bt0.0_MChi-1.0.root".format(model = model))
    JetBins = (

[
    "smsScan_ge0b_le3j_AlphaT55_275_325",
    "smsScan_ge0b_le3j_AlphaT55_325_375",
    "smsScan_ge0b_le3j_AlphaT55_375_475",
    "smsScan_ge0b_le3j_AlphaT55_475_575",
    "smsScan_ge0b_le3j_AlphaT55_575_675",
    "smsScan_ge0b_le3j_AlphaT55_675_775",
    "smsScan_ge0b_le3j_AlphaT55_775_875",
    "smsScan_ge0b_le3j_AlphaT55_875",
],
[
    "smsScan_ge0b_ge4j_AlphaT55_275_325",
    "smsScan_ge0b_ge4j_AlphaT55_325_375",
    "smsScan_eq0b_ge4j_AlphaT55_375_475",
    "smsScan_eq0b_ge4j_AlphaT55_475_575",
    "smsScan_eq0b_ge4j_AlphaT55_575_675",
    "smsScan_eq0b_ge4j_AlphaT55_675_775",
    "smsScan_eq0b_ge4j_AlphaT55_775_875",
    "smsScan_eq0b_ge4j_AlphaT55_875",
],


[
        "smsScan_ge0b_le3j_AlphaT55_275_325",
        "smsScan_ge0b_le3j_AlphaT55_325_375",
        "smsScan_ge0b_le3j_AlphaT55_375_475",
        "smsScan_ge0b_le3j_AlphaT55_475_575",
        "smsScan_ge0b_le3j_AlphaT55_575_675",
        "smsScan_ge0b_le3j_AlphaT55_675_775",
        "smsScan_ge0b_le3j_AlphaT55_775_875",
        "smsScan_ge0b_le3j_AlphaT55_875",
        "smsScan_ge0b_ge4j_AlphaT55_275_325",
        "smsScan_ge0b_ge4j_AlphaT55_325_375",
        "smsScan_eq0b_ge4j_AlphaT55_375_475",
        "smsScan_eq0b_ge4j_AlphaT55_475_575",
        "smsScan_eq0b_ge4j_AlphaT55_575_675",
        "smsScan_eq0b_ge4j_AlphaT55_675_775",
        "smsScan_eq0b_ge4j_AlphaT55_775_875",
        "smsScan_eq0b_ge4j_AlphaT55_875",

]


)
    noCuts = rFile.Get('smsScan_before/m0_m12_mChi_noweight')
    noCuts = threeToTwo(noCuts)
    noCuts.Draw("COLZ")
    # c1.Print()
    for jetBin in JetBins:
        Cuts = None
        for bin in jetBin:
            print bin
            oBin = "100.0"
            if "275_325" in bin: oBin = "73.7"
            if "325_375" in bin: oBin = "86.7"
            File = r.TFile.Open("rootFiles/sigScan_{model}_had_2012_{bin}_bt0.0_MChi-1.0.root".format(bin = oBin, model = model))
            h = File.Get(bin+"/m0_m12_mChi_noweight")
            h = threeToTwo(h)
            # h.Draw("COLZ")
            # c1.Print()
            if Cuts is None:
                Cuts = h.Clone()
            else:Cuts.Add(h) 
            # Cuts.Draw("COLZ")
            # c1.Print()
        Cuts.Draw("COLZ")
        Cuts.SetTitle(jetBin[0])
        # c1.Print()
        Cuts.Divide(noCuts)
        # Cuts = threeToTwo(Cuts)
        Cuts.Draw("COLZ")
        c1.Print()
    c1.close()
    rFile.Close()


















    
if __name__ == '__main__':
    main()