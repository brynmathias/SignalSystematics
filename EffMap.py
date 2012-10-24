#!/usr/bin/env python

from plottingUtils import GetSumHist, Print, threeToTwo
import ROOT as r
models = [ 'T2' ] #'T2bb', 'T1bbbb', 'T1']

histoDirPrefix = 'smsScan'
histoDirFmt = '{prefix}_{sel}_AlphaT55_{htRange}'
rootFileFmt = "rootFiles/sigScan_{model}_had_2012_{jPt}_bt0.0_MChi-1.0.root"

processStamps =  {
    'T2'     : "pp #rightarrow #tilde{q} #tilde{q}, "
        "#tilde{q} #rightarrow q + LSP, m(#tilde{g})>>m(#tilde{q})",
    'T2bb'   : "pp #rightarrow #tilde{b} #tilde{b}, "
        "#tilde{b} #rightarrow b + LSP, m(#tilde{g})>>m(#tilde{b})",
    'T2tt'   : "pp #rightarrow #tilde{t} #tilde{t}, "
        "#tilde{t} #rightarrow t + LSP, m(#tilde{g})>>m(#tilde{t})",
    'T1'     : "pp #rightarrow #tilde{g} #tilde{g}, "
        "#tilde{g} #rightarrow 2q + LSP, m(#tilde{q})>>m(#tilde{g})",
    'T1bbbb' : "pp #rightarrow #tilde{g} #tilde{g}, "
        "#tilde{g} #rightarrow 2b + LSP, m(#tilde{b})>>m(#tilde{g})",
    'T1tttt' : "pp #rightarrow #tilde{g} #tilde{g}, "
        "#tilde{g} #rightarrow 2t + LSP, m(#tilde{t})>>m(#tilde{g})",
    }

axes = {
        'T1': ('m_{gluino} (GeV)', 'm_{LSP} (GeV)'),
        'T2': ('m_{squark} (GeV)', 'm_{LSP} (GeV)'),
        'T1bbbb': ('m_{gluino} (GeV)', 'm_{LSP} (GeV)'),
        'T2bb': ('m_{sbottom} (GeV)', 'm_{LSP} (GeV)'),
    }

htBins = [ ('275','325'), ('325','375'), ('375','475'), ('475','575'),
        ('575','675'), ('675','775'), ('775','875'), ('875',) ]

JetBins = {
    'T2': {
        'n_{b}=1, n_{j}#geq4': ('eq1b','ge4j'),
        'n_{b}=2, n_{j}#geq4': ('eq2b','ge4j'),
        },
    'T2bb': {
        'n_{b}=1, n_{j}#leq3': ('eq1b', 'le3j'),
        },
    'T1': {
        'n_{b}=0, n_{j}#geq4': ('eq0b', 'ge4j'),
        'n_{b}=1, n_{j}#geq4': ('eq1b', 'ge4j'),
        },
    'T1bbbb': {
        'n_{b}=2, n_{j}#geq4': ('eq2b', 'ge4j'),
        'n_{b}=3, n_{j}#geq4': ('eq3b', 'ge4j'),
        'n_{b}=4, n_{j}#geq4': ('eq4b', 'ge4j'),
        },
    }

def makeEffMap(model):
    print "------------------------------"
    print "     {m}".format(m=model)
    print "------------------------------"
    c1 = Print("{model}_EffMap.pdf".format(model = model))
    c1.canvas.SetLeftMargin(0.12)
    c1.canvas.SetTopMargin(0.08)
    noCutFile = r.TFile.Open(rootFileFmt.format(model=model, jPt='100.0'))
    noCuts = noCutFile.Get('smsScan_before/m0_m12_mChi_noweight')
    noCuts = threeToTwo(noCuts)
    noCuts.Draw("COLZ")
    # c1.Print()
    for title, jetBin in JetBins[model].iteritems():
        Cuts = None
        for htBin in htBins:
            oBin = "100.0"
            if htBin is (275,325): oBin = "73.7"
            if htBin is (325,375): oBin = "86.7"
            histoFile = r.TFile.Open(rootFileFmt.format(model=model, jPt=oBin))
            htBin_dir = histoDirFmt.format(prefix=histoDirPrefix,
                    sel='_'.join(jetBin), htRange='_'.join(htBin))
            h = histoFile.Get(htBin_dir+"/m0_m12_mChi_noweight")
            h = threeToTwo(h)
            if Cuts is None:
                Cuts = h.Clone()
            else:Cuts.Add(h)

        Cuts.SetTitle('{p}, {t}'.format(p=processStamps[model], t=title))
        Cuts.GetXaxis().SetTitle(axes[model][0])
        Cuts.GetYaxis().SetTitle(axes[model][1])
        Cuts.GetZaxis().SetTitle('')
        Cuts.GetXaxis().SetTitleOffset(Cuts.GetXaxis().GetTitleOffset()*1.3)
        Cuts.GetYaxis().SetTitleOffset(Cuts.GetYaxis().GetTitleOffset()*1.3)
        Cuts.GetXaxis().SetTitleSize(0.04)
        Cuts.GetYaxis().SetTitleSize(0.04)
        Cuts.GetXaxis().SetLabelSize(Cuts.GetXaxis().GetLabelSize()*0.7)
        Cuts.GetYaxis().SetLabelSize(Cuts.GetYaxis().GetLabelSize()*0.7)
        Cuts.GetZaxis().SetLabelSize(Cuts.GetZaxis().GetLabelSize()*0.7)
        Cuts.Divide(noCuts)
        Cuts.Draw("COLZ")
        c1.Print()
    c1.close()
    noCutFile.Close()

if __name__ == '__main__':
    for model in models:
        makeEffMap(model)
