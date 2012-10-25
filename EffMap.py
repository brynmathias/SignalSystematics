#!/usr/bin/env python

from plottingUtils import GetSumHist, Print, threeToTwo
import ROOT as r
models = [ 'T2' ] #'T2bb', 'T1bbbb', 'T1']


selection = [ 'had', 'muon' ][0]
histoDirPrefix = 'smsScan'
histoDirFmt = '{prefix}_{jb}_AlphaT55_{htRange}'
rootFileFmt = "rootFiles/sigScan_{model}_{sel}_2012_{jPt}_bt0.0_MChi-1.0.root"

bTagBins = [
        ('n_{b}=0', 'eq0b'),
        ('n_{b}=1', 'eq1b'),
        ('n_{b}=2', 'eq2b'),
        ('n_{b}=3', 'eq3b'),
        ('n_{b}#geq4', 'ge4b'),
        ]

jetBins = [
        ('n_{j}#leq3', 'le3j'),
        ('n_{j}#geq4', 'ge4j'),
        ]

j1 = [(', '.join([b[0], jetBins[0][0]]), '_'.join([b[1],jetBins[0][1]])) for b in bTagBins[:-1]]
j2 = [(', '.join([b[0], jetBins[1][0]]), '_'.join([b[1],jetBins[1][1]])) for b in bTagBins]
binCombinations = j1+j2

htBins = [ ('375'), ('275','325'), ('325','375'), ('375','475'), ('475','575'),
        ('575','675'), ('675','775'), ('775','875'), ('875',) ]

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

def setHistoOpts(h, model, title):
    h.SetTitle('{p}, {t}'.format(p=processStamps[model], t=title))
    h.GetXaxis().SetTitle(axes[model][0])
    h.GetYaxis().SetTitle(axes[model][1])
    h.GetZaxis().SetTitle('')
    h.GetXaxis().SetTitleOffset(h.GetXaxis().GetTitleOffset()*1.3)
    h.GetYaxis().SetTitleOffset(h.GetYaxis().GetTitleOffset()*1.3)
    h.GetXaxis().SetTitleSize(0.04)
    h.GetYaxis().SetTitleSize(0.04)
    h.GetXaxis().SetLabelSize(h.GetXaxis().GetLabelSize()*0.7)
    h.GetYaxis().SetLabelSize(h.GetYaxis().GetLabelSize()*0.7)
    h.GetZaxis().SetLabelSize(h.GetZaxis().GetLabelSize()*0.7)

def makeEffMap(model):
    print "------------------------------"
    print "     {m}".format(m=model)
    print "------------------------------"

    canvas = r.TCanvas()
    #c1 = Print("{model}_EffMap.pdf".format(model = model))
    canvas.SetLeftMargin(0.12)
    canvas.SetTopMargin(0.08)
    noCutFile = r.TFile.Open(rootFileFmt.format(model=model, jPt='100.0',
        sel=selection))
    noCuts = noCutFile.Get('smsScan_before/m0_m12_mChi_noweight')
    noCuts = threeToTwo(noCuts)
    noCuts.Draw("COLZ")

    for (title, selBin) in binCombinations:
        Cuts = None
        for htBin in htBins:
            oBin = "100.0"
            if htBin is ('275','325'): oBin = "73.7"
            if htBin is ('325','375'): oBin = "86.7"
            histoFile = r.TFile.Open(rootFileFmt.format(model=model, jPt=oBin,
                sel=selection))
            htBin_dir = histoDirFmt.format(prefix=histoDirPrefix,
                    jb=selBin, htRange='_'.join(htBin))
            h = histoFile.Get(htBin_dir+"/m0_m12_mChi_noweight")
            h = threeToTwo(h)
            if Cuts is None:
                Cuts = h.Clone()
            else:Cuts.Add(h)

            setHistoOpts(h,model,title)
            h.Divide(noCuts)
            h.Draw('COLZ')
            canvas.Print('pdfs/{m}_{b}_{s}_{h}.pdf'.format(m=model, b=selBin,
                s=selection, h='_'.join(htBin)))
        setHistoOpts(Cuts,model,title)
        Cuts.Divide(noCuts)
        Cuts.Draw('COLZ')
        canvas.Print('pdfs/{m}_{b}_Eff{s}Sum.pdf'.format(m=model, b=selBin,
            s=selection))
    noCutFile.Close()

if __name__ == '__main__':
    for model in models:
        makeEffMap(model)
