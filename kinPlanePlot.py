#!/usr/bin/env python
import ROOT as r
from plottingstuff import *

r.gStyle.SetOptStat(0)
r.gROOT.SetBatch(r.kTRUE)

r.TH1.SetDefaultSumw2(1)
r.TH2.SetDefaultSumw2(1)

settings = {
    "model":["T2cc", "T1ttcc", "T2"][2],
    "mode":["JES", "ISR", "EffMap"][2],
    "inclHT":[False, True][0],
    "HTBins":["200_275", "275_325", "325_375", "375_475", "475_575", "575_675", "675_775", "775_875", "875_975", "975"],
    # "HTBins":["275_325", "325_375", "375_475", "475_575", "575_675", "675_775", "775_875", "875"],
    "deltaM":[False, True][0],
    "jMulti":["le3j", "ge4j", "eq2j", "eq3j", "ge2j"][1],
    "bMulti":["eq0b","eq1b","eq2b","ge0b"][3]
}

def threeToTwo(h3) :
    name = h3.GetName()
    binsz = h3.GetNbinsZ()
    # print binsz
    h2 = r.TH2D(name+"_2D",h3.GetTitle(),
                h3.GetNbinsX(), h3.GetXaxis().GetXmin(), h3.GetXaxis().GetXmax(),
                h3.GetNbinsY(), h3.GetYaxis().GetXmin(), h3.GetYaxis().GetXmax(),
                )
            
    # print h3.GetNbinsX(), h3.GetXaxis().GetXmin(), h3.GetXaxis().GetXmax()    
    for iX in range(1, 1+h3.GetNbinsX()) :
        for iY in range(1, 1+h3.GetNbinsY()) :
            content = h3.GetBinContent(iX, iY, 1) + h3.GetBinContent(iX, iY, 2)+ h3.GetBinContent(iX, iY, 0)
            h2.SetBinContent(iX, iY, content)
    h2.GetZaxis().SetTitle(h3.GetZaxis().GetTitle())
    return h2

def GetHist(File = None, folder = None, hist = None, Norm = None, rebinX = None, rebinY = None):
    h = None
    for f in folder:
        directory = File.Get(f)
        a = directory.Get(hist)
        if h is None:
            h = a.Clone()
        else: h.Add(a)

    return h

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
                dirs.append("smsScan_%s_%s_AlphaT0.51_%s" % (bM, jM, ht))

    return dirs

def getOutFile(model="", htbin="", format=""):

    bMulti = settings["bMulti"]
    jMulti = settings["jMulti"]

    ## convert individual string selections to lists for iteration
    if "list" not in str(type(bMulti)):
        bMulti = [bMulti]
    if "list" not in str(type(jMulti)):
        jMulti = [jMulti]

    if format == "txt":
        outName = "%s_%s_systOutput_%s_%s.txt" % (settings["mode"], model, "_".join(bMulti), "_".join(jMulti))
    elif "pdf" in format:
        outName = "%s_%s_%s_%s%s.pdf" % (settings["mode"], model, "_".join(bMulti), "_".join(jMulti),
                                            "_dM" if settings["deltaM"] else "")

    return "./out/"+outName

def printAxisSpec(hist=None):

	print "*** X-axis ***"
	print "NBins:", hist.GetNbinsX()
	print "Xmin:", hist.GetXaxis().GetXmin()
	print "Xmax:", hist.GetXaxis().GetXmax()

	print "\n*** Y-axis ***"
	print "NBins:", hist.GetNbinsY()
	print "Ymin:", hist.GetYaxis().GetXmin()
	print "Ymax:", hist.GetYaxis().GetXmax()
	print "\n"

model = settings["model"]

# centralRootFile100 = r.TFile.Open("./rootFiles/%s_v0/sigScan_%s_had_2012_100.0_bt0.0_MChi-1.0.root"%(model, model))
# centralRootFile87 = r.TFile.Open("./rootFiles/%s_v0/sigScan_%s_had_2012_86.7_bt0.0_MChi-1.0.root"%(model, model))
# centralRootFile73 = r.TFile.Open("./rootFiles/%s_v0/sigScan_%s_had_2012_73.7_bt0.0_MChi-1.0.root"%(model, model))

centralRootFile100 = r.TFile.Open("./rootFiles/test/sigScan_%s_had_2012_100.0_bt0.0_MChi-1.0.root"%(model))
centralRootFile87 = r.TFile.Open("./rootFiles/test/sigScan_%s_had_2012_86.7_bt0.0_MChi-1.0.root"%(model))
centralRootFile73 = r.TFile.Open("./rootFiles/test/sigScan_%s_had_2012_73.7_bt0.0_MChi-1.0.root"%(model))

c1 = r.TCanvas()

r.gPad.SetRightMargin(0.21)
r.gPad.SetLeftMargin(0.15)
r.gPad.SetTopMargin(0.08)
r.gPad.SetBottomMargin(0.15)

cutsNEvents = GetHist(File=centralRootFile73, folder=getRootDirs()[0:1], hist="m0_m12_noweight_tmp", Norm=None).Clone()
cutsNEvents.Add(GetHist(File=centralRootFile87, folder=getRootDirs()[1:2], hist="m0_m12_noweight_tmp", Norm=None))
cutsNEvents.Add(GetHist(File=centralRootFile100, folder=getRootDirs()[2:], hist="m0_m12_noweight_tmp", Norm=None))
# cutsNEvents = threeToTwo(cutsNEvents)

cutsHT = GetHist(File=centralRootFile73, folder=getRootDirs()[0:1], hist="m0_m12_noweight_HT", Norm=None).Clone()
cutsHT.Add(GetHist(File=centralRootFile87, folder=getRootDirs()[1:2], hist="m0_m12_noweight_HT", Norm=None))
cutsHT.Add(GetHist(File=centralRootFile100, folder=getRootDirs()[2:], hist="m0_m12_noweight_HT", Norm=None))

cutsMHT = GetHist(File=centralRootFile73, folder=getRootDirs()[0:1], hist="m0_m12_noweight_MHT", Norm=None).Clone()
cutsMHT.Add(GetHist(File=centralRootFile87, folder=getRootDirs()[1:2], hist="m0_m12_noweight_MHT", Norm=None))
cutsMHT.Add(GetHist(File=centralRootFile100, folder=getRootDirs()[2:], hist="m0_m12_noweight_MHT", Norm=None))

cutsalphaT = GetHist(File=centralRootFile73, folder=getRootDirs()[0:1], hist="m0_m12_noweight_alphaT", Norm=None).Clone()
cutsalphaT.Add(GetHist(File=centralRootFile87, folder=getRootDirs()[1:2], hist="m0_m12_noweight_alphaT", Norm=None))
cutsalphaT.Add(GetHist(File=centralRootFile100, folder=getRootDirs()[2:], hist="m0_m12_noweight_alphaT", Norm=None))

cutsjetM = GetHist(File=centralRootFile73, folder=getRootDirs()[0:1], hist="m0_m12_noweight_jetM", Norm=None).Clone()
cutsjetM.Add(GetHist(File=centralRootFile87, folder=getRootDirs()[1:2], hist="m0_m12_noweight_jetM", Norm=None))
cutsjetM.Add(GetHist(File=centralRootFile100, folder=getRootDirs()[2:], hist="m0_m12_noweight_jetM", Norm=None))

nX = cutsNEvents.GetNbinsX()
nY = cutsNEvents.GetNbinsY()


# print nX, nY
# print cutsHT.GetNbinsX(), cutsHT.GetNbinsY()
# print cutsNEvents.GetBinCenter(cutsNEvents.GetNbinsX()), cutsHT.GetBinCenter(cutsHT.GetNbinsX())
# print cutsNEvents.GetBinCenter(cutsNEvents.GetNbinsY()), cutsHT.GetBinCenter(cutsHT.GetNbinsY())
# print cutsNEvents.GetBinLowEdge(1), cutsHT.GetBinLowEdge(1)
# print cutsNEvents.GetXaxis()
# print type(cutsNEvents), type(cutsHT)

printAxisSpec(cutsNEvents)
printAxisSpec(cutsHT)

for x in range(nX):
	for y in range(nY):
		# print x, y
		# print cutsNEvents.GetBinContent(x,y), cutsHT.GetBinContent(x,y)
		pass

for i in range(nX*nY):
	if (cutsNEvents.GetBinCenter(i)-cutsHT.GetBinCenter(i) !=0):
		print "PANIX"

meanHT = cutsHT.Clone()
meanHT.Divide(cutsNEvents)

meanHT.Draw("colz")
c1.Print("test.pdf")