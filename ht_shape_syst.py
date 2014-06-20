#!/usr/bin/env python

import ROOT as r
from plottingUtils import Print, MakeCumu, SetBatch
from plotDetails import mapRanges, mapDMRanges, alphaTDict
from array import array
import math
import itertools


###-------------------------------------------------------------------###
###-------------------------------------------------------------------###

r.gROOT.SetBatch(r.kTRUE)
r.gStyle.SetOptStat(0)
r.TH1.SetDefaultSumw2(1)
r.TH2.SetDefaultSumw2(1)


settings = {
	"model":    ["T2cc", "T2", "T2_4body"][0],
	"version":  25,
	"mode":		["JES", "ISR"][0:1],
	"HTBins":   ["200_275", "275_325", "325_375", "375_475", "475_575", "575_675", "675_775", "775_875", "875_975", "975_1075", "1075"],
	"jMulti":   ["le3j", "ge4j", "ge2j"][-1:],
	"bMulti":   ["eq0b", "eq1b", "eq2b", "ge0b"][-1:],
}

###-------------------------------------------------------------------###


def GetHist(File=None, folder=None, hist=None, Norm=None, rebinX=None, rebinY=None):
	h = None

	for f in folder:
		directory = File.Get(f)
		a = directory.Get(hist)
		if h is None:
			h = a.Clone()
		else:
			h.Add(a)

	if rebinX:
		h.RebinX(rebinX)
	if rebinY:
		h.RebinY(rebinY)

	return h

###-------------------------------------------------------------------###


def ensure_list(var = None, tag = ""):

	if not var:
		return

	if type(var) != list:
		print "%s must be list." % tag
		exit()

	

###-------------------------------------------------------------------###


def getRootDirs(bMulti_="", jMulti_="", sitv_=False):

	dirs = []


	for ht in settings["HTBins"]:
		thisHT = ht.split("_")
		if len(thisHT)<2: thisHT.append(None)
		aT = alphaTDict["%s,%s" % (thisHT[0], thisHT[1])][0]
		aTString = "AlphaT%s" % (str(aT[0])+"_"+str(aT[1]) if aT[1] else str(aT[0]))
		dirs.append("smsScan_%s_%s_%s%s_%s" % (bMulti_, jMulti_,
					"SITV_" if sitv_ else "", aTString, ht))

	return dirs

###-------------------------------------------------------------------###

def getOutFile(model="", htbin="", format="", jMulti_ = "", bMulti_ = "", mode_ = ""):

	if format == "txt":
		outName = "%s_%s_systOutput_%s_%s.txt" % (mode_, model, bMulti_, jMulti_)
	elif "pdf" in format:
		outName = "%s_%s_v%d_%s_%s_%s%s.pdf" % (mode_, model, settings['version'], bMulti_, jMulti_,
											htbin, "_dM" if settings["deltaM"] else "")

	return "./out/"+outName

###-------------------------------------------------------------------###

def make_ht_shape_plots(bMulti = "", jMulti = "", mode = ""):


	if mode == "ISR":
		ins = "isr_"
	else:
		ins = ""

	in_files = {
				'central':	{},
				'plus':		{},
				'minus':	{}
				}

	for thresh in ['73.7', '86.7', '100.0']:
		in_files['central'][thresh] = r.TFile.Open("./rootFiles/%s_v%d/sigScan_%s_had_2012_%s_bt0.0_MChi-1.0.root" %
									(settings["model"], settings["version"], settings["model"], thresh))
		in_files['plus'][thresh]	= r.TFile.Open("./rootFiles/%s_v%d/%s/sigScan_%s_had_2012_%s_%s+ve_bt0.0_MChi-1.0.root" %
									(settings["model"], settings["version"], mode, settings["model"], thresh, ins))
		in_files['minus'][thresh]	= r.TFile.Open("./rootFiles/%s_v%d/%s/sigScan_%s_had_2012_%s_%s-ve_bt0.0_MChi-1.0.root" %
									(settings["model"], settings["version"], mode, settings["model"], thresh, ins))

	for i in in_files:
		print i
		print "\t", in_files[i]



if __name__ == "__main__":

	ensure_list(var = settings["jMulti"], tag = "jMulti settings")
	ensure_list(var = settings["bMulti"], tag = "bMulti settings")
	ensure_list(var = settings["mode"], tag = "mode settings")

	for opts in itertools.product(settings['bMulti'], settings['jMulti'], settings['mode']):
		print opts
		make_ht_shape_plots(bMulti = opts[0], jMulti = opts[1], mode = opts[2])
