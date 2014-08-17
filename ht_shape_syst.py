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
	"model":    ["T2cc", "T2", "T2_4body", "T2cc_250_240", "T2cc_250_170", "T2bw_0p75"][0],
	"version":  25,
	"mode":		["JES", "ISR"][1:],
	"HTBins":   ["200_275", "275_325", "325_375", "375_475", "475_575", "575_675", "675_775", "775_875", "875_975", "975_1075", "1075"],
	"jMulti":   ["le3j", "ge4j", "ge2j"],
	"bMulti":   ["eq0b", "eq1b", "eq2b", "ge0b"],
	"Logy":		[False, True][1]
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


def getRootDirs(bMulti_="", jMulti_="", sitv_=False, thresh_ = ""):

	dirs = []


	for ht in settings["HTBins"]:
		thisHT = ht.split("_")
		if len(thisHT)<2: thisHT.append(None)
		aT = alphaTDict["%s,%s" % (thisHT[0], thisHT[1])][0]
		aTString = "AlphaT%s" % (str(aT[0])+"_"+str(aT[1]) if aT[1] else str(aT[0]))
		dirs.append("smsScan_%s_%s_%s%s_%s" % (bMulti_, jMulti_,
					"SITV_" if sitv_ else "", aTString, ht))

	if thresh_ == "73.7":
		return dirs[0:2]
	elif thresh_ == "86.7":
		return dirs[2:3]
	elif thresh_ == "100.0":
		return dirs[3:]
	else:
		return dirs

###-------------------------------------------------------------------###

def getOutFile(format="", jMulti_ = "", bMulti_ = "", mode_ = ""):

	if format == "txt":
		outName = "%s_%s_HTshape_systOutput_%s_%s.txt" % (mode_, settings['model'], bMulti_, jMulti_)
	elif "pdf" in format:
		outName = "%s_HTshape_%s_v%d_%s_%s.pdf" % (mode_, settings['model'], settings['version'], bMulti_,
												jMulti_)

	return "./out/"+outName

###-------------------------------------------------------------------###

def clone_hist_content(source = None, target = None):

	if not source:
		print "func clone_hist_content: No source histo specified."
	if not target:
		print "func clone_hist_content: No target histo specified."

	last_bin = 1
	target_bin_val = 0
	target_bin_err = 0
	for i in range(source.GetNbinsX()):
		bin = target.FindBin( source.GetBinCenter(i) )
		bin_val = source.GetBinContent(i)
		bin_err = source.GetBinError(i)
		
		if last_bin == bin:
			# if still in same bin, add val, add err in quadrature
			target_bin_val += bin_val
			target_bin_err += math.pow(bin_err, 2)
		else:
			# new bin, so fill the previous bin with the summed val and err
			target.SetBinContent( last_bin, target_bin_val )
			target.SetBinError(last_bin, math.sqrt(target_bin_err))
			target_bin_val = 0 #reset val calc
			target_bin_err = 0 #reset error calc
		last_bin = bin

	return target

###-------------------------------------------------------------------###

def make_ht_shape_plots(bMulti = "", jMulti = "", mode = ""):

	print "\n>>> Running shape systs for:", mode, jMulti, bMulti 

	htbins_  = [200., 275., 325.] + [375.+100.*i for i in range(8)] + [1200.]
	variations = ['central', 'plus', 'minus']
	# c1 = r.TCanvas("shape_syst", "shape_syst", 2000, 1500)
	c1 = r.TCanvas()

	# set some global colours for plots
	v_colours = dict.fromkeys(variations)
	v_colours['central'] = r.kBlack
	v_colours['plus'] = r.kBlue
	v_colours['minus'] = r.kOrange+7

	if mode == "ISR":
		ins = "isr_"
	else:
		ins = ""

	# create empty histo dicts
	in_plots = dict.fromkeys(variations, None)
	ana_plots = dict.fromkeys(variations, None)

	# add a noweight entry (for mc stat err calculation)
	in_plots['noweight'] = None
	ana_plots['noweight'] = None

	# open rootfiles and extract basic HT plots
	for thresh in ['73.7', '86.7', '100.0']:
		f_tmp_central 	= r.TFile.Open("./rootFiles/%s_v%d/sigScan_%s_had_2012_%s_bt0.0_MChi-1.0.root" %
									(settings["model"], settings["version"], settings["model"], thresh))
		f_tmp_plus		= r.TFile.Open("./rootFiles/%s_v%d/%s/sigScan_%s_had_2012_%s_%s+ve_bt0.0_MChi-1.0.root" %
									(settings["model"], settings["version"], mode, settings["model"], thresh, ins))
		f_tmp_minus		= r.TFile.Open("./rootFiles/%s_v%d/%s/sigScan_%s_had_2012_%s_%s-ve_bt0.0_MChi-1.0.root" %
									(settings["model"], settings["version"], mode, settings["model"], thresh, ins))
		# f_tmp_plus = f_tmp_central
		# f_tmp_minus = f_tmp_central

		# get the correct root directories according to the threshold
		dirs = getRootDirs(bMulti_ = bMulti, jMulti_ = jMulti, sitv_ = True, thresh_ = thresh)

		for vari in variations:
			if not in_plots[vari]:
				in_plots[vari] = GetHist(File=eval("f_tmp_%s"%vari), folder=dirs, hist="commonHT_weight", Norm=None).Clone()
			else:
				in_plots[vari].Add(GetHist(File=eval("f_tmp_%s"%vari), folder=dirs, hist="commonHT_weight", Norm=None))

		if not in_plots['noweight']:
			in_plots['noweight'] = GetHist(File=eval("f_tmp_central"), folder=dirs, hist="commonHT_noweight", Norm=None).Clone()
		else:
			in_plots['noweight'].Add(GetHist(File=eval("f_tmp_central"), folder=dirs, hist="commonHT_noweight", Norm=None))

		# clean up tmp TFiles
		del f_tmp_central, f_tmp_plus, f_tmp_minus

	# create some empty plots with ana binning
	for v in variations+['noweight']:
		h_tmp_ana_ht = r.TH1D("ana_ht_%s"%v, "ana_ht_%s"%v, len(htbins_)-1, array('d', htbins_))
		ana_plots[v] = clone_hist_content(source = in_plots[v], target = h_tmp_ana_ht)
		del h_tmp_ana_ht
	
	# get n events per plot
	nevents_central = ana_plots['central'].Integral()
	nevents_plus 	= ana_plots['plus'].Integral()
	nevents_minus 	= ana_plots['minus'].Integral()


	######################
	# draw up some plots #
	######################

	r.gStyle.SetOptStat(0)

	lg = r.TLegend(0.75, 0.65, 0.9, 0.85)
	lg.SetFillColor(0)
	lg.SetFillStyle(0)
	lg.SetLineColor(0)
	lg.SetLineStyle(0)

	pd1 = r.TPad("pd1", "pd1", 0., 0.3, 1., 1.)
	pd1.SetBottomMargin(0.005)
	pd1.SetRightMargin(0.05)
	pd1.Draw()

	pd2 = r.TPad("pd2", "pd2", 0., 0.02, 1., 0.3)
	pd2.SetTopMargin(0.05)
	pd2.SetBottomMargin(0.26)
	pd2.SetRightMargin(0.05)
	pd2.SetGridx(1)
	pd2.SetGridy(1)
	pd2.Draw()

	pd1.cd()

	for n, v in enumerate(variations):
		if n==0:
			ana_plots[v].Draw("hist")
			ana_plots[v].GetXaxis().SetRangeUser(200., 1200.)
			ana_plots[v].GetXaxis().SetTitle("CommonHT (GeV)")
			ana_plots[v].GetYaxis().SetTitle("# count")
			ana_plots[v].SetTitle("HT %s Systematic Variations" % mode)
		else:
			ana_plots[v].Draw("histsame")
		ana_plots[v].SetLineColor( v_colours[v] )
		ana_plots[v].SetLineWidth(2)
		if 'central' not in v:
			ana_plots[v].SetLineStyle(2)
			# print "Before:", ana_plots[v].GetBinError(3)
			ana_plots[v].Scale( float(nevents_central/eval("nevents_%s"%v)) )
			print "> %s normalisation: %.3f" % (v, float(nevents_central/eval("nevents_%s"%v)))
	
		lg.AddEntry(ana_plots[v], v, "L")
	lg.Draw()

	# for i in range(ana_plots['central'].GetNbinsX()):
	# 	val = ana_plots['central'].GetBinContent(i)
	# 	# err = math.sqrt(ana_plots['central'].GetBinContent(i))
	# 	err = ana_plots['central'].GetBinError(i)
		
	# 	try:
	# 		up_frac_err = (val+err)/val
	# 	except ZeroDivisionError:
	# 		up_frac_err = 0
		
	# 	print "%f\t\t%.3f\t\t%.3f\t\t%.3f" % (ana_plots['central'].GetBinLowEdge(i), val, err, up_frac_err)

	ana_plots['central'].Draw('histsame') #cheeky hack so central line is on top (probably a better way...)

	ana_plots['central'].SetTitleSize(.12)
	ana_plots['central'].SetTitleSize(.05, "Y")
	ana_plots['central'].SetTitleOffset(.8, "Y")

	if settings["Logy"]:
		pd1.SetLogy(1)

	pd2.cd()

	hrat_central = ana_plots['central'].Clone()
	hrat_central.Divide(ana_plots['central'])
	# hrat_central.Scale(1./math.sqrt(2))
	# hrat_central.SetFillColor(15)
	hrat_central.SetFillColor(16)
	# hrat_central.SetFillStyle(3001)
	hrat_central.SetLabelSize(0.12, "X")
	hrat_central.SetLabelSize(0.07, "Y")
	hrat_central.SetTitleSize(0.13, "X")
	hrat_central.SetTitleSize(0.11, "Y")
	hrat_central.SetTitleOffset(0.25, "Y")
	hrat_central.SetTitleOffset(.9, "X")
	
	# scale all errors by sqrt(2) due to division by itself
	for i in range(hrat_central.GetNbinsX()):
		bin_err = hrat_central.GetBinError(i+1)
		hrat_central.SetBinError(i+1, bin_err/math.sqrt(2))

	# # calculate the mc stat uncert from unweighted distribution
	# for i in range(hrat_central.GetNbinsX()):
	# 	original_bin_val = ana_plots['noweight'].GetBinContent(i+1)
	# 	try:
	# 		new_bin_err = math.sqrt(original_bin_val)/original_bin_val
	# 	except ZeroDivisionError:
	# 		new_bin_err = 0.
	# 	hrat_central.SetBinError(i+1, new_bin_err)

	hrat_central.Draw("E5")

	hrat_plus = ana_plots['central'].Clone()
	hrat_plus.Divide(ana_plots['plus'])
	hrat_plus.SetMarkerStyle(4)
	hrat_plus.SetMarkerSize(.7)
	hrat_plus.SetLineWidth(1)
	hrat_plus.SetMarkerColor(v_colours['plus'])
	hrat_plus.SetLineColor(v_colours['plus'])
	hrat_plus.Draw("pe0same")

	# fit_plus = r.TF1("fit_plus","pol1", 200., 1200.)
	# hrat_plus.Fit(fit_plus, "R+")   

	hrat_minus = ana_plots['central'].Clone()
	hrat_minus.Divide(ana_plots['minus'])
	hrat_minus.SetMarkerStyle(4)
	hrat_minus.SetMarkerSize(.7)
	hrat_minus.SetLineWidth(1)
	hrat_minus.SetMarkerColor(v_colours['minus'])
	hrat_minus.SetLineColor(v_colours['minus'])
	hrat_minus.Draw("pe0same")

	# fit_minus = r.TF1("fit_minus","pol1", 200., 1200.)
	# hrat_minus.Fit(fit_minus, "R+")

	hrat_central.GetYaxis().SetTitle("Ratio (Nom/Var)")
	hrat_central.SetTitleOffset(.45, "Y")
	hrat_central.SetTitleSize(.08, "Y")
	hrat_central.GetYaxis().SetRangeUser(0.5, 1.5)
	hrat_central.SetTitle("")

	c1.Print( getOutFile(format="pdf", jMulti_ = jMulti, bMulti_ = bMulti, mode_ = mode) )


	#########################
	# make some text output #
	#########################

	out_text = ""
	out_text += "*"*40
	out_text += "\nHT bin\t\tPlus\t\tMinus\n\n"
	
	for i in range(hrat_plus.GetNbinsX()):
		out_text += "%d-%d\t" % (hrat_plus.GetBinLowEdge(i+1),
									hrat_plus.GetBinLowEdge(i+1) + hrat_plus.GetBinWidth(i+1))
		if hrat_plus.GetBinLowEdge(i+1)<900:
			out_text += "\t"
		out_text += "%.3f\t\t" % hrat_plus.GetBinContent(i+1)
		out_text += "%.3f\n" % hrat_minus.GetBinContent(i+1)
	out_text += "\n"

	text_file = file( getOutFile(format="txt", jMulti_ = jMulti, bMulti_ = bMulti, mode_ = mode), 'w' )
	text_file.write(out_text)
	text_file.close()

if __name__ == "__main__":

	ensure_list(var = settings["jMulti"], tag = "jMulti settings")
	ensure_list(var = settings["bMulti"], tag = "bMulti settings")
	ensure_list(var = settings["mode"], tag = "mode settings")

	for opts in itertools.product(settings['bMulti'], settings['jMulti'], settings['mode']):
		make_ht_shape_plots(bMulti = opts[0], jMulti = opts[1], mode = opts[2])
