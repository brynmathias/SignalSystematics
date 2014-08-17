# import plottingUtils as putils
import signalUtils as sutils
import numpy as np
import ROOT as r
import math as ma

###---------------------------------------------------------------------------------------------###
###---------------------------------------------------------------------------------------------###
# set some global root variables for the instance imported above

r.gROOT.SetBatch(r.kTRUE)
r.gStyle.SetOptStat(0)
r.TH1.SetDefaultSumw2(1)
r.TH2.SetDefaultSumw2(1)

###---------------------------------------------------------------------------------------------###
###---------------------------------------------------------------------------------------------###

class effMap(object):
	'''container for efficiency map'''
	def __init__(self, nom = None, denom = None, nom_err = None, denom_err = None):
		self._hist = None
		self._errorHist = None
		self._nom = nom
		self._nomerr = nom_err
		self._denom = denom
		self._denomerr = denom_err
		self._mean = 0.
		self._max = 0.
		self._min = 1.
		self._rms = 0.

		if not denom:
			self._hist = nom

		self.process()

		###
		# 1. Add check that both hists are self-consistent

	def process(self):
		'''make eff hist and calculate values'''

		if self._denom:
			self._hist = self._nom.Clone()
			self._hist.Divide(self._denom)

		self._nBins = self._hist.GetNbinsX() * self._hist.GetNbinsY() + 1000
		vals = []
		for i in range(1, self._nBins):
			val = self._hist.GetBinContent(i)
			if val > 0:
				vals.append(val)
		
		val_array = np.array(vals)
		
		if len(vals):
			self._rms 	= np.var(val_array)
			self._mean 	= np.mean(val_array)
			self._min 	= np.amin(val_array)
			self._max 	= np.amax(val_array)
		else:
			self._rms 	= 0.
			self._mean 	= 0.
			self._min 	= 0.
			self._max 	= 0.
		# create a histogram of eff errors
		self.makeErrorHist()

	def makeErrorHist(self):
		'''caclulate bin by bin stat error'''

		# clone a hist to fill with errors
		self._errorHist = self._nom.Clone()

		# loop all bins
		for i in range(1, self._nBins):

			# skip points that aren't in the scan
			# note this skips points with no num, a denom, and assoc errors...
			if self._nom.GetBinContent(i) == 0:
				self._errorHist.SetBinContent(i, 0.)
				continue

			# if no existing error hists are passed, calc from scratch
			if not self._nomerr and not self._denomerr:
				# if denom hist, sum stat errors
				if self._denom:
					nom_val 	= self._nom.GetBinContent(i)
					denom_val 	= self._denom.GetBinContent(i)

					if nom_val != 0. and denom_val != 0.:
						
						# calculate the statistical error for each yield
						nom_err 	= 1./ma.sqrt(nom_val)
						denom_err	= 1./ma.sqrt(denom_val)

						nom_ratio = nom_err/nom_val
						denom_ratio = denom_err/denom_val
					
						eff_val = nom_val/denom_val
						eff_err = eff_val*ma.sqrt(nom_ratio+denom_ratio)
					
					else:
						# set to zero if either hist has no entries
						eff_err = 0.
				else:
					# if no denom hist, just use error cald'd from nom hist
					eff_val = self._nom.GetBinContent(i)
					if eff_val > 0.:
						eff_err = 1./ma.sqrt(eff_val)
					else:
						eff_err = 0.
					
			else:
				# if error hists are passed, combine
				nom_err 	= self._nomerr.GetBinContent(i)
				denom_err 	= self._denomerr.GetBinContent(i)
				eff_err 	= ma.sqrt(ma.pow(nom_err, 2) + ma.pow(denom_err, 2))

				if self._denom:
					if self._denom.GetBinContent(i) > 0.:
						eff_val 	= self._nom.GetBinContent(i)/self._denom.GetBinContent(i)
				else:
					eff_val = 0.

			self._errorHist.SetBinContent(i, eff_err)
			# if eff_err and eff_val:
			# 	pass
				# print eff_err/eff_val
				# if float(eff_err/eff_val)>1.:
				# 	xbin, ybin, zbin = r.Long(0.), r.Long(0.), r.Long(0.)
				# 	self._nom.GetBinXYZ(i, xbin, ybin, zbin)
				# 	print "mStop: %f, mLSP: %f" % (self._nom.GetXaxis().GetBinCenter(xbin), self._nom.GetYaxis().GetBinCenter(ybin))
				# 	print "Nom_val: %f, Nom_err: %f, denom_val: %f, denom_err: %f" % (nom_val, nom_err, denom_val, denom_err)
			# if nom_val and denom_val:
			# 	print "Nom_val: %f, Nom_err: %f, denom_val: %f, denom_err: %f" % (nom_val, nom_err, denom_val, denom_err)

	def shiftCentre(self, shift = 0.):
		'''Shift all values by 1, to centre around zero'''
		self._mean 	+= shift
		self._min 	+= shift
		self._max 	+= shift

		for i in range(1, self._nBins+1):
			val = self._hist.GetBinContent(i)
			if val > 0.:
				self._hist.SetBinContent(i, val+shift)
			else:
				# set null points to -1000.
				self._hist.SetBinContent(i, -666.)

	def invertHist(self):
		''' invert an effMap object and all attributes'''
		self._mean 	= 1./self._mean
		self._min 	= 1./self._min
		self._max 	= 1./self._max

		for i in range(1, self._hist.GetNbinsX()*
						self._hist.GetNbinsY()+1000):
			val = self._hist.GetBinContent(i)
			if val > 0.:
				self._hist.SetBinContent(i, 1./val)


###---------------------------------------------------------------------------------------------###

class systMap(object):
	'''Simple systematic plotting'''
	
	### TO-DO ###
	# 1. implement deltaM plotting
	# 2. top-level pageNum plotting switch?
	# 3. text plotting (add plotString variable) - DONE
	# 4. Add stats print out to each plot - DONE
	# 5. Cap values in each plot (maxVal from run_syst_mapping.py?)
	# 6. implement logging module
	# 7. invert efficiency change for cut_systs
	# 8. cut systs still plotting down (and 3jet syst) - DONE

	def __init__(self, up = None, down = None, central = None, nocuts = None, test = "",
					model = ""):

		self._yieldPlots = {
							"up":		up,
							"down":		down,
							"central":	central,
							"nocuts":	nocuts
							}
		self._effs = {
						"up":			None,
						"down":			None,
						"central":		None,
						"up_change": 	None,
						"down_change": 	None,
						}
		self._syst = None
		self._model = model
		self._test = test
		self._cutSyst = True if self._test in ["MHT_MET", "DeadECAL", "LeptonVeto"][:2] else False
		self._ranges = {
					'x':[],
					'y':[],
					'z':[],
		}
		self._logZ = False
		self._plotString = "colz"
		self._nBins = self._yieldPlots['central'].GetNbinsX()* \
						self._yieldPlots['central'].GetNbinsY() + 1000

		# get model specific plot details
		import plotDetails as pdets
		try:
			self._plotSpec = pdets.modelPlotDetails[self._model]
		except KeyError:
			# print ">>> Warning: systMap: Model details not found. Using default values."
			self._plotSpec = {
					'xRange': [0., 1000.],
					'yRange': [0., 1000.],
					'xDMRange': [0., 1000.],
					'yDMRange': [0., 1000.],
					'xTitle': "m_{mother} (GeV)",
					'yTitle': "m_{LSP} (GeV)",
			}

		# get z ranges (varies for each test, not model)
		try:
			self._plotSpec['zRange'] = pdets.systZRanges[self._test]
		except KeyError:
			# print ">>> Warning: systMap: Test zRange not found. Using default values."
			self._plotSpec['zRange'] = [0.9,1.1]

		self.makeSystPlots()

	def makeSystPlots(self):
		'''process all yield plots into effs and systs'''

		# zero out the above-diagonal region
		for a in range(1, self._yieldPlots['central'].GetNbinsX()+1):
			xval = self._yieldPlots['central'].GetXaxis().GetBinCenter(a)
			for b in range(1, self._yieldPlots['central'].GetNbinsY()+1):
				yval = self._yieldPlots['central'].GetYaxis().GetBinCenter(b)
				if xval - yval < 0. or a < 0.:
					bin = self._yieldPlots['central'].FindBin(float(a), float(b))
					self._yieldPlots['central'].SetBinContent(bin, 0.)
					self._yieldPlots['up'].SetBinContent(bin, 0.)
					self._yieldPlots['nocuts'].SetBinContent(bin, 0.)
					if self._yieldPlots['down']: self._yieldPlots['down'].SetBinContent(bin, 0.)
		
		# create a load of effMap objects for each variation
		self._effs['central'] 	= effMap(self._yieldPlots['central'], self._yieldPlots['nocuts'])
		self._effs['up'] 		= effMap(self._yieldPlots['up'], self._yieldPlots['nocuts'])
		self._effs['up_change'] = effMap(self._effs['up']._hist, self._effs['central']._hist,
											self._effs['up']._errorHist,
											self._effs['central']._errorHist)

		if not self._cutSyst:# and self._test != "LeptonVeto":
			# only shift to centre around 0 if it's not a cut systematic
			self._effs['up_change'].shiftCentre(-1)

		# if 2-way syst vari exists, do down variation
		if self._yieldPlots['down']:
			self._effs['down'] 			= effMap(self._yieldPlots['down'],
													self._yieldPlots['nocuts'])
			self._effs['down_change'] 	= effMap(self._effs['down']._hist,
													self._effs['central']._hist,
													self._effs['down']._errorHist,
													self._effs['central']._errorHist)
			self._effs['down_change'].shiftCentre(-1)

		if self._cutSyst:
			### cut hist scenario ###
			# invert up_change to represent a cut efficiency
			# FIXME: THIS MAY NOT WORK BY CHANGING UP HERE
			self._effs['up_change'].invertHist()

		# if self._test == "LeptonVeto":
		# 	self._effs['up_change'].shiftCentre(-1)

		# calculate overall syst value
		tmp_hist = self._effs['central']._hist.Clone() # cheeky copy of some old histo
		for i in range(1, self._nBins+1):
			
			# skip null points
			if abs(self._effs['up_change']._hist.GetBinContent(i)) == 666.:
				tmp_hist.SetBinContent(i, -666.)
				continue

			pos_syst = abs(self._effs['up_change']._hist.GetBinContent(i))
			if self._effs['down_change']:
				neg_syst = abs(self._effs['down_change']._hist.GetBinContent(i))
			else:
				neg_syst = -1.

			# pick the bigger of both variations
			this_syst = pos_syst if pos_syst >= neg_syst else neg_syst
			
			tmp_hist.SetBinContent(i, this_syst)
		
		# create effMap from new total syst hist
		self._syst = effMap(tmp_hist)
		del tmp_hist

	def print_all(self, label = "", plotText = False):
		'''print syst output plots'''

		if plotText:
			if self._model not in ["T2cc", "T2_4body"]:
				print ">>> Warning: systMap: print_all: Text plot will look shit for this model,"
				print "    so it will be turned off. Remember, looks are everything in HEP."
				self._plotString = "colz"
			else:
				self._plotString = "colztext"
				r.gStyle.SetPaintTextFormat("0.4f");

		# create instance of a multi page PDF
		pdf0 = multiPagePDF(outFileName = "out/%s_%s_%s.pdf" % (self._model, self._test, label),
							title = "Systematics - %s - %s" % (self._test, label))

		# setup fine grain z-axis
		sutils.set_palette()

		# draw each variation
		for key in ["central", "up", "up_change", "down", "down_change"]:
			if "down" not in key or self._effs['down']:
				self.draw_plot(self._effs[key],
								pdf0,
								"Efficiency %s %s - %s" % (self._test, key, label),
								"Acceptance",
								1. if "change" in key else 0.)

		if not self._cutSyst:
			# draw total systematic (don't plot for cut systematics)
			self.draw_plot(self._syst,
							pdf0,
							"%s Systematic - %s" % (self._test, label), 
							"Systematic Value",
							0.)

		pdf0.close()

	def draw_plot(self, effMap = None, pdfFile = None, title = "", zTitle = "", shiftZ = 0.):
		'''draw a single plot on a single page of pdfFile'''

		hist = effMap._hist.Clone()
		hist.Draw(self._plotString)
		hist.SetTitle(title)
		if "change" in title:
			hist.GetZaxis().SetTitle(zTitle + " change")
		else:
			hist.GetZaxis().SetTitle(zTitle)

		if self._cutSyst:
			# remove the z-axis shift for cut-based systematics
			shiftZ = 0.	
		self.setDetails(hist, shiftZ)

		#hack to fix systematic zRange
		if "Syst" in title and not self._cutSyst:
			hist.GetZaxis().SetRangeUser(0., self._plotSpec['zRange'][1]-1.)

		if "text" in self._plotString:
			hist.SetMarkerSize(0.8)

		# draw all the stats numbers
		num0 = r.TLatex(0.151,0.8,"#scale[0.6]{avg = %.4f}" % effMap._mean)
		num0.SetNDC()
		num0.Draw("same")

		num1 = r.TLatex(0.15,0.77,"#scale[0.6]{RMS= %.4f}" % effMap._rms)
		num1.SetNDC()
		num1.Draw("same")

		num2 = r.TLatex(0.15,0.74,"#scale[0.6]{min = %.4f}" % effMap._min)
		num2.SetNDC()
		num2.Draw("same")

		num3 = r.TLatex(0.15,0.71,"#scale[0.6]{max = %.4f}" % effMap._max)
		num3.SetNDC()
		num3.Draw("same")

		pdfFile.AddPage()

	def setDetails(self, hist=None, zoffset=0.):
		'''set the ranges, titles, and title text specs of a plot'''

		hist.GetXaxis().SetRangeUser(self._plotSpec['xRange'][0], self._plotSpec['xRange'][1])
		hist.GetXaxis().SetTitle(self._plotSpec['xTitle'])
		hist.GetXaxis().SetTitleOffset(1.1)
		hist.GetXaxis().SetTitleSize(0.04)

		hist.GetYaxis().SetRangeUser(self._plotSpec['yRange'][0], self._plotSpec['yRange'][1])
		hist.GetYaxis().SetTitle(self._plotSpec['yTitle'])
		hist.GetYaxis().SetTitleOffset(1.1)
		hist.GetYaxis().SetTitleSize(0.04)

		if zoffset:
			hist.GetZaxis().SetRangeUser(self._plotSpec['zRange'][0]-zoffset,
											self._plotSpec['zRange'][1]-zoffset)
		# if self._cutSyst:
			# hist.GetZaxis().SetRangeUser(self., 1.2)
		if self._cutSyst:
			hist.SetMaximum(1.)

		hist.GetZaxis().SetTitleOffset(0.95)
		hist.GetZaxis().SetTitleSize(0.05)

	def __del__(self):
		'''destructor to deal with effMap objects'''

		del self._effs
		del self._yieldPlots


###---------------------------------------------------------------------------------------------###

class multiPagePDF(object):
	'''class for producing multipage PDF file'''

	def __init__(self, outFileName = "", title = "_NoTitle_"):
		self._canv = r.TCanvas() # FIXME: pick good size
		self._fName = outFileName
		self._pageNums = True
		self._doName = True
		self._analyst = "Chris Lucas"
		self._title = title
		self._pageCtr = 1 #initiate page numbers
		self.makeTitlePage()

	def makeTitlePage(self):
		'''make a title page with title, name, date'''
		from time import strftime

		title_text = r.TLatex(0.07,0.8,"#scale[1.5]{%s}" % self._title)
		title_text.SetNDC()
		title_text.Draw("same")

		if self._doName:
			name_text = r.TLatex(0.07, 0.2, "Analyst: %s" % self._analyst)
			name_text.SetNDC()
			name_text.Draw("same")

		date_text = r.TLatex(0.07, 0.15, "#scale[0.7]{%s}" %
								strftime("%a, %d %b %Y %H:%M:%S")) # FIXME: make this text smaller
		date_text.SetNDC()
		date_text.Draw("same")

		self._canv.Print(self._fName+"(")
		r.gPad.SetRightMargin(0.15)
		r.gPad.SetLeftMargin(0.10)
		r.gPad.SetTopMargin(0.08)
		r.gPad.SetBottomMargin(0.12)

	def AddPage(self):
		'''Add page with canvas drawn and page number'''
		num = r.TLatex(0.97,0.025,"%d"%(self._pageCtr))
		num.SetNDC()
		if self._pageNums: num.Draw("same")
		self._canv.Print(self._fName)
		self._pageCtr += 1
		pass

	def close(self):
		'''close the pdf file'''
		self._canv.Print(self._fName+"]")
		pass
