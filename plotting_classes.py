# import plottingUtils as putils
import signalUtils as sutils
import numpy as np
import ROOT as r

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
	def __init__(self, nom = None, denom = None):
		self._hist = None
		self._nom = nom
		self._denom = denom
		self._mean = 0.
		self._max = 0.
		self._min = 1.
		self._rms = 0.

		if nom and denom:
			self.process()
		elif nom:
			print "> EffMap: No denom defined. Assigning nom as eff."
			self._hist = nom
			self.process()

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
		
		self._rms 	= np.var(val_array)
		self._mean 	= np.mean(val_array)
		self._min 	= np.amin(val_array)
		self._max 	= np.amax(val_array)

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
				self._hist.SetBinContent(i, -1000.)


###---------------------------------------------------------------------------------------------###

class systMap(object):
	'''Simple systematic plotting'''
	
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
		self._ranges = {
					'x':[],
					'y':[],
					'z':[],
		}
		self._logZ = False
		self._fileTag = "jmulti_bmulti"
		self._nBins = self._yieldPlots['central'].GetNbinsX()* \
						self._yieldPlots['central'].GetNbinsY() + 1000

		# get model specific plot details
		import plotDetails as pdets
		try:
			self._plotSpec = pdets.modelPlotDetails[self._model]
		except KeyError:
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
			self._plotSpec['zRange'] = [0.9,1.1]

		self.makeSystPlots()

	def makeSystPlots(self):
		'''process all yield plots into effs and systs'''

		# zero out the above-diagonal region
		for a in range(1, self._yieldPlots['central'].GetNbinsX()+1):
			xval = self._yieldPlots['central'].GetXaxis().GetBinCenter(a)
			for b in range(1, self._yieldPlots['central'].GetNbinsY()+1):
				ybinval = self._yieldPlots['central'].GetYaxis().GetBinCenter(b)
				if xval - ybinval < 0. or a < 0.:
					bin = self._yieldPlots['central'].FindBin(float(a), float(b))
					self._yieldPlots['central'].SetBinContent(bin, 0.)
					self._yieldPlots['up'].SetBinContent(bin, 0.)
					self._yieldPlots['nocuts'].SetBinContent(bin, 0.)
					if self._yieldPlots['down']: self._yieldPlots['down'].SetBinContent(bin, 0.)
		
		# create a load of effMap objects for each variation
		self._effs['central'] 	= effMap(self._yieldPlots['central'], self._yieldPlots['nocuts'])
		self._effs['up'] 		= effMap(self._yieldPlots['up'], self._yieldPlots['nocuts'])
		self._effs['up_change'] = effMap(self._effs['up']._hist, self._effs['central']._hist)
		self._effs['up_change'].shiftCentre(-1)
		if self._yieldPlots['down']:
			self._effs['down'] 			= effMap(self._yieldPlots['down'],
													self._yieldPlots['nocuts'])
			self._effs['down_change'] 	= effMap(self._effs['down']._hist,
													self._effs['central']._hist)
			self._effs['down_change'].shiftCentre(-1)
	
		# calculate overall syst value
		tmp_hist = self._effs['central']._hist.Clone() # cheeky copy of some old histo
		for i in range(1, self._nBins+1):
			
			# skip null points
			if abs(self._effs['up_change']._hist.GetBinContent(i)) == 1000.:
				tmp_hist.SetBinContent(i, -1000.)
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


	def print_all(self, label = ""):
		'''print syst output plots'''

		# create instance of a multi page PDF
		pdf0 = multiPagePDF(outFileName = "banana.pdf", title = "Systematics - %s" % label)

		# setup fine grain z-axis
		sutils.set_palette()

		# draw each variation
		for key in ["central", "up", "up_change", "down", "down_change"]:
			self.draw_plot(self._effs[key]._hist, pdf0, "Efficiency %s %s - %s" % (self._test,
							key, label), "Acceptance change" if "change" in key else "Acceptance",
							1. if "change" in key else 0.)

		# draw total systematic
		self.draw_plot(self._syst._hist, pdf0, "%s Systematic - %s" % (self._test, label), 
						"Systematic Value", True)

		pdf0.close()

	def draw_plot(self, hist = None, pdfFile = None, title = "", zTitle = "", shiftZ = 0.):
		'''draw a single plot on a single page of pdfFile'''

		hist.Draw("colz")
		hist.SetTitle(title)
		hist.GetZaxis().SetTitle(zTitle)
		self.setDetails(hist, shiftZ)

		pdfFile.AddPage()



	def setDetails(self, hist=None, zoffset=0.):

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
		hist.GetZaxis().SetTitleOffset(1.0)
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
		self._title = title
		self._pageCtr = 1 #initiate page numbers
		self.makeTitlePage()

	def makeTitlePage(self):
		'''make a title page with title, name, date'''
		from time import strftime

		title_text = r.TLatex(0.09,0.8,"#scale[1.5]{%s}" % self._title)
		title_text.SetNDC()
		title_text.Draw("same")

		if self._doName:
			name_text = r.TLatex(0.09, 0.2, "Analyst: Chris Lucas")
			name_text.SetNDC()
			name_text.Draw("same")

		date_text = r.TLatex(0.09, 0.15, "#scale[0.7]{%s}" %
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
