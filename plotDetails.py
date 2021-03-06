def newDetails(xRange = [], yRange = [], xDMRange = [], yDMRange = [], xTitle = "", yTitle = ""):
	d = {}
	d['xRange'] 	= xRange
	d['yRange'] 	= yRange
	d['xDMRange'] 	= xDMRange
	d['yDMRange'] 	= yDMRange
	d['xTitle'] 	= xTitle
	d['yTitle'] 	= yTitle

	return d

mapRanges = {
		"T1":		[[0.,1000.],[0.,1000.]],
		"T2":		[[100.,1400.],[0.,1400.]],
		"T1ttcc":	[[250.,1000.],[250.,1000.]],
		"T2cc":		[[80.,400.],[0.,360.]],
		"T2bb":		[[80.,400.],[0.,400.]],
		"T2tt":		[[0.,1000.],[0.,1000.]],
		"T2bw_0p25":[[0.,1000.],[0.,1000.]],
		"T2bw_0p75":[[0.,1000.],[0.,1000.]],
		"T2_4body":	[[80.,425.],[0.,425.]],
}

mapDMRanges = {
		"T1":		[[0.,1000.],[0.,1000.]],
		"T2":		[[100.,1400.],[0.,1400.]],
		"T1ttcc":	[[250.,1000.],[250.,1000.]],
		"T2cc":		[[80.,400.],[-10.,360.]],
		"T2_4body":	[[80.,400.],[-10.,360.]],
		"T2bb":		[[0., 1200.],[0., 120.]]
}

alphaTDict = {'200,275':[[0.65, None]],
	            '275,325':[[0.60, None]],
	            '325,375':[[0.55, None]],
	            '375,475':[[0.55, None]],
	            '475,575':[[0.55, None]],
	            '575,675':[[0.55, None]],
	            '675,775':[[0.55, None]],
	            '775,875':[[0.55, None]],
	            '875,975':[[0.55, None]],
	            '975,1075':[[0.55, None],],
	            '975,None':[[0.55, None],],
	            '375,None':[[0.55, None],],
	            '875,None':[[0.55, None],],
	            '1075,None':[[0.55, None],],
            }

systZRanges = {
				"JES": 			[0.7,1.3],
				"ISR": 			[0.7,1.3],
				"MHT_MET": 		[0.7,1.3],
				"DeadECAL": 	[0.7,1.3],
				"LeptonVeto": 	[0.7,1.3],
			}

modelPlotDetails = {
	"T2cc": newDetails([80.,400.],[0.,360.],[80.,400.],[0.,360.],
						"m_{Stop} (GeV)","m_{LSP} (GeV)"),
	"T2_4body": newDetails([80.,425.],[0.,425.],[80.,400.],[0.,360.],
						"m_{Stop} (GeV)","m_{LSP} (GeV)"),
	"T2tt": newDetails([0.,1000.],[0.,1000.],[0.,1000.],[0.,1000.],
						"m_{Stop} (GeV)","m_{LSP} (GeV)"),
	"T2bw_0p25": newDetails([0.,1000.],[0.,1000.],[0.,1000.],[0.,1000.],
						"m_{Stop} (GeV)","m_{LSP} (GeV)"),
	"T2bw_0p75": newDetails([0.,1000.],[0.,1000.],[0.,1000.],[0.,1000.],
						"m_{Stop} (GeV)","m_{LSP} (GeV)"),
	"T2tt": newDetails([0.,1000.],[0.,1000.],[0.,1000.],[0.,1000.],
		# "T2tt": newDetails([400., 550.],[200., 350.],[400.,450.],[400.,450.],
						"m_{Stop} (GeV)","m_{LSP} (GeV)"),
}