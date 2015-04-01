import ROOT as r
from sys import exit


# pdf_file = r.TFile.Open("envCvRelHist.root", 'READ')
# pdf_hist = pdf_file.Get("acc_cvRel_m0_m12")
pdf_file = r.TFile.Open("T2_4body_shifted.root", 'READ')
pdf_hist = pdf_file.Get("acc_cvRel_m0_m12_incl_ra1cat")
# new hist to match other T2cc binning
new_file = r.TFile.Open("T2_4body_systematics_shifted.root", 'RECREATE')
# new_hist = r.TH2D("new_pdf", "new_pdf", 13, 62.5, 387.5, 70, 7.5, 357.5)
new_hist = r.TH2D("new_pdf", "new_pdf", 38, 87.5, 1037.5, 39, -12.5, 962.5)

vals_dict = {}

for xbin in range(1, pdf_hist.GetXaxis().GetNbins()+1):
	
	xloed = pdf_hist.GetXaxis().GetBinLowEdge(xbin)

	vals_dict[xloed] = {}

	for ybin in range(1, pdf_hist.GetYaxis().GetNbins()+1):
		val = pdf_hist.GetBinContent(xbin, ybin)
		if val:
			print pdf_hist.GetYaxis().GetBinLowEdge(ybin)
			vals_dict[xloed][pdf_hist.GetYaxis().GetBinLowEdge(ybin)] = val

for xkey in vals_dict:
	# print xkey, len(vals_dict[xkey].keys())
	for ykey in vals_dict[xkey]:
		# print "  %.1f: %f" % (ykey, vals_dict[xkey][ykey])
		# print "pdfhist:", xkey, ykey
		target_bin = new_hist.FindBin(xkey, ykey)
		new_hist.SetBinContent(target_bin, vals_dict[xkey][ykey])

for xbin in range(1, new_hist.GetXaxis().GetNbins()+1):
	last_val = 0.
	for ybin in range(1, new_hist.GetYaxis().GetNbins()+1):
		val = new_hist.GetBinContent(xbin, ybin)
		if not val:
			if last_val:
				new_hist.SetBinContent(xbin, ybin, last_val)
		last_val = val

new_hist.Write()
new_file.Close()