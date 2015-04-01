import signalUtils as sutils
import ROOT as r

lut_file = open("/Users/chrislucas/SUSY/SignalScans/effStudies/SignalSystematics/stopVectWeights/stop_weight_lut.txt")

weights = {}

for line in lut_file:
    
    line_split = (line.strip()).split("\t")
    
    if len(line_split) != 3:
        continue

    stop_mass = float(line_split[0])
    boost = float(line_split[1])
    weight = float(line_split[2])

    if stop_mass not in weights.keys():
        weights[stop_mass] = {}

    weights[stop_mass][boost] = weight

# print the weights dict
# sutils.dict_printer(weights)

r.gStyle.SetOptStat(0)

weight_distro = r.TH2D("Weights Summary", "Weights Summary", 12, 87.5, 387.5, 9, -50., 850.)
c1 = r.TCanvas()

for mass in weights:
    for boost in weights[mass]:
        weight_distro.Fill(mass, boost, weights[mass][boost])

sutils.set_palette(ncontours=100)
r.gStyle.SetPaintTextFormat("0.3f")
weight_distro.GetXaxis().SetTitle("m_{stop} (GeV)")
weight_distro.GetYaxis().SetTitle("Boost p_{T} (GeV)")
weight_distro.Draw("colztext")

c1.Print("out/stop_lut_weights.pdf")
