#!/usr/bin/env python
import sys
import os
import ROOT as r
from plottingstuff import *
from plottingUtils import Print, MakeCumu
import math

###-------------------------------------------------------------------###
###-------------------------------------------------------------------###

r.gStyle.SetOptStat(0)
r.TH1.SetDefaultSumw2(1)
r.TH2.SetDefaultSumw2(1)

settings = {
    "mode":["JES", "ISR"][0],
    "inclHT":[False, True][0],
    "HTBins":["275_325", "325_375", "375_475", "475_575", "575_675", "675_775", "775_875", "875"],
    "deltaM":[False, True][0],
    "jMulti":["le3j", "ge4j", "eq2j", "eq3j"][0],
    "bMulti":["eq0b", "eq1b"][0:2]
}

###-------------------------------------------------------------------###

def threeToTwo(h3) :
    name = h3.GetName()
    binsz = h3.GetNbinsZ()
    # print binsz
    h2 = r.TH2D(name+"_2D",h3.GetTitle(),
                h3.GetNbinsX(), h3.GetXaxis().GetXmin(), h3.GetXaxis().GetXmax(),
                h3.GetNbinsY(), h3.GetYaxis().GetXmin(), h3.GetYaxis().GetXmax(),
                )
                
    for iX in range(1, 1+h3.GetNbinsX()) :
        for iY in range(1, 1+h3.GetNbinsY()) :
            content = h3.GetBinContent(iX, iY, 1) + h3.GetBinContent(iX, iY, 2)+ h3.GetBinContent(iX, iY, 0)
            h2.SetBinContent(iX, iY, content)
    h2.GetZaxis().SetTitle(h3.GetZaxis().GetTitle())

    #h2.RebinX(2)
    h2.RebinY(2)

    return h2

###-------------------------------------------------------------------###

def GetHist(File = None, folder = None, hist = None, Norm = None, rebinX = None, rebinY = None):
    h = None

    for f in folder:
        directory = File.Get(f)
        a = directory.Get(hist)
        if h is None:
            h = a.Clone()
        else: h.Add(a)

    return h  

###-------------------------------------------------------------------###
    
def nloTotalXsecMaker(individualXSecs = None):
    out = None
    for h in individualXSecs:
        if out is None: out = h.Clone()
        else: out.Add(h)
    return out

###-------------------------------------------------------------------###

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
                dirs.append("smsScan_%s_%s_AlphaT55_%s"%(bM, jM, ht))

    return dirs

###-------------------------------------------------------------------###

def getOutFile(model = "", htbin="", format=""):

    bMulti = settings["bMulti"]
    jMulti = settings["jMulti"]

    ## convert individual string selections to lists for iteration
    if "str" in str(type(bMulti)):
        bMulti = [bMulti]
    if "str" in str(type(jMulti)):
        jMulti = [jMulti]

    if format == "txt":
        outName = "%s_%s_systOutput_%s_%s.txt"%(settings["mode"], model, "_".join(bMulti), "_".join(jMulti))
    elif format == "pdf":
        outName = "%s_%s_%s_%s_%s%s.pdf"%(settings["mode"], model, "_".join(bMulti), "_".join(jMulti), 
            htbin, "_dM" if settings["deltaM"] else "")

    return "./out/"+outName

###-------------------------------------------------------------------###
###-------------------------------------------------------------------###

models = ["T2cc"]#,"T2bb","T1bbbb","T1","T2"]

for model in models:

    xTitle = None
    yTitle = None
    if model == "T2cc":
      if settings["deltaM"]:
        yTitle = "deltaM (GeV)"
      else:
        yTitle = "m_{LSP} (GeV)"
      xTitle = "m_{stop} (GeV)"    
    if model == "T2bb":
      xTitle = "m_{sbottom} (GeV)"
      yTitle = "m_{LSP} (GeV)"
    if model == "T2tt":
      xTitle = "m_{stop} (GeV)"
      yTitle = "m_{LSP} (GeV)"
    if model == "T1":
      xTitle = "m_{gluino} (GeV)"
      yTitle = "m_{LSP} (GeV)"
    if model == "T2":
      xTitle = "m_{squark} (GeV)"
      yTitle = "m_{LSP} (GeV)"
    if model == "T1tttt":
      xTitle = "m_{gluino} (GeV)"
      yTitle = "m_{LSP} (GeV)"
    if model == "T1bbbb":
      xTitle = "m_{gluino} (GeV)"
      yTitle = "m_{LSP} (GeV)"
    if model == "T1ttttProto":
      xTitle = "m_{gluino} (GeV)"
      yTitle = "m_{LSP} (GeV)"

    if settings["mode"]=="ISR":
      ins = "isr_"
    else:
      ins = ""

    centalRootFile73 = r.TFile.Open("./rootFiles/%s/%s/sigScan_%s_had_2012_73.7_%sbt0.0_MChi-1.0.root"%(model+"_v5", settings["mode"], model, ins))
    centalRootFile86 = r.TFile.Open("./rootFiles/%s/%s/sigScan_%s_had_2012_86.7_%sbt0.0_MChi-1.0.root"%(model+"_v5", settings["mode"], model, ins))
    centalRootFile100 = r.TFile.Open("./rootFiles/%s/%s/sigScan_%s_had_2012_100.0_%sbt0.0_MChi-1.0.root"%(model+"_v5", settings["mode"], model, ins))
    

    jesPlusRootFile73 = r.TFile.Open("./rootFiles/%s/%s/sigScan_%s_had_2012_73.7_%s+ve_bt0.0_MChi-1.0.root"%(model+"_v5", settings["mode"], model, ins))
    jesPlusRootFile86 = r.TFile.Open("./rootFiles/%s/%s/sigScan_%s_had_2012_86.7_%s+ve_bt0.0_MChi-1.0.root"%(model+"_v5", settings["mode"], model, ins))
    jesPlusRootFile100 = r.TFile.Open("./rootFiles/%s/%s/sigScan_%s_had_2012_100.0_%s+ve_bt0.0_MChi-1.0.root"%(model+"_v5", settings["mode"], model, ins))
    

    jesNegRootFile73 = r.TFile.Open("./rootFiles/%s/%s/sigScan_%s_had_2012_73.7_%s-ve_bt0.0_MChi-1.0.root"%(model+"_v5", settings["mode"], model, ins))
    jesNegRootFile86 = r.TFile.Open("./rootFiles/%s/%s/sigScan_%s_had_2012_86.7_%s-ve_bt0.0_MChi-1.0.root"%(model+"_v5", settings["mode"], model, ins))
    jesNegRootFile100 = r.TFile.Open("./rootFiles/%s/%s/sigScan_%s_had_2012_100.0_%s-ve_bt0.0_MChi-1.0.root"%(model+"_v5", settings["mode"], model, ins))
      
    HTList = []  
    if settings["inclHT"]:
      HTList = ["incl"]
    else:
      HTList = settings["HTBins"]
      oF = open(getOutFile(model=model, format="txt"), 'w')

    for htbin in HTList:

        processCrossSections = []
        cuts = []
        cutsJESNeg = []
        cutsJESPlus = []
        # cutsJESRan = []
        nocuts = []

        suf = ""
        if "275_" in htbin:
          suf="73"
        elif "325_" in htbin:
          suf="86"
        else:
          suf="100"


        c1 = Print( getOutFile(model=model, htbin=htbin, format="pdf") )
        c1.DoPageNum = False

        r.gPad.SetRightMargin(0.175)
        r.gPad.SetLeftMargin(0.15)
        r.gPad.SetTopMargin(0.08)
        r.gPad.SetBottomMargin(0.15)

        nocuts = GetHist(File = centalRootFile100,folder = ["smsScan_before",],hist = "m0_m12_mChi_noweight", Norm = None ,rebinX= 4)
        nocuts = threeToTwo(nocuts)
        
        if settings["inclHT"]:

          cutsHist = GetHist(File = centalRootFile73,folder = getRootDirs()[0:1],hist = "m0_m12_mChi_noweight", Norm = None ,rebinX= 4).Clone()
          cutsHist.Add(GetHist(File = centalRootFile86,folder = getRootDirs()[1:2],hist = "m0_m12_mChi_noweight", Norm = None ,rebinX= 4))
          cutsHist.Add(GetHist(File = centalRootFile100,folder = getRootDirs()[2:],hist = "m0_m12_mChi_noweight", Norm = None ,rebinX= 2))

          cutsJESPlusHist = GetHist(File =   jesPlusRootFile73,folder = getRootDirs()[0:1],hist = "m0_m12_mChi_noweight", Norm = None ,rebinX= 4).Clone()
          cutsJESPlusHist.Add(GetHist(File = jesPlusRootFile86,folder = getRootDirs()[1:2],hist = "m0_m12_mChi_noweight", Norm = None ,rebinX= 4))
          cutsJESPlusHist.Add(GetHist(File = jesPlusRootFile100,folder = getRootDirs()[2:],hist = "m0_m12_mChi_noweight", Norm = None ,rebinX= 4))
          
          cutsJESNegHist = GetHist(File =   jesNegRootFile73,folder = getRootDirs()[0:1],hist = "m0_m12_mChi_noweight", Norm = None ,rebinX= 4).Clone()
          cutsJESNegHist.Add(GetHist(File = jesNegRootFile86,folder = getRootDirs()[1:2],hist = "m0_m12_mChi_noweight", Norm = None ,rebinX= 4))
          cutsJESNegHist.Add(GetHist(File = jesNegRootFile100,folder = getRootDirs()[2:],hist = "m0_m12_mChi_noweight", Norm = None ,rebinX= 4))
        else:
          d = [getRootDirs()[i] for i,x in enumerate(getRootDirs()) if htbin in x]
          cutsHist = GetHist(File = eval("centalRootFile%s"%suf),folder = d, hist = "m0_m12_mChi_noweight", Norm = None ,rebinX= 4).Clone()
          cutsJESPlusHist = GetHist(File =   eval("jesPlusRootFile%s"%suf),folder = d, hist = "m0_m12_mChi_noweight", Norm = None ,rebinX= 4).Clone()
          cutsJESNegHist = GetHist(File =   eval("jesNegRootFile%s"%suf),folder = d, hist = "m0_m12_mChi_noweight", Norm = None ,rebinX= 4).Clone()

        # convert to TH2 plots
        cutsHist = threeToTwo(cutsHist)
        cutsJESPlusHist = threeToTwo(cutsJESPlusHist)
        cutsJESNegHist = threeToTwo(cutsJESNegHist)                                       
        
        if settings["deltaM"]:
            nocuts = deltaM(nocuts)
            cutsHist = deltaM(cutsHist)
            cutsJESPlusHist = deltaM(cutsJESPlusHist)
            cutsJESNegHist = deltaM(cutsJESNegHist)
        
        #cutsHist.GetXaxis().SetRangeUser(50.,300.)
        #cutsJESPlusHist.GetXaxis().SetRangeUser(50.,300.)
        #cutsJESNegHist.GetXaxis().SetRangeUser(50.,300.)
        #cutsHist.GetYaxis().SetRangeUser(0.,300.)
        #cutsJESPlusHist.GetYaxis().SetRangeUser(0.,300.)
        #cutsJESNegHist.GetYaxis().SetRangeUser(0.,300.)


        l =  [i for i in range(301)]
        for a in l:
          xbinval = cutsHist.GetXaxis().GetBinCenter(a)
          for b in l:
            ybinval = cutsHist.GetYaxis().GetBinCenter(b)
            if  xbinval - ybinval < 0. or a < 0. :
              bin = cutsHist.FindBin(float(a),float(b))
              cutsHist.SetBinContent(bin,0.)
              cutsJESPlusHist.SetBinContent(bin,0.)
              cutsJESNegHist.SetBinContent(bin,0.)
          

        if settings["mode"]=="JES":
            mini = 0.85
            maxi =1.15
        else:
            mini = 0.75
            maxi =1.25            

        c1.canvas.SetLogz()
        offset = 1.1
        # c1.Print()    
        c1.canvas.SetLogz(False)
        TotalEff = cutsHist.Clone()
        TotalEff.GetZaxis().SetTitle("Fraction of expected signal yield")
        TotalEff.GetZaxis().SetTitleOffset(offset)
        TotalEff.GetZaxis().SetTitleSize(0.05)
        TotalEff.GetXaxis().SetTitle(xTitle)

        TotalEff.SetTitle("Total Efficiency")

        # TotalEff.GetYaxis().SetLabelSize(0.04)
        TotalEff.GetYaxis().SetTitleOffset(1.3)
        TotalEff.GetYaxis().SetTitleSize(0.05)        
        TotalEff.GetYaxis().SetTitle(yTitle)
        # TotalEff.SetTitle("Efficiency (No JES)")
        TotalEff.Divide(nocuts)
        #TotalEff.GetZaxis().SetRangeUser(0., 1.)
        maxVal = TotalEff.GetMaximum()
        TotalEff.SetMaximum(0.0025)
        

        #TotalEff.Scale(100.)
        #r.gStyle.SetPaintTextFormat("0.2f %%");
        #TotalEff.SetMarkerSize(1.4)
        #r.gStyle.SetPalette(3)
        #TotalEff.Draw("COLZ TEXT40")
        TotalEff.Draw("COLZ")
        # TotalEff.Scale(100.)
        
        tot = 0.
        ctr = 0
        for i in range( TotalEff.GetNbinsX()*TotalEff.GetNbinsY() ):
            val = TotalEff.GetBinContent(i)
            if val>0.:
                tot+=val
                ctr+=1

        num0 = r.TLatex(0.17,0.85,"Average efficiency: %.3f%%"%(float(tot/ctr)*100))
        num0.SetNDC()
        num0.Draw("same")

        c1.Print()
        # TotalEff.Scale(1./100.)


        TotalEffPlus = cutsJESPlusHist.Clone()
        TotalEffPlus.GetZaxis().SetTitle("Relative change in efficiency")
        TotalEffPlus.GetZaxis().SetTitleOffset(offset)
        TotalEffPlus.SetTitle("Total Up Efficiency")
        TotalEffPlus.Divide(nocuts)
        TotalEffPlus.GetXaxis().SetTitle(xTitle)
        TotalEffPlus.GetYaxis().SetTitle(yTitle)
        TotalEffPlus.GetZaxis().SetTitleSize(0.05)
        
        # TotalEffPlus.GetYaxis().SetLabelSize(0.04)
        TotalEffPlus.GetYaxis().SetTitleOffset(1.3)
        TotalEffPlus.GetYaxis().SetTitleSize(0.05)
        TotalEffPlus.SetMaximum(maxVal) 
        #TotalEffPlus.GetZaxis().SetRangeUser(0., 1.)
        TotalEffPlus.Draw("COLZ")
        c1.Print()    


        TotalEffNeg = cutsJESNegHist.Clone()
        TotalEffNeg.GetZaxis().SetTitle("Relative change in efficiency")
        TotalEffNeg.GetZaxis().SetTitleOffset(offset)
        TotalEffNeg.SetTitle("Total Down Efficiency")
        TotalEffNeg.Divide(nocuts)
        TotalEffNeg.GetZaxis().SetTitleSize(0.05)
        
        TotalEffNeg.GetXaxis().SetTitle(xTitle)
        TotalEffNeg.GetYaxis().SetTitle(yTitle)
        # TotalEffNeg.GetYaxis().SetLabelSize(0.04)
        TotalEffNeg.GetYaxis().SetTitleOffset(1.3)
        TotalEffNeg.GetYaxis().SetTitleSize(0.05)
        TotalEffNeg.SetMaximum(maxVal)
        #TotalEffNeg.GetZaxis().SetRangeUser(0., 1.)
        TotalEffNeg.Draw("COLZ")
        c1.Print()
        

        EffOverJESPlus = TotalEffPlus.Clone()
        EffOverJESPlus.SetTitle("Up deltaEff (Up / Total)")
        EffOverJESPlus.Divide(TotalEff)
        EffOverJESPlus.Draw("COLZ")
        EffOverJESPlus.GetZaxis().SetRangeUser(0.5, 1.5)
        #c1.Print()      

        EffOverJESNeg = TotalEffNeg.Clone()
        EffOverJESNeg.SetTitle("Down deltaEff (Down / Total)")
        EffOverJESNeg.Divide(TotalEff)
        EffOverJESNeg.Draw("COLZ")
        EffOverJESNeg.GetZaxis().SetRangeUser(0.5, 1.5)
        #c1.Print()

        r.gStyle.SetOptStat(0)

        EffOverJESNegClone = EffOverJESNeg.Clone()
        EffOverJESPlusClone = EffOverJESPlus.Clone()

        #force the zRange accross all bins
        for bin in range(EffOverJESNegClone.GetNbinsX()*EffOverJESNegClone.GetNbinsY()):
          if EffOverJESNegClone.GetBinContent(bin) > 0:
            if EffOverJESNegClone.GetBinContent(bin) < mini: EffOverJESNegClone.SetBinContent(bin,mini)
            if EffOverJESNegClone.GetBinContent(bin) > maxi: EffOverJESNegClone.SetBinContent(bin,maxi)
        for bin in range(EffOverJESPlusClone.GetNbinsX()*EffOverJESPlusClone.GetNbinsY()):
          if EffOverJESPlusClone.GetBinContent(bin) > 0:
            if EffOverJESPlusClone.GetBinContent(bin) < mini: EffOverJESPlusClone.SetBinContent(bin,mini)
            if EffOverJESPlusClone.GetBinContent(bin) > maxi: EffOverJESPlusClone.SetBinContent(bin,maxi)
            
        # centre distribution around 0., set all others to 1000
        for bin in range(EffOverJESPlusClone.GetXaxis().GetNbins()*EffOverJESPlusClone.GetYaxis().GetNbins()+10000):
          if EffOverJESNegClone.GetBinContent(bin) > 0:
            EffOverJESNegClone.SetBinContent(bin,EffOverJESNegClone.GetBinContent(bin)-1.)
          else:EffOverJESNegClone.SetBinContent(bin,-1000)
          if EffOverJESPlusClone.GetBinContent(bin) > 0:
            EffOverJESPlusClone.SetBinContent(bin,EffOverJESPlusClone.GetBinContent(bin)-1.)
          else:EffOverJESPlusClone.SetBinContent(bin,-1000)
        
        #force range
        EffOverJESNegClone.SetMinimum(mini-1.)
        EffOverJESNegClone.SetMaximum(maxi-1.)
        EffOverJESPlusClone.SetMinimum(mini-1.)
        EffOverJESPlusClone.SetMaximum(maxi-1.)

        EffOverJESPlusClone.Draw("COLZ")
        #EffOverJESPlusClone.GetZaxis().SetRangeUser(-.14, .14)
        c1.Print()

        EffOverJESNegClone.Draw("COLZ")
        #EffOverJESNegClone.GetZaxis().SetRangeUser(-.14, .14)
        c1.Print()


        oneDJesMinus = r.TH1D("oneDJesMinus","",500,0.,0.5)
        oneDJesPlus = r.TH1D("oneDJesPlus","",500,0.,0.5)
        # oneDJesRan = r.TH1D("oneDJesRan","JES variation for signal efficiency 1D Projection",1000,-5.,5.)
        nEvents = r.TH1D("totEv","totEv",2500,0,2500)
        minf=0.9
        maxf=1.1
        fit =r.TF1("Gaussian","gaus",minf,maxf)
    
        ## loop over all bins
        totalBins = TotalEff.GetXaxis().GetNbins()*TotalEff.GetYaxis().GetNbins()
        for bin in range(totalBins):
            ##get deviations from "no change"
            contentMinus  = math.fabs(EffOverJESNeg.GetBinContent(bin)-1.)
            errMinus      = EffOverJESNeg.GetBinError(bin)
            contentPlus   = math.fabs(EffOverJESPlus.GetBinContent(bin)-1.)
            errPlus       = EffOverJESPlus.GetBinError(bin)
            #print contentMinus, contentPlus
            # contentRan =   EffOverJESRan.GetBinContent(bin)
            content = nocuts.GetBinContent(bin)
            ## skip if no original events, smaller than .01
            if content == 0: continue
            if EffOverJESPlus.GetBinContent(bin) < 0.01: continue
            if EffOverJESPlus.GetBinContent(bin) < 0.  or EffOverJESPlus.GetBinContent(bin) > 100.:
               nEvents.Fill(math.fabs(cutsJESPlusHist.GetBinContent(bin)-cutsHist.GetBinContent(bin)))
            ## fill the variation in 1d
            if contentMinus > 0.:
                #oneDJesMinus.Fill(math.fabs(contentMinus))
                this = oneDJesMinus.FindBin(contentMinus)
                oneDJesMinus.SetBinContent(this, contentMinus)
                oneDJesMinus.SetBinError(this, errMinus)
            if contentPlus > 0.:
                #oneDJesPlus.Fill(math.fabs(contentPlus))
                this = oneDJesPlus.FindBin(contentPlus)
                oneDJesPlus.SetBinContent(this, contentPlus)
                oneDJesPlus.SetBinError(this, errPlus)

        l =  [i * 25 for i in range(60)]


        binlist = []
        closeBins = []
        farBins = []
        closeToLine = r.TH1D("OneD_Projection_closeToLine","",400,0,0.4)
        farToLine = r.TH1D("OneD_Projection_farToLine","",400,0,0.4)
        closeToLine.GetXaxis().SetTitle("Relative change in efficiency")
        farToLine.GetXaxis().SetTitle("Relative change in efficiency")
        line = 350.
        cutoff = 450.
        for a in l:
          for b in l:
            # m_sq (m_gl) - m_LSP >= 175 && m_sq (m_gl) >= 300
            bin = EffOverJESPlus.FindBin(float(a),float(b))
            if  a - b >= line and a > cutoff :                
              # print bin
              farBins.append(bin)
            else:
              closeBins.append(bin)


        # 
        closeBins = set(closeBins)
        farBins = set(farBins)
        for bin in farBins:
          if bin in closeBins: 
            closeBins.remove(bin)
            if "T2cc" not in model: print "Bin %d is in both near and Far!!!!!!!"%bin
        TestMe = r.TH2D(EffOverJESPlus)
        TestMe.SetMinimum(0.)
        TestMe.SetMaximum(1.)
        for b in farBins: 
          if EffOverJESPlus.GetBinContent(b) >0.:farToLine.Fill(math.fabs(EffOverJESPlus.GetBinContent(b)-1.))
          if EffOverJESNeg.GetBinContent(b) >0.: farToLine.Fill(math.fabs(EffOverJESNeg.GetBinContent(b)-1.))
        for b in closeBins: 
          if EffOverJESPlus.GetBinContent(b) >0.:closeToLine.Fill(math.fabs(EffOverJESPlus.GetBinContent(b)-1.))
          if EffOverJESNeg.GetBinContent(b) >0.: closeToLine.Fill(math.fabs(EffOverJESNeg.GetBinContent(b)-1.))

        closeBinNom = closeToLine.Integral()
        farBinNorm = farToLine.Integral()
        closeToLineClone = MakeCumu(closeToLine)
        farToLineClone = MakeCumu(farToLine)
        if closeBinNom > 0.:
            closeToLineClone.Scale(1./closeBinNom)
        if farBinNorm > 0.:
            farToLineClone.Scale(1./farBinNorm)
        bin68Close = 0
        for bin in range(closeToLine.GetNbinsX()):
          # print model, closeToLine.GetBinContent(bin) , bin
          if closeToLineClone.GetBinContent(bin)  <= 0.68:
            bin68Close = bin
        bin68far = 0
        for bin in range(farToLine.GetNbinsX()):
          # print model, farToLine.GetBinContent(bin) , bin
          if farToLineClone.GetBinContent(bin)  <= 0.68:
            bin68far = bin

        # r.gStyle.SetOptStat(0)
        
        closeToLine68 = r.TH1D(closeToLine)
        for bin in range(closeToLine.GetNbinsX()):
           if bin > bin68Close:
             closeToLine68.SetBinContent(bin,0.)
             closeToLine68.SetBinError(bin,0.)
        closeToLine.Draw("hist")
        closeToLine68.SetFillColor(r.kRed)
        closeToLine68.SetLineColor(r.kRed)
        closeToLine68.Draw("samehist")
        num = r.TLatex(0.4,0.8,"68%% of events below %.3f"%(closeToLineClone.GetBinLowEdge(bin68Close)))
        num.SetNDC()
        num.Draw("same")
        
        ##c1.Print()


        farToLine68 = r.TH1D(farToLine)
        for bin in range(farToLine.GetNbinsX()):
           if bin > bin68far:
             farToLine68.SetBinContent(bin,0.)
             farToLine68.SetBinError(bin,0.)
        
        farToLine.Draw("hist")
        farToLine68.SetFillColor(r.kRed)
        farToLine68.SetLineColor(r.kRed)
        farToLine68.Draw("samehist")
        num = r.TLatex(0.4,0.8,"68%% of events below %.3f"%(farToLineClone.GetBinLowEdge(bin68far)))
        num.SetNDC()
        num.Draw("same")
        
        ##c1.Print()

        r.gStyle.SetOptStat(0)
        

        for b in range(TestMe.GetXaxis().GetNbins()*TestMe.GetYaxis().GetNbins()):
          # if b in closeBins: TestMe.SetBinContent(b,1.)
          if b in farBins: 
            TestMe.SetBinContent(b,0.5)
          elif b in closeBins: TestMe.SetBinContent(b,1.0)
          else:
            TestMe.SetBinContent(b,0.)
          if EffOverJESNeg.GetBinContent(b) < 0.01: TestMe.SetBinContent(b,0.)
        TestMe.SetTitle("TestMe")
        TestMe.Draw("COLZ")
        r.gStyle.SetOptStat(0)
        #c1.Print()


        # with the 1d distro of variations from 1., plot a normalised cumulative distro
        scalValPlus = oneDJesPlus.Integral()
        oneDJesPlus = MakeCumu(oneDJesPlus)
        oneDJesPlus.Scale(1./scalValPlus)
        oneDJesPlus.SetTitle("Up deltaEff")
        oneDJesPlus.Draw("hist")
        oneDJesPlus.GetXaxis().SetTitle("Relative change in efficiency")
        oneDJesPlus.GetXaxis().SetTitleSize(0.05)
        oneDJesPlus.Fit("Gaussian","RQ")
        fit.Draw("lsame")

        c1.Print()

        r.gStyle.SetOptStat(0)
        scalValMinus = oneDJesMinus.Integral()
        oneDJesMinus = MakeCumu(oneDJesMinus)
        oneDJesMinus.Scale(1./scalValMinus)
        oneDJesMinus.SetTitle("Down deltaEff")
        oneDJesMinus.Draw("hist")
        oneDJesMinus.GetXaxis().SetTitle("Relative change in efficiency")
        oneDJesMinus.GetXaxis().SetTitleSize(0.05)
        oneDJesMinus.Fit("Gaussian","RQ")
        fit.Draw("lsame")

        c1.Print()
        
        
        JesTotal = oneDJesPlus.Clone()
        JesTotal.Add(oneDJesMinus)
        JesTotal.Scale(1./2.)
        JesTotal.SetName("JESTotal")
        bin68 = 0
        for bin in range(JesTotal.GetNbinsX()):
          # print model, JesTotal.GetBinContent(bin) , bin
          if JesTotal.GetBinContent(bin) <= 0.68:
            bin68 = bin
        JesTotClone = r.TH1D(JesTotal)
        for bin in range(JesTotal.GetNbinsX()):
          if bin > bin68: 
            JesTotClone.SetBinContent(bin,0.)
            JesTotClone.SetBinError(bin,0.)
        JesTotal.SetTitle("Total deltaEff (Sum of Up and Down)")
        JesTotal.Draw("h")
        JesTotClone.SetFillColor(r.kRed)
        JesTotClone.Draw("sameh")
        
        num = r.TLatex(0.4,0.3,"68%% of events below %.3f"%(JesTotClone.GetBinLowEdge(bin68)))
        num.SetNDC()
        num.Draw("same")

        c1.Print()
      
        print "\n\nSYST: %f\n\n"%JesTotClone.GetBinLowEdge(bin68)
        if not settings["inclHT"]:
          oF.write("%s\t\t%.3f +\\- %.3f\n"%(htbin, JesTotClone.GetBinLowEdge(bin68), JesTotClone.GetBinError(bin68)))


        c1.close()

    if not settings["inclHT"]:
      oF.close()
