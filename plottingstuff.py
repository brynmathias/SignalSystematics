#!/usr/bin/env python
# encoding: utf-8
"""
plottingstuff.py

Created by Bryn Mathias on 2011-08-22.
Copyright (c) 2011 Imperial College. All rights reserved.
"""

import sys
import os
import ROOT as r


leg = r.TLegend(0.5, 0.5, 0.8, 0.8)
leg.SetShadowColor(0)
leg.SetBorderSize(0)
leg.SetFillStyle(4100)
leg.SetFillColor(0)
leg.SetLineColor(0)
leg.SetShadowColor(0)
leg.SetBorderSize(0)
leg.SetFillStyle(4100)
leg.SetFillColor(0)
leg.SetLineColor(0)
closeList= []
def GetHist(DataSetName = None,folder = None ,hist = "myHist",col = 0,norm = None ,Legend = "hist", rebin = None):
    # print type(DataSetName)
    if ".root" in DataSetName:
      a = r.TFile.Open(DataSetName) #open the file
    else: a = DataSetName
    if folder != None:
      b = a.Get(folder) #open the directory in the root file
      Hist = b.Get(hist)
    else:
      Hist = a.Get(hist) # get your histogram by name
    if Hist == None : Hist = r.TH1D()
    if Legend != 0:
      leg.AddEntry(Hist,Legend,"LP") # add a legend entry
    Hist.SetLineWidth(3)
    Hist.SetLineColor(col) #set colour
    if norm != None:
       Hist.Scale(norm) #if not data normilse to the data by lumi, MC is by default weighted to 100pb-1, if you have changed this change here!
    if rebin != None:
       Hist.Rebin2D(rebin,rebin)
    return Hist

def GetHistFromOpenFile(DataSetName = None,folder = None ,hist = "myHist",col = 0,norm = None ,Legend = "hist", rebin = None):
    a = DataSetName
    if folder != None:
      b = a.Get(folder) #open the directory in the root file
      Hist = b.Get(hist)
    else:
      Hist = a.Get(hist) # get your histogram by name
    if Hist == None : Hist = r.TH1D()
    if Legend != 0:
      leg.AddEntry(Hist,Legend,"LP") # add a legend entry
    Hist.SetLineWidth(3)
    Hist.SetLineColor(col) #set colour
    if norm != None:
       Hist.Scale(norm) #if not data normilse to the data by lumi, MC is by default weighted to 100pb-1, if you have changed this change here!
    if rebin != None:
       Hist.Rebin2D(rebin,rebin)
    return Hist




def Adder(hist):
  out = None
  tot = 0
  for h in hist:
    tot += h.GetBinContent(62,36)
    # print "Totting up"  ,tot
    if out == None: out = h.Clone()
    else: out.Add(h)
  """docstring for Adder"""
  return out
  pass

def weightedAdder(hist,weight):
  if len(hist) != len(weight): print "Weight list and hist list are not the same lenght"
  out = None
  for h,w in zip(hist,weight):
    if out == None:
      out = h.Clone()
      out.Multiply(w)
    else :
      h.Multiply(w)
      out.Add(h)
  return out
  pass

def nloTotalXsecMaker(weighted,notweighted):
  """Take two sets of histograms and make the NLO total cross section from it"""
  out = None
  for nom,denom in zip(weighted,notweighted):
    a = nom.Clone()
    a.Divide(denom)
    if out is None: out = a.Clone()
    else: out.Add(a)
  # print "no events per bin:", denom.GetBinContent(62,36)
  return out


def NloEffHisto(aftercuts,beforecuts,processCrossSections,TotalXsec):
  """Make CMSSM NLO plots from input histograms"""
  out = None
  for after,before,processCrossSection in zip(aftercuts,beforecuts,processCrossSections):#gives us €_p
      h = after.Clone()
      h.Divide(before)
      # print type(h),type(processCrossSection)
      h.Multiply(processCrossSection)
      if out is None: out = h.Clone()
      else:           out.Add(h)
  out.Divide(TotalXsec)
  return out
  pass





# def rebinScan(Hist):
#   """docstring for rebinScan"""
#   # find lowest bin
#   for bin in range(Hist.GetNBinsX())
#     if Hist.GetBinContent(bin) != 0:
#       low = Hist.GetBinLowEdge(bin - 100)
#       break
#   def HistogramMaxX(H):
#   Nbins = H.GetNbinsX()
#   BackItr = range(0,Nbins)
#   BackItr.reverse()
#   for x in BackItr :
#     if H.GetBinContent(x) != 0:
#       return H.GetBinLowEdge(x+1)
#   pass