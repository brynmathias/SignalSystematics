#!/usr/bin/env python

import math
import ROOT as r
import array

def ExclusionCurve(xVals = None, yVals = None):
    """Define an exclusion curve based on an array of xvalues and y values. Smooth with a tspline3"""
    if len(xVals) is not len(yVals): assert "N XVals != N YVals, please correct"
    graph = r.TGraph(len(xVals),array.array('d',xVals),array.array('d',yVals))
    graph.SetLineWidth(3)
    # return graph
    spline = r.TSpline3("spline",graph,"",yVals[0],yVals[len(yVals)-1])
    spline.SetLineWidth(3)
    return spline
    pass



def returnBinsAroundLimit(histogram = None, limitX = None, limitY = None, nBinsAroundLimit = None):
  """docstring for returnBinsAroundLimit"""
  out = []
  # First check that we have a correct input line
  if len(limitX) != len(limitY): assert "Need to have the same number of x co-ords as y co-ords"
  spline = ExclusionCurve(limitX,limitY)
  xWidth = histogram.GetXaxis().GetBinWidth(1)
  yWidth = histogram.GetYaxis().GetBinWidth(1)
  for x in range(histogram.GetNbinsX()):
    if x*xWidth < 1.:continue
    y = spline.Eval(x*xWidth)
    vals = [y + yWidth*i for i in range(nBinsAroundLimit)]+[y - yWidth*i for i in range(nBinsAroundLimit)]
    vals = set(vals)
    # vals = [y]
    for val in vals:
      # print (x*xWidth)-1, val
      bin = histogram.FindBin((x*xWidth)-1, val)
      out.append(bin)
  return out
  pass