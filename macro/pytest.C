#!/usr/bin/env python

import array

x = array.array('f', [1.3, 3.0, 6.8, 9.0])
y = array.array('f', [1.2, 4.1, 5.8, 8.5])
xerr_h = array.array('f', [0.2, 0.1, 0.5, 0.2])
xerr_l = array.array('f', [0.4, 0.2, 0.3, 0.3])
yerr_h = array.array('f', [0.3, 0.2, 0.3, 0.5])
yerr_l = array.array('f', [0.2, 0.2, 0.4, 0.6])

xmin = -2
xmax = 12
ymin = -2
ymax = 12

# ----------------------------------
#  ROOT 
# ----------------------------------
# http://root.cern.ch/root/html/TGraphAsymmErrors.html

import ROOT

ROOT.gStyle.SetTitleBorderSize(0)
can = ROOT.TCanvas("can", "can", 10, 10, 500, 500)
can.cd()
graph = ROOT.TGraphAsymmErrors(len(x), x, y, xerr_l, xerr_h, yerr_l, yerr_h)
graph.SetTitle("Title;X;Y")
graph.SetMarkerStyle(20)
graph.SetMarkerSize(1.1)
graph.SetMarkerColor(2)
graph.SetLineColor(2)
graph.Draw("AP")

tmp = ROOT.gPad.DrawFrame(xmin, ymin, xmax, ymax, 
                          "%s;%s;%s" % (graph.GetTitle(), 
                                        graph.GetXaxis().GetTitle(), 
                                        graph.GetYaxis().GetTitle()))
graph.Draw("P same")
can.Print("root_grapherror.pdf")

# ----------------------------------
#  Matplotlib
# ----------------------------------
# http://matplotlib.sourceforge.net/examples/pylab_examples/errorbar_demo.html

import pylab
from matplotlib import rc
rc('text', usetex=True)
rc('font',**{'family':'sans-serif','sans-serif':['Helvetica'], 
             'size':'14'})

fig = pylab.figure(figsize=(7,7))
ax = fig.add_subplot(111)
ax.set_xlim(xmin, xmax)
ax.set_ylim(ymin, ymax)
ax.set_xlabel('X')
ax.set_ylabel('Y')
ax.set_title('title')
ax.errorbar(x, y, xerr=[xerr_l, xerr_h], yerr=[yerr_l, yerr_h], fmt='ro')

fig.savefig("pylab_grapherror.pdf")
