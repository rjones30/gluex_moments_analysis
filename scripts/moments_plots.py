#!/usr/bin/env python3
#
# moments_plots.py - useful functions for making plots of the
#                    output from DSelector_etapi0_moments.
#
# author: richard.t.jones at uconn.edu
# version: june 19, 2023
#
# Assumes you have already filled the 144 histograms in moments_plots.root
# and all you want to do is to make pretty multi-pane plots of the profiles.

import ROOT

def usage():
   print("usage: python3")
   print(" >>> import moments_plots")
   print(" >>> moments_plots.plot_raw_moments()")

def plot_raw_moments(inputfile="moments_plots.root"):
   rootfile = ROOT.TFile(rootfile)
   hprof = []
   for im in range(0,49):
      h = [inputfile.Get(f"h{it}m{im}") for it in range(3)]
      for it in range(3):
         h[it].GetXaxis().SetTitle("mass(Eta,Pi0) [GeV]")
         h[it].GetYaxis().SetTitle("average moment")
      h[0].SetLineColor(1)
      h[1].SetLineColor(2)
      h[2].SetLineColor(3)
      hprof.append(h)
   c1 = ROOT.TCanvas("c1", "c1", 1000, 1000);
   c1.Divide(3,4)
   for im in range(0,48):
      c1.cd((im % 12) + 1)
      hmax = max([hprof[im][it].GetMaximum() for it in range(3)])
      hmin = min([hprof[im][it].GetMinimum() for it in range(3)])
      hprof[im][0].SetMaximum(hmax + (hmax - hmin) * 0.1)
      hprof[im][0].SetMinimum(hmin - (hmax - hmin) * 0.1)
      hprof[im][0].SetTitle(f"moment {im}, |t|=0.2(black), 0.5(red), 0.9(green)")
      hprof[im][0].Draw("hist")
      hprof[im][1].Draw("hist same")
      hprof[im][2].Draw("hist same")
      if (im % 12) == 11:
         c1.Update()
         ans = input("press enter to continue")
   return hprof
