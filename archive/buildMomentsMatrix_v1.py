#!/usr/bin/env python3
#
# buildMomentsMatrix.py - reads from a ReactionFilter containing the selected
#                         events from analysis of a Monte Carlo signal sample
#                         with all angles generated uniform in CM angles, saved
#                         as a flat tree with the moments functions computed on
#                         the generated (with trailing underscores) as well as
#                         the reconstructed angles (without trailing underscores)
#                         of the final state by a special DSelector for the given
#                         reaction, eg. see DSelector_etapi0_moments.[Ch]. The
#                         MC events are read in chunks from from an input ROOT
#                         tree, and the square M matrix is saved in hdf5 files.
#

import sys
import ROOT
import h5py
import numpy as np

treename = "etapi0_moments"
treefile = "etapi0_moments.root"
xrootdurl = "root://cn440.storrs.hpc.uconn.edu/Gluex/resilient/"
treedir = "simulation/moments-6-2023/"
finput = ROOT.TFile.Open(xrootdurl + treedir + treefile)
mcevents = finput.Get(treename)
ninput = mcevents.GetEntries()

nev = 0
YmomGJ, YmomGJ_ = 0,0
YmomEta, YmomEta_ = 0,0
YmomPi0, YmomPi0_ = 0,0
for mcevent in mcevents:
   try:
      YmomGJ[nev,:] = np.array(mcevent.YmomGJ)
   except:
      YmomGJ = np.ndarray([ninput, mcevent.momentsGJ], dtype=float)
      YmomGJ_ = np.ndarray([ninput, mcevent.momentsGJ], dtype=float)
      YmomEta = np.ndarray([ninput, mcevent.momentsEta], dtype=float)
      YmomEta_ = np.ndarray([ninput, mcevent.momentsEta], dtype=float)
      YmomPi0 = np.ndarray([ninput, mcevent.momentsPi0], dtype=float)
      YmomPi0_ = np.ndarray([ninput, mcevent.momentsPi0], dtype=float)
      YmomGJ[nev,:] = np.array(mcevent.YmomGJ)
   YmomGJ_[nev,:] = np.array(mcevent.YmomGJ_)
   YmomEta[nev,:] = np.array(mcevent.YmomEta)
   YmomEta_[nev,:] = np.array(mcevent.YmomEta_)
   YmomPi0[nev,:] = np.array(mcevent.YmomPi0)
   YmomPi0_[nev,:] = np.array(mcevent.YmomPi0_)
   nev += 1
   if nev % 10000 == 0:
      print("read", nev, "events from input sample")
