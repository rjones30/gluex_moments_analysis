#!/usr/bin/env python3
#
# buildMomentsMatrix.py - reads from a ReactionFilter containing the selected
#                         events from analysis of a Monte Carlo signal sample
#                         with all angles generated uniform in CM angles, saved
#                         as a flat tree with the moments functions computed on
#                         the generated (with trailing underscores) as well as
#                         the reconstructed angles (without trailing underscores)
#                         of the final state by a special DSelector for the given
#                         reaction, eg. see DSelector_etapi0_moments.[Ch]. A
#                         subset of the full count of MC events is read in from
#                         from the input ROOT tree, and the square M matrix for
#                         this subset of the event sum is saved in hdf5 files.
#
# author: richard.t.jones at uconn.edu
# version: june 7, 2023

use_c_extension_library = True
threading_split_level = 1

import uproot
import numpy as np
import awkward as ak
import C_buildMomentsMatrix
import threading
import h5py
import sys

def usage():
  print("usage: buildMomentsMatrix.py <reaction_tree> <event_count> <event_offset>")
  print("  <reaction_tree> is the input tree written in uproot format (eg. treefile.root:treename)")
  print("  <event_count> is the number of accepted MC events to be summed in this subset")
  print("  <event_offset> is the multiplier of <event_count> skipped before the start of this subset")
  sys.exit(1)

def usage2():
  print("usage: >>> import buildMomentsMatrix as bmm")
  print("       >>> bmm.open('etapi0_moments.root:etapi0_moments')")
  print("       >>> events=range(100)")
  print("       >>> M,Mvar = bmm.buildMomentsMatrix_sequential()")
  print("       >>> M2,Mvar2 = bmm.buildMomentsMatrix_threaded()")
  print("       >>> outfile = 'save_moments_matrix.h5')")
  print("       >>> save_output(M, Mvar, len(events), events[0], outfile)")

if len(sys.argv) == 4:
   try:
      intree_name = sys.argv[1]
      event_count = int(sys.argv[2])
      sequence_no = int(sys.argv[3])
      skip_count = sequence_no * event_count
   except:
      usage()
else:
   mPi0,mPi0 = 0,0
   mEta,mEta_ = 0,0
   mGJ,mGJ_ = 0,0
   usage2()

def open(intree_name):
   intree = uproot.open(intree_name)
   global YmomGJ
   YmomGJ = intree["YmomGJ"].array()[:].to_numpy()
   global YmomGJ_
   YmomGJ_ = intree["YmomGJ_"].array()[:].to_numpy()
   global YmomEta
   YmomEta = intree["YmomEta"].array()[:].to_numpy()
   global YmomEta_
   YmomEta_ = intree["YmomEta_"].array()[:].to_numpy()
   global YmomPi0
   YmomPi0 = intree["YmomPi0"].array()[:].to_numpy()
   global YmomPi0_
   YmomPi0_ = intree["YmomPi0_"].array()[:].to_numpy()
   global mGJ
   mGJ = YmomGJ[0].shape[0]
   global mGJ_
   mGJ_ = YmomGJ_[0].shape[0]
   global mEta
   mEta = YmomEta[0].shape[0]
   global mEta_
   mEta_ = YmomEta_[0].shape[0]
   global mPi0
   mPi0 = YmomPi0[0].shape[0]
   global mPi0_
   mPi0_ = YmomPi0_[0].shape[0]
   return intree

def buildMomentsMatrix_sequential(events, mPi0=mPi0, mEta=mEta, mGJ=mGJ, use_c_extension_library=False):
   """
   Single-threaded implementation, for small matrices and checks
   """
   mPi0_, mEta_, mGJ_ = mPi0, mEta, mGJ
   M = np.zeros([mGJ * mEta * mPi0, mGJ_ * mEta_ * mPi0_], dtype=float)
   Mvar = np.zeros([mGJ * mEta * mPi0, mGJ_ * mEta_ * mPi0_], dtype=float)
   for iev in events:
      if use_c_extension_library:
         C_buildMomentsMatrix.add_event(M2, 
           [YmomPi0[iev], YmomEta[iev], YmomGJ[iev]],
           [YmomPi0_[iev], YmomEta_[iev], YmomGJ_[iev]])
         C_buildMomentsMatrix.add_event(Mvar2,
           [np.square(YmomPi0[iev]), np.square(YmomEta[iev]), np.square(YmomGJ[iev])],
           [np.square(YmomPi0_[iev]), np.square(YmomEta_[iev]), np.square(YmomGJ_[iev])])
      else:
         M_GJ = np.array([YmomGJ[iev][iGJ] * YmomGJ_[iev] for iGJ in range(mGJ)], dtype=float)
         Mvar_GJ = np.array([np.square(YmomGJ[iev])[iGJ] * np.square(YmomGJ_[iev]) for iGJ in range(mGJ)], dtype=float)
         for iPi0 in range(mPi0):
            for iPi0_ in range(mPi0_):
               M_Pi0 = YmomPi0[iev][iPi0] * YmomPi0_[iev][iPi0_]
               for iEta in range(mEta):
                  m = iPi0 * mEta + iEta
                  for iEta_ in range(mEta_):
                     M_Pi0_Eta = M_Pi0 * YmomEta[iev][iEta] * YmomEta_[iev][iEta_]
                     m_ = iPi0_ *mEta_ + iEta_
                     M[m * mGJ : (m+1) * mGJ, m_ * mGJ_ : (m_+1) * mGJ_] += M_Pi0_Eta * M_GJ
                     Mvar[m * mGJ : (m+1) * mGJ, m_ * mGJ_ : (m_+1) * mGJ_] += M_Pi0_Eta**2 * Mvar_GJ
      print(f"did event {iev}")
   return M, Mvar

def buildMomentsMatrixSlice1(events, M, Mvar, iPi0, mPi0, mEta, mGJ, mPi0_, mEta_, mGJ_):
   """
   Builds a slice of M and Mvar, for a single Pi0 moment.
   """
   mstart = iPi0 * mEta * mGJ
   mend = mstart + mEta * mGJ
   M_slice = M[mstart:mend,:]
   Mvar_slice = Mvar[mstart:mend,:]
   for iev in events:
      Ymom = [YmomPi0[iev][iPi0:iPi0+1], YmomEta[iev][:mEta], YmomGJ[iev][:mGJ]]
      Ymom_ = [YmomPi0_[iev][:mPi0_], YmomEta_[iev][:mEta_], YmomGJ_[iev][:mGJ_]]
      Ysqr = [np.square(YmomPi0[iev][iPi0:iPi0+1]), np.square(YmomEta[iev][:mEta]), np.square(YmomGJ[iev][:mGJ])]
      Ysqr_ = [np.square(YmomPi0_[iev][:mPi0_]), np.square(YmomEta_[iev][:mEta_]), np.square(YmomGJ_[iev][:mGJ_])]
      C_buildMomentsMatrix.add_event(M_slice, Ymom, Ymom_)
      C_buildMomentsMatrix.add_event(Mvar_slice, Ysqr, Ysqr_)

def buildMomentsMatrixSlice2(events, M, Mvar, iPi0, iEta, mPi0, mEta, mGJ, mPi0_, mEta_, mGJ_):
   """
   Builds a smaller slice of M and Mvar, for a single Pi0,Eta moment pair.
   """
   mstart = (iPi0 * mEta + iEta) * mGJ
   mend = mstart + mGJ
   M_slice = M[mstart:mend,:]
   Mvar_slice = Mvar[mstart:mend,:]
   for iev in events:
      Ymom = [YmomPi0[iev][iPi0:iPi0+1], YmomEta[iev][iEta:iEta+1], YmomGJ[iev][:mGJ]]
      Ymom_ = [YmomPi0_[iev][:mPi0_], YmomEta_[iev][:mEta_], YmomGJ_[iev][:mGJ_]]
      Ysqr = [np.square(YmomPi0[iev][iPi0:iPi0+1]), np.square(YmomEta[iev][iEta:iEta+1]), np.square(YmomGJ[iev][:mGJ])]
      Ysqr_ = [np.square(YmomPi0_[iev][:mPi0_]), np.square(YmomEta_[iev][:mEta_]), np.square(YmomGJ_[iev][:mGJ_])]
      C_buildMomentsMatrix.add_event(M_slice, Ymom, Ymom_)
      C_buildMomentsMatrix.add_event(Mvar_slice, Ysqr, Ysqr_)

def buildMomentsMatrix_threaded(events, mPi0=mPi0, mEta=mEta, mGJ=mGJ, threading_split_level=1):
   """
   Multi-threaded implementation, for large matrices
   """
   mPi0_, mEta_, mGJ_ = mPi0, mEta, mGJ
   M = np.zeros([mGJ * mEta * mPi0, mGJ_ * mEta_ * mPi0_], dtype=float)
   Mvar = np.zeros([mGJ * mEta * mPi0, mGJ_ * mEta_ * mPi0_], dtype=float)

   print("starting threads")
   threads = []
   if threading_split_level == 1:
      for iPi0 in range(mPi0):
         args = (events, M, Mvar, iPi0, mPi0, mEta, mGJ, mPi0_, mEta_, mGJ_)
         t = threading.Thread(target=buildMomentsMatrixSlice1, args=args)
         threads.append(t)
         t.start()
   else:
      for iPi0 in range(mPi0):
         for iEta in range(mEta):
            args = (events, M, Mvar, iPi0, iEta, mPi0, mEta, mGJ, mPi0_, mEta_, mGJ_)
            t = threading.Thread(target=buildMomentsMatrixSlice2, args=args)
            threads.append(t)
            t.start()
   print("threads joining")
   for t in threads:
      t.join()
   print("threads joined")
   return M, Mvar

def save_output(M, Mvar, event_count, skip_count, outfile):
   h5out = h5py.File(outfile, 'w')
   h5out.create_dataset("Moments", data=M)
   h5out.create_dataset("Moments_variance", data=Mvar)
   h5out.create_dataset("accepted_subset", data=event_count)
   h5out.create_dataset("skipped_accepted", data=skip_count)
   h5out.close()

if len(sys.argv) == 4:
   if skip_count >= len(YmomGJ):
      print("error in buildMomentsMatrix.py - skip count >= number of acceptence MC events available")
      print(skip_count, ">=", len(YmomGJ))
      usage()
   elif skip_count + event_count > len(YmomGJ):
      event_count = len(YmomGJ) - skip_count
      print("warning in buildMomentsMatrix.py - event_count reduced to number of acceptence MC events available")
      print(event_count, "remaining events in the accepted MC sample")

   intree = open(intree_name)
   events = range(skip_count, skip_count + event_count)
   M,Mvar = buildMomentsMatrix_threaded(events)
   outfile = "acceptance_moments_{0}_{1}.h5".format(event_count, sequence_no)
   save_output(M, Mvar, event_count, skip_count, outfile)
