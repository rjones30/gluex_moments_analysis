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

from concurrent.futures import ThreadPoolExecutor

def usage():
  print("usage: >>> import buildMomentsMatrix as bmm")
  print("  >>> bmm.open('etapi0_moments_0.root:etapi0_moments')")
  print("  >>> acc_events = bmm.select_events(massEtaPi0_limits=(1.2,1.5), abst_limits=(0,0.2))")
  print("  >>> M,Mvar = bmm.buildMomentsMatrix_sequential(acc_events)")
  print("  >>> bmm.open('etapi0_moments_1.root:etapi0_moments')")
  print("  >>> acc_events2 = bmm.select_events(massEtaPi0_limits=(1.2,1.5), abst_limits=(0,0.2))")
  print("  >>> M2,Mvar2 = bmm.buildMomentsMatrix_threaded(acc_events2)")
  print("  >>> M += M2")
  print("  >>> Mvar += Mvar2")
  print("  >>> acc_events += acc_events2")
  print("  >>> # more of the above, till all accepted moments trees are read ...")
  print("  >>> bmm.open('generated_moments.root:etapi0_moments')")
  print("  >>> gen_events=bmm.select_events(massEtaPi0_limits=(1.2,1.5), abst_limits=(0,0.2))")
  print("  >>> outfile = 'save_moments_matrix.h5')")
  print("  >>> bmm.save_output(M, Mvar, acc_events, gen_events, outfile)")

upcache = {}
mPi0,mPi0 = 0,0
mEta,mEta_ = 0,0
mGJ,mGJ_ = 0,0
usage()

def open(intree_name):
   global intree
   try:
      intree.close()
   except:
      pass
   upcache.clear()
   intree = uproot.open(intree_name)
   global mGJ
   global mEta
   global mPi0
   mGJ = intree["momentsGJ"].array()[0]
   mEta = intree["momentsEta"].array()[0]
   mPi0 = intree['momentsPi0'].array()[0]
   return intree

def select_events(massEtaPi0_limits=(0,99), abst_limits=(0,99), model=0, 
                  start=0, stop=None, maxweight=0, use_generated=0):
   """
   Select a subset of events from the most recently loaded tree, respecting the
   the specified limits in massEtaPi0 and abst, starting at event start and continuing
   until event stop. If maxweight > 0 then perform accept-reject on the model weight,
   otherwise all events are accepted. Return a list of selected event numbers and
   their weights. If maxweight=0 then all weights are zero.
   """
   events = []
   weights = []
   columns = {'mX': "massEtaPi0", '|t|': "abst"}
   if use_generated and "massEtaPi0_" in intree:
      columns['mX'] = "massEtaPi0_"
      columns['|t|'] = "abst_"
   if model == 1:
      columns['mom'] = "model1moment"
      if "YmomGJ_" in intree:
         columns['Y'] = "YmomGJ_"
      else:
         columns['Y'] = "YmomGJ"

   decomp = ThreadPoolExecutor(32)
   arrays = intree.arrays(columns.values(), entry_start=start, entry_stop=stop,
                          decompression_executor=decomp, array_cache=upcache,
                          library="np")

   massEtaPi0 = arrays[columns['mX']]
   abst = arrays[columns['|t|']]
   if model == 1:
      YmomGJ = arrays[columns['Y']]
      model1moment = arrays[columns['mom']]
   for i in range(len(massEtaPi0)):
      if massEtaPi0[i] > massEtaPi0_limits[0] and massEtaPi0[i] < massEtaPi0_limits[1]:
         if abst[i] > abst_limits[0] and abst[i] < abst_limits[1]:
            weight = 0
            if model == 1:
               weight = np.dot(YmomGJ[i], model1moment[i])
               if weight <= 0:
                  print("error in select_events: non-positive event weight found",
                        weight, model1moment[i], YmomGJ[i])
                  return [],[]
               elif weight > maxweight and maxweight > 0:
                  print("error in select_events: event with weight > maxweight found",
                        weight, model1moment[i], YmomGJ[i])
                  return [],[]
               elif weight < maxweight * np.random.uniform():
                  continue
            events.append(i)
            weights.append(weight)
   return events,weights

def compute_moments(events, mPi0=0, mEta=0, mGJ=0, weights=[], use_generated=0):
   if mPi0 == 0:
      mPi0 = globals()['mPi0']
   if mEta == 0:
      mEta = globals()['mEta']
   if mGJ == 0:
      mGJ = globals()['mGJ']
   columns = {}
   if use_generated and "YmomGJ_" in intree:
      columns['Y'] = "YmomGJ_"
   else:
      columns['Y'] = "YmomGJ"

   decomp = ThreadPoolExecutor(32)
   arrays = intree.arrays(columns.values(), 
                          decompression_executor=decomp, array_cache=upcache,
                          library="np")

   YmomGJ = arrays[columns['Y']]
   moments = np.zeros([mPi0 * mEta * mGJ], dtype=float)
   momentsvar = np.zeros([mPi0 * mEta * mGJ], dtype=float)
   for i in range(len(events)):
      iev = events[i]
      if weights:
         moments += YmomGJ[iev] * weights[i]
         momentsvar += np.square(YmomGJ[iev] * weights[i])
      else:
         moments += YmomGJ[iev]
         momentsvar += np.square(YmomGJ[iev])
   return moments,momentsvar

def buildMomentsMatrix_sequential(events, mPi0=0, mEta=0, mGJ=0, use_c_extension_library=False):
   """
   Single-threaded implementation, for small matrices and checks
   """
   if mPi0 == 0:
      mPi0 = globals()['mPi0']
   if mEta == 0:
      mEta = globals()['mEta']
   if mGJ == 0:
      mGJ = globals()['mGJ']
   if mPi0 * mEta * mGJ > 1e5:
      print("This M matrix is too large to fit in memory, total size",
            (mPi0 * mEta * mGJ)**2 * 8 / 1024 / 1024 / 1024., "GB")
      print("Reduce the number of moments or split M into subblocks.")
      return 0,0
   columns = {'YmomPi0': "YmomPi0", 'YmomEta': "YmomEta", 'YmomGJ': "YmomGJ"}
   if "YmomPi0_" in intree:
      columns['YmomPi0_'] = "YmomPi0_"
      columns['YmomEta_'] = "YmomEta_"
      columns['YmomGJ_'] = "YmomGJ_"

   decomp = ThreadPoolExecutor(32)
   arrays = intree.arrays(columns.values(),
                          decompression_executor=decomp, array_cache=upcache,
                          library="np")

   YmomPi0 = arrays[columns['YmomPi0']]
   YmomEta = arrays[columns['YmomEta']]
   YmomGJ = arrays[columns['YmomGJ']]
   if 'YmomPi0_' in columns:
      YmomPi0_ = arrays[columns['YmomPi0_']]
      YmomEta_ = arrays[columns['YmomEta_']]
      YmomGJ_ = arrays[columns['YmomGJ_']]
   mPi0_, mEta_, mGJ_ = mPi0, mEta, mGJ
   M = np.zeros([mGJ * mEta * mPi0, mGJ_ * mEta_ * mPi0_], dtype=float)
   Mvar = np.zeros([mGJ * mEta * mPi0, mGJ_ * mEta_ * mPi0_], dtype=float)
   for iev in events:
      if use_c_extension_library:
         C_buildMomentsMatrix.add_event(M, 
           [YmomPi0[iev], YmomEta[iev], YmomGJ[iev]],
           [YmomPi0_[iev], YmomEta_[iev], YmomGJ_[iev]])
         C_buildMomentsMatrix.add_event(Mvar,
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
                     m_ = iPi0_ * mEta_ + iEta_
                     M[m * mGJ : (m+1) * mGJ, m_ * mGJ_ : (m_+1) * mGJ_] += M_Pi0_Eta * M_GJ
                     Mvar[m * mGJ : (m+1) * mGJ, m_ * mGJ_ : (m_+1) * mGJ_] += M_Pi0_Eta**2 * Mvar_GJ
   return M, Mvar

def buildMomentsMatrixSlice1(events, M, Mvar, iPi0, mPi0, mEta, mGJ, mPi0_, mEta_, mGJ_,
                             YmomPi0, YmomPi0_, YmomEta, YmomEta_, YmomGJ, YmomGJ_):
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

def buildMomentsMatrixSlice2(events, M, Mvar, iPi0, iEta, mPi0, mEta, mGJ, mPi0_, mEta_, mGJ_,
                             YmomPi0, YmomPi0_, YmomEta, YmomEta_, YmomGJ, YmomGJ_):
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

def buildMomentsMatrix_threaded(events, mPi0=0, mEta=0, mGJ=0, threading_split_level=1):
   """
   Multi-threaded implementation, for large matrices
   """
   if mPi0 == 0:
      mPi0 = globals()['mPi0']
   if mEta == 0:
      mEta = globals()['mEta']
   if mGJ == 0:
      mGJ = globals()['mGJ']
   if mPi0 * mEta * mGJ > 1e5:
      print("This M matrix is too large to fit in memory, total size",
            (mPi0 * mEta * mGJ)**2 * 8 / 1024 / 1024 / 1024., "GB")
      print("Reduce the number of moments or split M into subblocks.")
      return 0,0
   columns = {'YmomPi0': "YmomPi0", 'YmomEta': "YmomEta", 'YmomGJ': "YmomGJ"}
   if "YmomPi0_" in intree:
      columns['YmomPi0_'] = "YmomPi0_"
      columns['YmomEta_'] = "YmomEta_"
      columns['YmomGJ_'] = "YmomGJ_"

   decomp = ThreadPoolExecutor(32)
   arrays = intree.arrays(columns.values(),
                          decompression_executor=decomp, array_cache=upcache,
                          library="np")

   YmomPi0 = arrays[columns['YmomPi0']]
   YmomEta = arrays[columns['YmomEta']]
   YmomGJ = arrays[columns['YmomGJ']]
   if 'YmomPi0_' in columns:
      YmomPi0_ = arrays[columns['YmomPi0_']]
      YmomEta_ = arrays[columns['YmomEta_']]
      YmomGJ_ = arrays[columns['YmomGJ_']]
   mPi0_, mEta_, mGJ_ = mPi0, mEta, mGJ
   M = np.zeros([mGJ * mEta * mPi0, mGJ_ * mEta_ * mPi0_], dtype=float)
   Mvar = np.zeros([mGJ * mEta * mPi0, mGJ_ * mEta_ * mPi0_], dtype=float)

   print("starting threads")
   threads = []
   if threading_split_level == 1:
      for iPi0 in range(mPi0):
         args = (events, M, Mvar, iPi0, mPi0, mEta, mGJ, mPi0_, mEta_, mGJ_,
                 YmomPi0, YmomPi0_, YmomEta, YmomEta_, YmomGJ, YmomGJ_)
         t = threading.Thread(target=buildMomentsMatrixSlice1, args=args)
         threads.append(t)
         t.start()
   else:
      for iPi0 in range(mPi0):
         for iEta in range(mEta):
            args = (events, M, Mvar, iPi0, iEta, mPi0, mEta, mGJ, mPi0_, mEta_, mGJ_,
                    YmomPi0, YmomPi0_, YmomEta, YmomEta_, YmomGJ, YmomGJ_)
            t = threading.Thread(target=buildMomentsMatrixSlice2, args=args)
            threads.append(t)
            t.start()
   print("threads joining")
   for t in threads:
      t.join()
   print("threads joined")
   return M, Mvar

def save_output(M, Mvar, acc_events, gen_events, outfile):
   h5out = h5py.File(outfile, 'w')
   h5out.create_dataset("Moments", data=M)
   h5out.create_dataset("Moments_variance", data=Mvar)
   h5out.create_dataset("accepted_subset", data=len(acc_events))
   h5out.create_dataset("generated_subset", data=len(gen_events))
   h5out.create_dataset("accepted_start", data=acc_events[0])
   h5out.create_dataset("generated_start", data=gen_events[0])
   h5out.close()
