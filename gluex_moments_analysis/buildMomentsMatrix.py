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
  print("  >>> acc_events,wgt = bmm.select_events(massEtaPi0_limits=(1.2,1.5), abst_limits=(0,0.2))")
  print("  >>> M,Mvar = bmm.buildMomentsMatrix(acc_events)")
  print("  >>> bmm.open('etapi0_moments_1.root:etapi0_moments')")
  print("  >>> acc_events2,wgt2 = bmm.select_events(massEtaPi0_limits=(1.2,1.5), abst_limits=(0,0.2))")
  print("  >>> M2,Mvar2 = bmm.buildMomentsMatrix_threaded(acc_events2)")
  print("  >>> M += M2")
  print("  >>> Mvar += Mvar2")
  print("  >>> acc_events += acc_events2")
  print("  >>> # more of the above, till all accepted moments trees are read ...")
  print("  >>> bmm.open('generated_moments.root:etapi0_moments')")
  print("  >>> gen_events,wgt = bmm.select_events(massEtaPi0_limits=(1.2,1.5), abst_limits=(0,0.2))")
  print("  >>> outfile = 'save_moments_matrix.h5')")
  print("  >>> bmm.save_output(M, Mvar, acc_events, gen_events, outfile)")

upcache = {}
mPi0,mPi0 = 0,0
mEta,mEta_ = 0,0
mGJ,mGJ_ = 0,0
usage()

def open(intree_name):
   """
   Opens a ROOT file containing a moments TTree, and prepare for subsequent
   analysis actions on the contents of the intree. The intree_name argument
   should be in the format "<rootfile_pathname>:<tree_name>. Return value
   is the uproot.TTree object constructed from intree_name. Moments are
   stored in the intree in the order of increasing L,M with M=[-L,...,+L]
   up to a maximum number of mPi0 moments for the pi0 subsystem, mEta for
   the eta subsystem, and mGJ for the Gottfried-Jackson decay angles.
   Elsewhere in this package, products of these three spherical moments
   are indexed by a mash-up of the pi0, eta, and GJ moments indices,

     moments_index = [[[(iPi0 * mEta +iEta) * mGJ + iGJ for iPi0 in range(mPi0)]
                                                        for iEta in range(mEta)]
                                                        for iGJ in range(mGJ)]

   To restrict the number of moments included in an analysis, it is often
   possible to truncate the series by specifying smaller values for mPi0,
   mEta, and mGJ in subsequent function calls that access the tree than
   are provided in the input tree itself, in which case the higher order
   moments are simply dropped from the sums and ignored. This is why most
   of the functions that read from the intree have optional arguments
   mPi0, mEta, mGJ whose defaults select the full list of what appears
   in the intree.
   """
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
                  start=0, stop=None, maxweight=0, use_generated=0,
                  use_generated_angles=1):
   """
   Select a subset of events from the most recently loaded tree, respecting the
   the specified limits in massEtaPi0 and abst, starting at event start and continuing
   until event stop. If maxweight > 0 then perform accept-reject on the model weight
   using maxweight as the cutoff, otherwise all events are accepted. Return value is
   a pair of lists containing the selected event numbers and their weights.
   If maxweight=0 then all weights are returned zero.
   """
   events = []
   weights = []
   columns = {'mX': "massEtaPi0", '|t|': "abst"}
   if use_generated and "massEtaPi0_" in intree:
      columns['mX'] = "massEtaPi0_"
      columns['|t|'] = "abst_"
   if model == 1:
      columns['mom'] = "model1moment"
      if use_generated_angles and "YmomGJ_" in intree:
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

def compute_moments(events, mPi0=0, mEta=0, mGJ=0,
                    weights=[], use_generated_angles=0):
   """
   Read TTree data from the moments tree initialized by the latest call to open(),
   and compute the sum of moments from rows listed in events. If weights is not
   empty then it must have the same length as events, and the moment sums are
   weighted by weights. If use_generated is nonzero then the moments computed
   at the generated decay angles are used instead of the reconstructed. Return
   value is a 2-tuple containing the [weighted] moments vector sum followed by
   their statistical covariance matrix.
   """
   if mPi0 == 0:
      mPi0 = globals()['mPi0']
   if mEta == 0:
      mEta = globals()['mEta']
   if mGJ == 0:
      mGJ = globals()['mGJ']
   columns = {}
   if use_generated_angles and "YmomGJ_" in intree:
      columns['Y'] = "YmomGJ_"
   else:
      columns['Y'] = "YmomGJ"

   decomp = ThreadPoolExecutor(32)
   arrays = intree.arrays(columns.values(), 
                          decompression_executor=decomp, array_cache=upcache,
                          library="np")

   YmomGJ = arrays[columns['Y']]
   Nmoments = mPi0 * mEta * mGJ
   moments = np.zeros([Nmoments], dtype=float)
   momentscov = np.zeros([Nmoments, Nmoments], dtype=float)
   for i in range(len(events)):
      iev = events[i]
      Y = np.copy(YmomGJ[iev])
      if weights:
         Y *= weights[i]
      moments += Y
      momentscov += np.outer(Y, Y)
   return moments,momentscov

def histogram_moments(events, histograms, mPi0=0, mEta=0, mGJ=0,
                      weights=[], svectors=[], use_generated=0,
                      use_generated_angles=0):
   """
   Convenience function for filling histograms of angular moments from
   the moments intree accessed in the last call to open(). Only rows that
   are found in the events list are added to the histograms. If supplied
   then the weights argument must be a list of weights of the same length
   as events, to be used when filling the histograms. The argument named
   histogram must be a list of predefined 2D histograms with mass(eta,pi0)
   on the x axis and |t| on the y axis. The length of the histogram list
   must be at least the number of moments mPi0*mEta*mGJ. If svectors is not
   empty, it must be a square matrix of dimension equal to Pi0*mEta*mGJ.
   It is applied as a matrix product to the moments vector in intree,
      new_moments = svectors @ old_moments
   to form the linear combination new_moments, and the new_moments are
   used to increment the histograms. The use_generated argument is a
   True/False flat that selects whether generated or reconstructed
   angles are used to compute the moments being histogrammed. Return
   value is the total number of histograms updated.
   """
   if mPi0 == 0:
      mPi0 = globals()['mPi0']
   if mEta == 0:
      mEta = globals()['mEta']
   if mGJ == 0:
      mGJ = globals()['mGJ']
   columns = {}
   columns = {'mX': "massEtaPi0", '|t|': "abst", 'Y': "YmomGJ"}
   if use_generated and "massEtaPi0_" in intree:
      columns['mX'] = "massEtaPi0_"
      columns['|t|'] = "abst_"
   if use_generated_angles and "YmomGJ_" in intree:
      columns['Y'] = "YmomGJ_"
   else:
      columns['Y'] = "YmomGJ"

   decomp = ThreadPoolExecutor(32)
   arrays = intree.arrays(columns.values(), 
                          decompression_executor=decomp, array_cache=upcache,
                          library="np")

   massEtaPi0 = arrays[columns['mX']]
   abst = arrays[columns['|t|']]
   YmomGJ = arrays[columns['Y']]
   for i in range(len(events)):
      iev = events[i]
      if len(svectors) > 0:
         Y = YmomGJ[iev] @ svectors
      else:
         Y = YmomGJ[iev]
      nmom = 0
      if weights:
         w = weights[i]
      else:
         w = 1
      for iPi0 in range(mPi0):
         for iEta in range(mEta):
            for iGJ in range(mGJ):
               histograms[nmom].Fill(massEtaPi0[iev], abst[iev], w * Y[iGJ])
               nmom += 1
   return nmom

def histogram_acceptance(events, histograms, mPi0=0, mEta=0, mGJ=0,
                         weights=[], svectors=[], use_generated=0,
                         use_generated_angles=1):
   """
   Computes the moments of the acceptance function using the moments intree
   accessed in the last call to open(). Only rows that are found in the 
   events list are included in the calculation. If weights are supplied
   then the weights argument must be a list of weights of the same length
   as events, to be used when computing the acceptance. The argument named
   histogram must be a list of predefined 2D histograms with mass(eta,pi0)
   on the x axis and |t| on the y axis. The length of the histogram list
   must be at least the number of moments mPi0*mEta*mGJ. If svectors is not
   empty, it must be a square matrix of dimension equal to Pi0*mEta*mGJ.
   It is applied as a matrix product to the moments vector in intree,
      new_moments = svectors @ old_moments
   to form the linear combination new_moments, and the new_moments are
   used to increment the histograms. The use_generated argument is a
   True/False flat that selects whether generated or reconstructed
   angles are used to compute the moments being histogrammed. Return value
   is the total number of acceptance histograms filled.
   """
   if mPi0 == 0:
      mPi0 = globals()['mPi0']
   if mEta == 0:
      mEta = globals()['mEta']
   if mGJ == 0:
      mGJ = globals()['mGJ']
   columns = {}
   columns = {'mX': "massEtaPi0", '|t|': "abst", 'Y': "YmomGJ"}
   if use_generated and "massEtaPi0_" in intree:
      columns['mX'] = "massEtaPi0_"
      columns['|t|'] = "abst_"
   if use_generated_angles and "YmomGJ_" in intree:
      columns['Y_'] = "YmomGJ_"

   decomp = ThreadPoolExecutor(32)
   arrays = intree.arrays(columns.values(), 
                          decompression_executor=decomp, array_cache=upcache,
                          library="np")

   massEtaPi0 = arrays[columns['mX']]
   abst = arrays[columns['|t|']]
   YmomGJ = arrays[columns['Y']]
   if "Y_" in columns:
      YmomGJ_ = arrays[columns['Y_']]
   else:
      YmomGJ_ = YmomGJ
   for i in range(len(events)):
      iev = events[i]
      if len(svectors) > 0:
         Y = YmomGJ[iev] @ svectors
         Y_ = YmomGJ_[iev] @ svectors
      else:
         Y = YmomGJ[iev]
         Y_ = YmomGJ_[iev]
      nmom = 0
      if weights:
         w = weights[i]
      else:
         w = 1
      for iPi0 in range(mPi0):
         for iEta in range(mEta):
            for iGJ in range(mGJ):
               histograms[nmom].Fill(massEtaPi0[iev], abst[iev], w * Y[iGJ] * Y_[iGJ])
               nmom += 1
   return nmom

def buildMomentsMatrix(events, mPi0=0, mEta=0, mGJ=0,
                       use_c_extension_library=False, use_generated_angles=1):
   """
   Single-threaded implementation of buildMomentsMatrix, for small matrices and checks.
   Constructs the M matrix from a subset of rows in intree selected by the events list.
   On some platforms, significant speed-up is obtained by setting use_c_extention_library
   to True, assuming it has already been built for the current platform. Return value is
   a 2-tuple consisting of the computed M matrix, followed by its covariance matrix.
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
   if use_generated_angles and "YmomPi0_" in intree:
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
   else:
      YmomPi0_ = YmomPi0
      YmomEta_ = YmomEta
      YmomGJ_ = YmomGJ
   mPi0_, mEta_, mGJ_ = mPi0, mEta, mGJ
   M = np.zeros([mGJ * mEta * mPi0, mGJ_ * mEta_ * mPi0_], dtype=float)
   Mvar = np.zeros([mGJ * mEta * mPi0, mGJ_ * mEta_ * mPi0_], dtype=float)
   for iev in events:
      if use_c_extension_library:
         C_buildMomentsMatrix.add_event(M, Mvar,
           [YmomPi0[iev][:mPi0], YmomEta[iev][:mEta], YmomGJ[iev][:mGJ]],
           [YmomPi0_[iev][:mPi0_], YmomEta_[iev][:mEta_], YmomGJ_[iev][:mGJ_]])
      else:
         M_GJ = np.array([YmomGJ[iev][iGJ] * YmomGJ_[iev] for iGJ in range(mGJ)], dtype=float)
         Mvar_GJ = np.array([YmomGJ_[iev][iGJ] * YmomGJ_[iev] for iGJ in range(mGJ)], dtype=float)
         Mvar_GJ *= sum(np.square(YmomGJ[iev])) * sum(np.square(YmomPi0[iev])) * sum(np.square(YmomEta[iev]))
         for iPi0 in range(mPi0):
            for iPi0_ in range(mPi0_):
               M_Pi0 = YmomPi0[iev][iPi0] * YmomPi0_[iev][iPi0_]
               Mvar_Pi0 = YmomPi0_[iev][iPi0] * YmomPi0_[iev][iPi0_]
               for iEta in range(mEta):
                  m = iPi0 * mEta + iEta
                  for iEta_ in range(mEta_):
                     M_Pi0_Eta = M_Pi0 * YmomEta[iev][iEta] * YmomEta_[iev][iEta_]
                     Mvar_Pi0_Eta = Mvar_Pi0 * YmomEta_[iev][iEta] * YmomEta_[iev][iEta_]
                     m_ = iPi0_ * mEta_ + iEta_
                     M[m * mGJ : (m+1) * mGJ, m_ * mGJ_ : (m_+1) * mGJ_] += M_Pi0_Eta * M_GJ
                     Mvar[m * mGJ : (m+1) * mGJ, m_ * mGJ_ : (m_+1) * mGJ_] += Mvar_Pi0_Eta * Mvar_GJ
   return M, Mvar

def _buildMomentsMatrixSlice1(events, M, Mvar, iPi0, mPi0, mEta, mGJ, mPi0_, mEta_, mGJ_,
                             YmomPi0, YmomPi0_, YmomEta, YmomEta_, YmomGJ, YmomGJ_):
   """
   Builds a slice of M and Mvar, for a single Pi0 moment. This is a helper function for
   internal use by buildMomentsMatrix_threaded, do not call directly.
   """
   mstart = iPi0 * mEta * mGJ
   mend = mstart + mEta * mGJ
   M_slice = M[mstart:mend,:]
   Mvar_slice = Mvar[mstart:mend,:]
   for iev in events:
      Ymom = [YmomPi0[iev][iPi0:iPi0+1], YmomEta[iev][:mEta], YmomGJ[iev][:mGJ]]
      Ymom_ = [YmomPi0_[iev][:mPi0_], YmomEta_[iev][:mEta_], YmomGJ_[iev][:mGJ_]]
      C_buildMomentsMatrix.add_event(M_slice, Mvar_slice, Ymom, Ymom_)

def _buildMomentsMatrixSlice2(events, M, Mvar, iPi0, iEta, mPi0, mEta, mGJ, mPi0_, mEta_, mGJ_,
                             YmomPi0, YmomPi0_, YmomEta, YmomEta_, YmomGJ, YmomGJ_):
   """
   Builds a smaller slice of M and Mvar, for a single Pi0,Eta moment pair.
   This is a helper function for internal use by buildMomentsMatrix_threaded,
   do not call directly.
   """
   mstart = (iPi0 * mEta + iEta) * mGJ
   mend = mstart + mGJ
   M_slice = M[mstart:mend,:]
   Mvar_slice = Mvar[mstart:mend,:]
   for iev in events:
      Ymom = [YmomPi0[iev][iPi0:iPi0+1], YmomEta[iev][iEta:iEta+1], YmomGJ[iev][:mGJ]]
      Ymom_ = [YmomPi0_[iev][:mPi0_], YmomEta_[iev][:mEta_], YmomGJ_[iev][:mGJ_]]
      C_buildMomentsMatrix.add_event(M_slice, Mvar_slice, Ymom, Ymom_)

def buildMomentsMatrix_threaded(events, mPi0=0, mEta=0, mGJ=0, threading_split_level=1,
                                use_generated_angles=1):
   """
   Multi-threaded implementation of buildMomentsMatrix, for large matrices.
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
   if use_generated_angles and "YmomPi0_" in intree:
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
   else:
      YmomPi0_ = YmomPi0
      YmomEta_ = YmomEta
      YmomGJ_ = YmomGJ
   mPi0_, mEta_, mGJ_ = mPi0, mEta, mGJ
   M = np.zeros([mGJ * mEta * mPi0, mGJ_ * mEta_ * mPi0_], dtype=float)
   Mvar = np.zeros([mGJ * mEta * mPi0, mGJ_ * mEta_ * mPi0_], dtype=float)

   print("starting threads")
   threads = []
   if threading_split_level == 1:
      for iPi0 in range(mPi0):
         args = (events, M, Mvar, iPi0, mPi0, mEta, mGJ, mPi0_, mEta_, mGJ_,
                 YmomPi0, YmomPi0_, YmomEta, YmomEta_, YmomGJ, YmomGJ_)
         t = threading.Thread(target=_buildMomentsMatrixSlice1, args=args)
         threads.append(t)
         t.start()
   else:
      for iPi0 in range(mPi0):
         for iEta in range(mEta):
            args = (events, M, Mvar, iPi0, iEta, mPi0, mEta, mGJ, mPi0_, mEta_, mGJ_,
                    YmomPi0, YmomPi0_, YmomEta, YmomEta_, YmomGJ, YmomGJ_)
            t = threading.Thread(target=_buildMomentsMatrixSlice2, args=args)
            threads.append(t)
            t.start()
   print("threads joining")
   for t in threads:
      t.join()
   print("threads joined")
   return M, Mvar

def save_output(M, Mvar, acc_events, gen_events, outfile):
   """
   Save results from buildMomentsMatrix() together with some properties
   of the dataset used to compute M which might be useful later on.
   """
   h5out = h5py.File(outfile, 'w')
   h5out.create_dataset("Moments", data=M)
   h5out.create_dataset("Moments_variance", data=Mvar)
   h5out.create_dataset("accepted_subset", data=len(acc_events))
   h5out.create_dataset("generated_subset", data=len(gen_events))
   h5out.create_dataset("accepted_start", data=acc_events[0])
   h5out.create_dataset("generated_start", data=gen_events[0])
   h5out.close()
