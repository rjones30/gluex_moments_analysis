#!/usr/bin/env python3
#
# analyzeMomentsMatrix.py - general utilities for examining and histogramming
#                           a Monte Carlo acceptance moments matrix.
#
# author: richard.t.jones at uconn.edu
# version: june 10, 2023

import numpy as np
import h5py
import sys
import ROOT
import math
import threading
import scipy.linalg
import time

import gluex_moments_analysis.buildMomentsMatrix as bmm

intreedir = ".."
workdir = "../trial1"
savedir = "."

def usage():
  print("usage: >>> import analyzeMomentsMatrix.py as ana")
  print(" >>> ana.open(\"filename.h5\")")
  print(" >>> h = ana.histogram_M_diagonal()")
  print(" >>> h = ana.histogram_M_diagonal_GJ0()")
  print(" >>> h = ana.histogram_M_diagonal_GJ0_Eta0()")
  print(" >>> h = ana.histogram_M_offdiagonal()")
  print(" >>> h = ana.histogram_eigenvalues(name, title, w)")
  print(" >>> h = ana.histogram_moments(name, title, moments, errors)")
  print(" >>> h = ana.analyze_moments(massEtaPi0_limits=(0,99), abst_limits=(0,99),")
  print("                             sample_subset=range(10), acceptance_subset=range(10),")
  print("                             model=1, Mmatrix=\"Msaved.h5\")")
  print(" >>> kinbins = ana.standard_kinematic_bins()")
  print(" >>> h = ana.model1_corrected_moments()")
  print(" >>> cmom,ccov,chists = ana.apply_constraints(S, moments, covariance, hists)")
  print(" >>> h = ana.histogram_moments_correlations()")

if len(sys.argv) > 1 and sys.argv[1][0] == '-':
  usage()

h5inputs = {}
h2dsupport = {}

def open(h5input):
  """
  Recover the M matrix saved in file pathname h5input,
  probably from a prior run of buildMomentsMatrix.save().
  Return value is the h5py.File object connected to h5input.
  """
  global f5
  global M
  f5 = h5py.File(h5input, 'r')
  M = f5['Moments']
  h5inputs[0] = f5
  return f5

def histogram_M_diagonal(h5input=0):
  """
  Read a M matrix from h5input (defaults to last M read from a file)
  and run some statistics over its values down the diagonal, which
  represent the acceptance eigenvalues in the ideal acceptance limit.
  Return value is a histogram of the diagonal values ordered by the
  moment index.
  """
  if not h5input in h5inputs:
    try:
      h5inputs[h5input] = open(h5input)
    except:
      print("unable to open input moments file", h5input)
      return 0
  MCsample_size = h5inputs[h5input]['generated_subset'][()]
  mdim = M.shape[0]
  hMdiag = ROOT.TH1D("hMdiag", "M diagonal elements", mdim, 0, mdim)
  hMdiag.GetXaxis().SetTitle("moment index")
  hMdiag.GetYaxis().SetTitle("M[i,i]")
  hMdiag.SetMinimum(0)
  normfact = (4 * math.pi)**3 / MCsample_size
  smallest = [1e9,0]
  largest = [0,0]
  for m in range(mdim):
    hMdiag.SetBinContent(m+1, M[m,m] * normfact)
    if M[m,m] < smallest[0]:
      smallest = M[m,m]
    elif M[m,m] > largest[0]:
      largest = M[m,m]
  print("smallest diagonal element is ", end='')
  print(smallest[0] * normfact, end='')
  print(", largest is ", end='')
  print(largest[0] * normfact)
  return hMdiag

def histogram_M_diagonal_GJ0(h5input=0):
  """
  Variant of histogram_M_diagonal that examines only the subset of
  diagonal elements corresponding to L=0,M=0 in the Gottfried-Jackson
  angles.  Return value is a histogram of the diagonal values ordered
  by the moment index with iGJ=0.
  """
  if not h5input in h5inputs:
    try:
      h5inputs[h5input] = open(h5input)
    except:
      print("unable to open input moments file", h5input)
      return 0
  MCsample_size = h5inputs[h5input]['generated_subset'][()]
  mdim = M.shape[0]
  hMdiag = ROOT.TH1D("hMdiagGJ0", "M diagonal elements with L(GJ)=0", 768, 0, 768)
  hMdiag.GetXaxis().SetTitle("moment index")
  hMdiag.GetYaxis().SetTitle("M[i,i]")
  hMdiag.SetMinimum(0)
  normfact = (4 * math.pi)**3 / MCsample_size
  for m in range(mdim):
    if (m % 49) == 0:
      hMdiag.SetBinContent(int(m/49)+1, M[m,m] * normfact)
  return hMdiag

def histogram_M_diagonal_GJ0_Eta0(h5input=0):
  """
  Variant of histogram_M_diagonal that examines only the subset of
  diagonal elements corresponding to L=0,M=0 in the eta decay angles.
  Return value is a histogram of the diagonal values ordered
  by the moment index with iEta=0.
  """
  if not h5input in h5inputs:
    try:
      h5inputs[h5input] = open(h5input)
    except:
      print("unable to open input moments file", h5input)
      return 0
  MCsample_size = h5inputs[h5input]['generated_subset'][()]
  mdim = M.shape[0]
  hMdiag = ROOT.TH1D("hMdiagGJ0Eta0", "M diagonal elements with L(GJ)=0, L(Eta)=0", 28, 0, 28)
  hMdiag.GetXaxis().SetTitle("moment index")
  hMdiag.GetYaxis().SetTitle("M[i,i]")
  hMdiag.SetMinimum(0)
  normfact = (4 * math.pi)**3 / MCsample_size
  for m in range(mdim):
    if (m % 1372) == 0:
      hMdiag.SetBinContent(int(m/1372)+1, M[m,m] * normfact)
  return hMdiag

def histogram_M_offdiagonal(h5input=0):
  """
  Read a M matrix from h5input (defaults to last M read from a file)
  and run some statistics over its values off the diagonal, which
  represent the departure from the ideal acceptance limit. Return
  value is a 4-tuple of 1D histograms displaying the distribution
  of the normalized off-diagonal matrix elements of M, ordered as
    [0] : elements of M above the diagonal
    [1] : elements of M below the diagonal
    [2] : statistical error on M elements above the diagonal
    [3] : statistical error on M elements below the diagonal
  """
  if not h5input in h5inputs:
    try:
      h5inputs[open(h5input)] = 1
    except:
      print("unable to open input moments file", h5input)
      return 0
  Marray = np.asarray(M)
  mdim = M.shape[0]
  nthreads = 10
  hMabove = []
  hMbelow = []
  for i in range(nthreads):
    hMabovediag = ROOT.TH1D(f"hMabovediag{i}", "M above-diagonal elements normalized to diagonal", 200,-0.1, 0.1)
    hMbelowdiag = ROOT.TH1D(f"hMbelowdiag{i}", "M below-diagonal elements normalized to diagonal", 200,-0.1, 0.1)
    hMabovediag.GetXaxis().SetTitle("normalized M matrix element")
    hMabovediag.GetYaxis().SetTitle("moments")
    hMbelowdiag.GetXaxis().SetTitle("normalized M matrix element")
    hMbelowdiag.GetYaxis().SetTitle("moments")
    hMabove.append(hMabovediag)
    hMbelow.append(hMbelowdiag)
  #Execute this line once at module load time
  #ROOT.gROOT.ProcessLine(".L fill_M_histograms.C+O")
  ROOT.fill_M_histograms.__release_gil__ = False
  ithread = 0
  threads = []
  for m in range(mdim):
    args = (Marray, m, mdim, hMabove[ithread], hMbelow[ithread])
    t = threading.Thread(target=ROOT.fill_M_histograms, args=args)
    if ithread < len(threads):
      ierr = threads[ithread].join()
      threads[ithread] = t
    else:
      threads.append(t)
    t.start()
    print(f"started row {m}", end='\r', flush=True)
    ithread = (ithread + 1) % nthreads
  hMabovediag = hMabove[0]
  hMbelowdiag = hMbelow[0]
  for i in range(1,nthreads):
    hMabovediag.Add(hMabove[i])
    hMbelowdiag.Add(hMbelow[i])
  return (hMabovediag,hMbelowdiag)

def histogram_offdiagonal(name, M):
  """
  Variant of histogram_M_offdiagonal() where the matrix M is
  supplied by argument rather than read from an hdf5 input file.
  """
  hMabovediag = ROOT.TH1D(f"hMabovediag{name}", "M above-diagonal elements normalized to diagonal", 10000,-0.1, 0.1)
  hMbelowdiag = ROOT.TH1D(f"hMbelowdiag{name}", "M below-diagonal elements normalized to diagonal", 10000,-0.1, 0.1)
  hMabovediag.GetXaxis().SetTitle("normalized M matrix element")
  hMabovediag.GetYaxis().SetTitle("moments")
  hMbelowdiag.GetXaxis().SetTitle("normalized M matrix element")
  hMbelowdiag.GetYaxis().SetTitle("moments")
  mdim = M.shape[0]
  for i in range(mdim):
    for j in range(i):
      hMbelowdiag.Fill(M[i,j] / (M[i,i] * M[j,j])**0.5)
    for j in range(i+1,mdim):
      hMabovediag.Fill(M[i,j] / (M[i,i] * M[j,j])**0.5)
  return (hMabovediag, hMbelowdiag)

def histogram_eigenvalues(name, title, w, normfactor=1):
  """
  Convenience function to form a histogram of a vector w, interpreted
  as an ordered list of the eigenvalues of the M matrix. Return value
  is a TH1D object with name,title given in the arguments. Optional
  argument normfactor is multiplied by w when it is histogrammed.
  """
  ndim = len(w)
  h = ROOT.TH1D(name, title, ndim, 0, ndim)
  for i in range(ndim):
    h.SetBinContent(ndim - i, w[i] * normfactor)
  h.GetXaxis().SetTitle("eigenvalue index")
  h.GetYaxis().SetTitle("acceptance eigenvalue")
  return h

def histogram_moments(name, title, moments, errors):
  """
  Convenience function to form a histogram of vector moments, treated
  as an ordered list of angular moments of a sample, with statistical
  uncertainty provided in the errors aregument. Return value is a
  new TH1D object with name,title given in the arguments.
  """
  ndim = len(moments)
  h = ROOT.TH1D(name, title, ndim, 0, ndim)
  for i in range(ndim):
    h.SetBinContent(i+1, moments[i])
    h.SetBinError(i+1, errors[i])
  h.GetXaxis().SetTitle("moment index")
  h.GetYaxis().SetTitle("sample moment")
  return h

def histograms_of_moments(Nmoments, basename, basetitle, 
                          atitles, nbins, alimits):
  """
  Returns a list of empty 2D histograms, one per moment 0..(Nmoments-1)
  with basename [basetitle] being string templates for the histogram
  names [titles] with embedded {0} strings subsituted for the moment
  index. Desired axis titles are supplied as atitles[0] for the x axis,
  atitles[1] for the y axis. Similarly, nbins[0] and nbins[1] are the
  number of histogram divisions along x and y, with alimits[0][0] and
  alimits[0][1] being the x axis limits, and alimits[1][0],alimits[1][1]
  being the corresponding limits along y.
  """
  histograms = []
  for i in range(Nmoments):
     name = basename.format(i)
     title = basetitle.format(i)
     h = ROOT.TH2D(name, title, nbins[0], alimits[0][0], alimits[0][1],
                                nbins[1], alimits[1][0], alimits[1][1])
     h.GetXaxis().SetTitle(atitles[0])
     h.GetYaxis().SetTitle(atitles[1])
     histograms.append(h)
  return histograms

def analyze_moments(massEtaPi0_limits=(0.6,2.5), abst_limits=(0.0,2.5), model=1,
                    sample_subset=range(5), acceptance_subset=range(5,10),
                    Mmatrix="Msaved.h5", use_generated_angles_for_acceptance=1,
                    use_generated_angles_for_sampling=1):
  """
  The primary functionality of the analyzeMomentsMatrix module is found within
  this one function. It carries out three tasks in sequence, as follows:
    1. generates a mock data sample on a single bin in (mX,t)
       a) scans over the Monte Carlo set devoted to sample creation
          and counts all of the generated events within the kinematic bin;
       b) uses the model weights from the scan in 1a to histogram the true
          values of the model moments within this kinematic bin, to be used
          later in a comparison with the acceptance-corrected moments;
       c) using the maximum model weight in the generated events scan, pass
          over the reconstructed events from the same generated sample and
          perform accept-reject to select the mock data sample within this
          kinematic bin.
    2. generate the M acceptance matrix on this same bin in (mX,t)
       a) compute the M matrix and its error on the unweighted Monte Carlo
          sample assigned to acceptance estimation, different from the one
          used to produce the mock data sample;
       b) diagonalize M to find the acceptance support moments vectors in
          this kinematic bin;
    3. map out the (mX,t) dependence of the acceptance eigenvalues within
       this kinematic bin for each support vector moment from step 2.
       a) histogram the acceptance in (mX,t) for each support vector moment
          as a sum over unweighted events in the acceptance subset;
       b) normalize the acceptance to the number of generated events in
          (mX,t) over the same kinematic bin;
  Even with multithreaded uproot processing on large-memory machine, the above
  steps take too long for interactive work because each call takes about 1 hour
  to complete, and it must be repeated for each kinematic bin. The shell script
  rebuildMoments.sh has been provided in the scripts directory to automate the
  original execution of analyze_moments over the full set of kinematic bins.
  After this has been done once, analyze_moments saves the results in the work
  directory so that subsequent calls to analyze_moments are fast.

  Input arguments:
    massEtaPi0_limits, abst_limits: defines the kinematic bin to be analyzed
    model: normally should be =1 for operation as defined above, if =0 then
           assign uniform weight to all events, simulating a model with H00=1
           and all other moments set to 0, useful for testing software.
    sample_subset, acceptance_subset: a list of indices in the range [0..9]
           assigning portions of the original 100 million event Monte Carlo
           dataset to formation of the mock data sample vs estimation of the
           experimental acceptance, in blocks of 10 million events.
    Mmatrix: don't assign this, expert use only.

  Return value: a 4-tuple consisting of
    [0]: TH2D histogram containing the experimental moments sum over the
         mock data sample, without any acceptance correction, as a function
         of (mX,t) within this kinematic bin;
    [1]: TH2D histogram containing the acceptance-corrected moments derived
         from the mock data sample by support-vector acceptance correction,
         as a function of (mX,t) within this kinematic bin;
    [2]: TH2D histogram of the generated model moments computed directly 
         from the generated event sample weighted by the trial model, scaled
         by maxweight for direct comparison with [1];
    [3]: TH2D difference between [1] and [2].
  """
  kinbins = (190, 25)
  kinbounds = ((0.6,2.5), (0.0,2.5)) # mass in GeV, |t| in GeV^2

  try:
    f5 = h5py.File("Msample.h5")
    samplemom = f5['sample_moment'][:]
    refermom = f5['reference_moment'][:]
    samplecov = f5['sample_covariance'][:]
    refercov = f5['reference_covariance'][:]
    mock_events = f5['sample_events'][:]
    gen_events = f5['generated_events'][:]
    maxweight = f5['sample_maxweight'][()]
  except:
    try:
      f5.close()
    except:
      pass
    print("hinhout regenerating Msample.h5")
    mock_events = []
    gen_events = []
    h2dmodel = histograms_of_moments(169, "model1_{0}", "model 1 moment {0}", 
                                     ("massEtaPi0 (GeV)", "|t| (GeV^2)"),
                                     kinbins, kinbounds)
    h2dsample = histograms_of_moments(169, "sample_{0}", "mock sample moment {0}", 
                                      ("massEtaPi0 (GeV)", "|t| (GeV^2)"),
                                      kinbins, kinbounds)
    for i in sample_subset:
      tstart1 = time.perf_counter()
      bmm.open(f"{intreedir}/generated_moments_x10_{i}.root:etapi0_moments")
      events,weights = bmm.select_events(massEtaPi0_limits=massEtaPi0_limits,
                                         abst_limits=abst_limits, model=1)
      gen_events += events
      maxweight = max(weights) * 1.2
      tstep1a = time.perf_counter()
      print(f"  select generated events: {tstep1a-tstart1:.3f}s")
      mom,cov = bmm.compute_moments(events, mPi0=1, mEta=1, weights=weights)
      try:
        refermom += mom
        refercov += cov
      except:
        refermom = mom
        refercov = cov
      tstep1b = time.perf_counter()
      print(f"  compute generated moments: {tstep1b-tstep1a:.3f}s")
      bmm.histogram_moments(events, h2dmodel, mPi0=1, mEta=1, weights=weights)
      tstep1c = time.perf_counter()
      print(f"  histogram generated moments: {tstep1c-tstep1b:.3f}s")
      bmm.open(f"{intreedir}/etapi0_moments_x10_{i}.root:etapi0_moments")
      events,weights = bmm.select_events(massEtaPi0_limits=massEtaPi0_limits,
                                         abst_limits=abst_limits, model=1,
                       use_generated_angles=use_generated_angles_for_sampling,
                                         maxweight=maxweight)
      mock_events += events
      tstep1d = time.perf_counter()
      print(f"  select mock data sample: {tstep1d-tstep1c:.3f}s")
      mom,cov = bmm.compute_moments(events, mPi0=1, mEta=1)
      try:
        samplemom += mom
        samplecov += cov
      except:
        samplemom = mom
        samplecov = cov
      tstep1e = time.perf_counter()
      print(f"  compute mock data moments: {tstep1e-tstep1d:.3f}s")
      bmm.histogram_moments(events, h2dsample, mPi0=1, mEta=1)
      tstop1 = time.perf_counter()
      print(f"  histogram mock data moments: {tstop1-tstep1e:.3f}s")
      print(f" *total sample creation time: {tstop1-tstart1:.3f}s")
    f5 = h5py.File(f"{savedir}/Msample.h5", "w")
    f5.create_dataset("sample_moment", data=samplemom)
    f5.create_dataset("reference_moment", data=refermom)
    f5.create_dataset("sample_covariance", data=samplecov)
    f5.create_dataset("reference_covariance", data=refercov)
    f5.create_dataset("sample_events", data=mock_events)
    f5.create_dataset("generated_events", data=gen_events)
    f5.create_dataset("sample_maxweight", data=maxweight)
    f5.close()
    f = ROOT.TFile(f"{savedir}/Msample.root", "recreate")
    [h.Write() for h in h2dmodel]
    [h.Write() for h in h2dsample]
    f.Close()

  try:
    f5 = h5py.File("Msaved.h5")
    M = f5['Moments'][:]
    Nacc = f5['accepted_subset'][()]
    Ngen = f5['generated_subset'][()]
  except:
    try:
      f5.close()
    except:
      pass
    print("hinhout regenerating Msaved.h5")
    f5 = h5py.File("Msaved.h5", 'w')
    acc_events = []
    for i in acceptance_subset:
      tstart2 = time.perf_counter()
      bmm.open(f"{intreedir}/etapi0_moments_x10_{i}.root:etapi0_moments")
      events,weights = bmm.select_events(massEtaPi0_limits=massEtaPi0_limits,
                                         abst_limits=abst_limits)
      f5.create_dataset(f"accepted_subset[{i}]", data=len(events))
      M_ = bmm.buildMomentsMatrix_threaded(events, mPi0=1, mEta=1,
                     use_generated_angles=use_generated_angles_for_acceptance)
      f5.create_dataset(f"Moments[{i}]", data=M_)
      acc_events += events
      try:
        M += M_
      except:
        M = M_
      tstop2 = time.perf_counter()
      print(f"  time to compute M matrix: {tstop2-tstart2:.3f}s")

      # disable this, unless you need to check the orthonormality of the moments
      """
      bmm.open(f"{intreedir}/generated_moments_x10_{i}.root:etapi0_moments")
      events,weights = bmm.select_events(massEtaPi0_limits=massEtaPi0_limits,
                                         abst_limits=abst_limits)
      f5.create_dataset(f"generated_subset[{i}]", data=len(events))
      M_ = bmm.buildMomentsMatrix_threaded(events, mPi0=1, mEta=1,
                     use_generated_angles=use_generated_angles_for_acceptance)
      f5.create_dataset(f"Mgenerated[{i}]", data=M_)
      try:
        Mgen += M_
      except:
        Mgen = M_
      """

    W,E,Vt = np.linalg.svd(M)

    h2daccepted = histograms_of_moments(169, "accept_{0}",
                                        "acceptance for support vector {0}", 
                                        ("massEtaPi0 (GeV)", "|t| (GeV^2)"),
                                        kinbins, kinbounds)
    h2dgenerated = histograms_of_moments(1, "generated", "generated spectrum", 
                                         ("massEtaPi0 (GeV)", "|t| (GeV^2)"),
                                         kinbins, kinbounds)
    gen_events = []
    for i in acceptance_subset:
      tstart3 = time.perf_counter()
      bmm.open(f"{intreedir}/etapi0_moments_x10_{i}.root:etapi0_moments")
      events,weights = bmm.select_events(massEtaPi0_limits=massEtaPi0_limits,
                                         abst_limits=abst_limits)
      bmm.histogram_acceptance(events, h2daccepted, mPi0=1, mEta=1, svectors=W,
                    use_generated_angles=use_generated_angles_for_acceptance)
      tstep3a = time.perf_counter()
      print(f"  time to histogram acceptance moments: {tstep3a-tstart3:.3f}s")
      bmm.open(f"{intreedir}/generated_moments_x10_{i}.root:etapi0_moments")
      events,weights = bmm.select_events(massEtaPi0_limits=massEtaPi0_limits,
                                         abst_limits=abst_limits)
      f5.create_dataset(f"generated_subset[{i}]", data=len(events))
      gen_events += events
      bmm.histogram_acceptance(events, h2dgenerated, mPi0=1, mEta=1, mGJ=1)
      tstop3 = time.perf_counter()
      print(f"  time to histogram generated kinematics: {tstop3-tstep3a:.3f}s")

    f5.create_dataset("Moments", data=M)
    f5.create_dataset("accepted_subset", data=len(acc_events))
    f5.create_dataset("generated_subset", data=len(gen_events))
    f5.close()

    Nacc = len(acc_events)
    Ngen = len(gen_events)
    f = ROOT.TFile(f"{savedir}/Msaved.root", "recreate")
    [h.Write() for h in h2daccepted]
    [h.Write() for h in h2dgenerated]
    f.Close()

  M *= (4 * math.pi)**3 / Ngen

  if Mmatrix != "Msaved.h5":
    f5 = h5py.File(Mmatrix)
    Mregion = f5['Moments']
  else:
    Mregion = M
  W,E,Vt = np.linalg.svd(Mregion)

  Minv = np.linalg.inv(M)
  correctmom = Minv @ samplemom
  correctcov = Minv @ samplecov @ Minv.T

  # transform from spherical basis to support vector moments
  if False:
    samplemom = W.T @ samplemom
    samplecov = W.T @ samplecov @ W
    correctmom = W.T @ correctmom
    correctcov = W.T @ correctcov @ W
    refermom = Vt @ refermom @ Vt.T
    refercov = Vt @ refercov @ Vt.T

  Nmom = len(samplemom)
  hsamp = ROOT.TH1D("hsamp", "sample moments, uncorrected", Nmom, 0, Nmom)
  hsamp.GetXaxis().SetTitle("moment index")
  hsamp.GetYaxis().SetTitle("sample moment")
  hcorr = ROOT.TH1D("hcorr", "sample moments, acceptance corrected", Nmom, 0, Nmom)
  hcorr.GetXaxis().SetTitle("moment index")
  hcorr.GetYaxis().SetTitle("sample moment")
  hmodel = ROOT.TH1D("hmodel", "model moments, before acceptance", Nmom, 0, Nmom)
  hmodel.GetXaxis().SetTitle("moment index")
  hmodel.GetYaxis().SetTitle("model moment")
  hdiff = ROOT.TH1D("hdiff", "corrected sample - model moments difference", Nmom, 0, Nmom)
  hdiff.GetXaxis().SetTitle("moment index")
  hdiff.GetYaxis().SetTitle("moment (corrected_sample - model")
  for i in range(Nmom):
    hsamp.SetBinContent(i+1, samplemom[i])
    hsamp.SetBinError(i+1, samplecov[i,i]**0.5)
    hcorr.SetBinContent(i+1, correctmom[i])
    hcorr.SetBinError(i+1, correctcov[i,i]**0.5)
    hmodel.SetBinContent(i+1, refermom[i])
    hmodel.SetBinError(i+1, refercov[i,i]**0.5)
    hdiff.SetBinContent(i+1, correctmom[i])
    hdiff.SetBinError(i+1, correctcov[i,i]**0.5)
  hmodel.Scale(1 / maxweight)
  hdiff.Add(hmodel, -1)
  return hsamp,hcorr,hmodel,hdiff

def standard_kinematic_bins(finebins=0):
  """
  Returns a standard subdivision of the (mX,t) plane into kinematic bins,
  according to a user-selected scheme. This can be extended as desired,
  but doing so will require a new pass over analyze_moments() for the
  new set of kinematic bins. Current choices are:
    finebins=0: 24 coarse bins, fastest option, sufficient for many things;
    finebins=1: 48 finer bins, finer segmentation in mX, same in t;
    finebins=2: 144 finer bins, finer segmentation in both mX and t;
  """
  if finebins == 2:
    tbins = ((0.0,0.2), (0.2,0.3), (0.3,0.45), (0.45,0.6),
             (0.6,0.9), (0.9,1.2), (1.2,1.6), (1.6,2.4))
    mbins = ((0.6,0.7), (0.7,0.8), (0.8,0.9), (0.9,1.0),
             (1.0,1.1), (1.1,1.2), (1.2,1.3), (1.3,1.4),
             (1.4,1.5), (1.5,1.6), (1.6,1.7), (1.7,1.8),
             (1.8,1.9), (1.9,2.0), (2.0,2.1), (2.1,2.2),
             (2.2,2.3), (2.3,2.4))
  elif finebins == 1:
    tbins = ((0.0,0.3), (0.3,0.6), (0.6,1.2), (1.2,2.5))
    mbins = ((0.6,0.75), (0.75,0.9), (0.9,1.05),
             (1.05,1.2), (1.2,1.35), (1.35,1.5),
             (1.5,1.65), (1.65,1.8), (1.8,1.95),
             (1.95,2.1), (2.1,2.25), (2.25,2.4))
  else:
    tbins = ((0.0,0.3), (0.3,0.6), (0.6,1.2), (1.2,2.5))
    mbins = ((0.6,0.9), (0.9,1.2), (1.2,1.5),
             (1.5,1.8), (1.8,2.1), (2.1,2.5))
  kinbins = []
  for tbin in tbins:
    for mbin in mbins:
      kinbins.append((tbin,mbin))
  return kinbins

def model1_corrected_moments(imoments=range(169), kinbins=[], finebins=0):
  """
  This function takes What analyze_moments() does for a single kinematic bin
  in (mX,t), and applies it to the full range kinematic range covered by the
  standard kinematic binning scheme finebins. It reads from the saved hdf5
  files in the standard bin directories that it assumes are located under the
  work directory, and assembles the complete results for the mock data sample,
  both with and without acceptance correction, together with the calculated
  moments from the trial model straight from the generator. The returned
  result is a 4-tuple consising of the following lists, each one of the
  same length as the imoments input list.
   [0]: TH2D histograms on (mX,t) of the sample moments listed in imoments,
        prior to acceptance correction;
   [1]: TH2D histograms on (mX,t) of the sample moments listed in imoments
        after the support-vector acceptance correction has been applied;
   [2]: TH2D histograms on (mX,t) of the trial model moments straight from
        the generator, averaged over the same binning scheme as [0] and [1];
   [3]: TH2D histogram on (mX,t) of the number of generated events in the
        sample used to estimate the acceptance, for acceptance normalization.
  """
  hsample = [0 for i in imoments]
  hcorrect = [0 for i in imoments]
  hmodel1 = [0 for i in imoments]
  if len(kinbins) == 0:
    kinbins = standard_kinematic_bins(finebins)
  for tbin,mbin in kinbins:
    datadir = f"{workdir}/etapi0_moments_{mbin[0]},{mbin[1]}_{tbin[0]},{tbin[1]}"
    f5saved = h5py.File(datadir + "/Msaved.h5")
    Ngen = f5saved['generated_subset'][()]
    M = f5saved['Moments'][:]
    M *= (4 * math.pi)**3 / Ngen
    W,E,Vt = np.linalg.svd(M)
    Minv = np.linalg.inv(M)
    Nmoments = M.shape[0]
    fsaved = ROOT.TFile(datadir + "/Msaved.root")
    h2dgenerated = fsaved.Get(f"generated")
    f5sample = h5py.File(datadir + "/Msample.h5")
    maxweight = f5sample['sample_maxweight'][()]
    fsample = ROOT.TFile(datadir + "/Msample.root")
    h2dsample = [fsample.Get(f"sample_{k}") for k in range(Nmoments)]
    svector = f"sample_{mbin}_{tbin}"
    support_rootfile = datadir + "/Msupport.root"
    try:
      fsupport = ROOT.TFile(support_rootfile)
      for k in range(Nmoments):
        svk = f"{svector}_s{k}"
        h2dsupport[svk] = fsupport.Get(svk)
        h2dsupport[svk].SetDirectory(0)
      fsupport.Close()
    except:
      print("building support vectors for", svector)
      fsupport = ROOT.TFile(support_rootfile, "recreate")
      for k in range(Nmoments):
        svk = f"{svector}_s{k}"
        h2dsupport[svk] = h2dsample[k].Clone(svk)
        h2dsupport[svk].SetTitle(f"sample support vector moment {k}")
        h2dsupport[svk].Reset()
      for k in range(Nmoments):
        svk = f"{svector}_s{k}"
        for j in range(Nmoments):
          h2dsupport[svk].Add(h2dsample[j], W[j,k])
        h2daccept = fsaved.Get(f"accept_{k}")
        h2daccept.Divide(h2dgenerated)
        h2dsupport[svk].Divide(h2daccept)
        h2dsupport[svk].Write()
        h2dsupport[svk].SetDirectory(0)
      fsupport.Close()
    for imoment in imoments:
      h2dmodel1 = fsample.Get(f"model1_{imoment}")
      h2dcorrect = h2dsample[imoment].Clone(f"correct_{imoment}")
      h2dcorrect.Reset()
      for k in range(Nmoments):
        if False:
          h2dcorrect.Add(h2dsample[k], Minv[imoment,k])
        else:
          svk = f"{svector}_s{k}"
          h2dcorrect.Add(h2dsupport[svk], Vt[k,imoment])
      try:
        hmodel1[imoment].Add(h2dmodel1)
        hsample[imoment].Add(h2dsample[imoment], maxweight)
        hcorrect[imoment].Add(h2dcorrect, maxweight)
      except:
        hmodel1[imoment] = h2dmodel1.Clone(f"hmodel1_{imoment}")
        hsample[imoment] = h2dsample[imoment].Clone(f"hsample_{imoment}")
        hcorrect[imoment] = h2dcorrect.Clone(f"hcorrect_{imoment}")
        hsample[imoment].Scale(maxweight)
        hcorrect[imoment].Scale(maxweight)
        hmodel1[imoment].SetDirectory(0)
        hsample[imoment].SetDirectory(0)
        hcorrect[imoment].SetDirectory(0)
    try:
      hgenerated.Add(h2dgenerated)
    except:
      hgenerated = h2dgenerated.Clone(f"hgenerated")
      hgenerated.SetDirectory(0)
  return hsample,hcorrect,hmodel1,hgenerated

def apply_constraints(S, moments, covariance, histograms=[]):
  """
  Takes in a set of acceptance-corrected sample moments together with
  their covariance matrix, and applies the constraints contained in
  argument S, and returns the updated moments and covariance matrix
  as the first two elements of a 3-tuple. If argument histograms is
  supplied in the argument list then it must have the same length
  as moments and as the dimension of square matrix covariance, and
  represent the acceptance-corrected sample moments as a TH1D or TH2D.
  Matrix constraint S is in the form of a linear condition among the
  moments vector h, expressed as S h = 0. The number of columns must
  be equal to the number of moments, while the number of rows is the
  number of constraint equations. If argument histograms is present
  in the input argument list, then the updated histograms following
  application of the constraints will be returned in the third 
  element of the return 3-tuple, otherwise an empty list.
  """
  C = covariance
  D = C @ S.T @ np.linalg.inv(S @ C @ S.T) @ S
  if D.shape != covariance.shape:
    print("error in apply_constraints - shape of S does not match moments")
    return S, moments, covariance
  h = moments - D @ moments
  CC = C - D @ C - C @ D.T + D @ C @ D.T
  histograms_cons = []
  Nmoments = len(moments)
  if len(histograms) > 0:
    for hist in histograms:
      name = hist.GetName() + "_cons"
      title = "constrained " + hist.GetTitle()
      histograms_cons.append(hist.Clone(name))
      histograms_cons[-1].SetTitle(title)
    for k in range(Nmoments):
      for j in range(Nmoments):
        histograms_cons[k].Add(histograms[k], -D[k,j])
  return h, CC, histograms_cons

def scan_em(corrected=1, scale=1, tcut=0, finebins=0, 
            hchisq=0, constraints=[], hchisq2=0, maxmoments=0,
            acceptance_subset=[], interactive=True):
  """
  Cycle through the kinematic bins and perform a systematic comparison
  between the acceptance-corrected moments and those from the original
  trial model. The distribution of chisquare values from a projection
  of each of the moments onto the mX axis for a given slice through |t|
  is returned.

  This method has been tweaked repeatedly to add the capability of studying
  systematic errors from a number of sources. These features are enabled by
  calling with a non-empty constraints array, even if all of the elements
  are zero (implies no constraints). In each case, the results are returned
  as a second chisquare histogram in the second element of the return value.
    1) acceptance_subset - select a subset of the MC event sample that was
            passed to analyze_moments in a prior call. This provides a hook
            for study of systematic error from a finite acceptance sample.
    2) maxmoments [=Q] - if non-zero, truncate the M matrix to QxQ before its
            inversion, and extract only Q moments from the sample. It must
            not be set to a value greater than 169.
  """
  if hchisq == 0:
    title = "#chi^{2} distribution for corrected-model moments, 95 dof"
    hchisq = ROOT.TH1D("hchisq", title, 100, 0, 200)
    hchisq.GetXaxis().SetTitle("#chi^{2}")
    hchisq.GetYaxis().SetTitle("moments")
    hchisq.SetDirectory(0)
  if hchisq2 == 0:
    title = "#chi^{2} distribution for constrained moments"
    hchisq2 = ROOT.TH1D("hchisq2", title, 100, 0, 500)
    hchisq2.GetXaxis().SetTitle("#chi^{2}")
    hchisq2.GetYaxis().SetTitle("kinematic bins")
    hchisq2.SetDirectory(0)
  hchisq3 = hchisq2.Clone("hchisq3")
  hchisq3.SetTitle("#chi^{2} distribution of constraints")
  hchisq3.SetDirectory(0)
  if maxmoments == 0:
    maxmoments = 169
  hsample = [0] * maxmoments
  hcorrect = [0] * maxmoments
  hmodel1 = [0] * maxmoments
  prompt = interactive
  itbin,tlimits = (0,(0,0))
  for tbin,mbin in standard_kinematic_bins(finebins):
    if tbin != tlimits:
      tlimits = tbin
      itbin += 1
      if itbin == tcut:
        break
  for tbin,mbin in standard_kinematic_bins(finebins):
    if tcut != 0 and tbin != tlimits:
      continue
    hs,hc,hm,hg = model1_corrected_moments(range(maxmoments),
                                           kinbins=[(tbin,mbin)])
    if len(constraints) > 0:
      datadir = f"{workdir}/etapi0_moments_{mbin[0]},{mbin[1]}_{tbin[0]},{tbin[1]}"
      f5saved = h5py.File(datadir + "/Msaved.h5")
      if len(acceptance_subset) > 0:
        M = 0
        Ngen = 0
        for i in acceptance_subset:
          M_ = f5saved[f"Moments[{i}]"][:maxmoments,:maxmoments]
          Ngen_ = f5saved[f"generated_subset[{i}]"][()]
          try:
            M += M_
            Ngen += Ngen_
          except:
            M = M_
            Ngen = Ngen_
      else:
        M = f5saved["Moments"][:maxmoments,:maxmoments]
        Ngen = f5saved['generated_subset'][()]
      M *= (4 * math.pi)**3 / Ngen
      Minv = np.linalg.inv(M)
      f5sample = h5py.File(datadir + "/Msample.h5")
      samplemom = f5sample['sample_moment'][:maxmoments]
      samplecov = f5sample['sample_covariance'][:maxmoments,:maxmoments]
      maxweight = f5sample['sample_maxweight'][()] 
      refermom = f5sample['reference_moment'][:maxmoments]
      refercov = f5sample['reference_covariance'][:maxmoments,:maxmoments]
      correctmom = Minv @ samplemom
      correctcov = Minv @ samplecov @ Minv.T
      if np.max(constraints) > np.min(constraints):
        cmom,ccov,chc = apply_constraints(constraints, correctmom, correctcov, hc)
      else:
        cmom,ccov,chc = (correctmom, correctcov, 0)
      title = f"moments for |t|={tbin}, m_{{X}}={mbin}"
      h1 = ROOT.TH1D("h1", "constrained " + title, maxmoments, 0, maxmoments)
      h1.GetXaxis().SetTitle("moments index")
      h1.GetYaxis().SetTitle("sample moment")
      h1.SetStats(0)
      h2 = ROOT.TH1D("h2", "model 1 " + title, maxmoments, 0, maxmoments)
      h2.GetXaxis().SetTitle("moments index")
      h2.GetYaxis().SetTitle("model moment")
      h2.SetStats(0)
      for i in range(maxmoments):
        h1.SetBinContent(i+1, cmom[i])
        h1.SetBinError(i+1, ccov[i,i]**0.5)
        h2.SetBinContent(i+1, refermom[i] / maxweight)
        h2.SetBinError(i+1, refercov[i,i]**0.5 / maxweight)
      diff_mom = cmom - refermom / maxweight
      netcov = ccov + refercov / maxweight**2
      chi2 = diff_mom @ np.linalg.inv(netcov) @ diff_mom
      ndof = len(diff_mom)
      hchisq2.Fill(chi2)
      diff_cons = cmom - correctmom
      chi2_cons = diff_cons @ np.linalg.inv(correctcov) @ diff_cons
      hchisq3.Fill(chi2_cons)
      h2.SetLineColor(2)
      h1.Draw()
      h2.Draw("same")
      ROOT.gROOT.FindObject("c1").Update()
      if interactive:
        print(f"chi2 for bin {mbin},{tbin} is {chi2:.1f} / {ndof}")
        if prompt != 'g':
          prompt = input("g to continue without prompt, " +
                         "q to quit? ")
          if prompt == 'q':
             return 0,0
    for m in range(maxmoments):
      try:
        hsample[m].Add(hs[m])
        hcorrect[m].Add(hc[m])
        hmodel1[m].Add(hm[m])
      except:
        hsample[m] = hs[m].Clone()
        hcorrect[m] = hc[m].Clone()
        hmodel1[m] = hm[m].Clone()
  prompt = interactive
  for m in range(maxmoments):
    if corrected:
      h = hcorrect[m]
    else:
      h = hsample[m]
    px1 = h.GetName() + "_px"
    h1 = h.ProjectionX(px1)
    h1.Scale(scale)
    h1.Rebin(2)
    px2 = hmodel1[m].GetName() + "_px"
    h2 = hmodel1[m].ProjectionX(px2)
    h2.SetLineColor(2)
    h2.Rebin(2)
    hmax = 0
    hmin = 0
    for i in range(h1.GetNbinsX()):
      y1 = h1.GetBinContent(i+1)
      e1 = h1.GetBinError(i+1)
      y2 = h2.GetBinContent(i+1)
      e2 = h2.GetBinError(i+1)
      if e1 < 1e5:
        hmax = max(hmax, y1 + e1)
        hmin = min(hmin, y1 - e1)
      if e2 < 1e5:
        hmax = max(hmax, y2 + e2)
        hmin = min(hmin, y2 - e2)
    h1.SetMaximum(hmax + (hmax - hmin) * 0.1)
    h1.SetMinimum(hmin + (hmin - hmax) * 0.1)
    chisq = 0
    ndof = 0
    for n in range(h1.GetNbinsX()):
       y1 = h1.GetBinContent(n+1)
       e1 = h1.GetBinError(n+1)
       y2 = h2.GetBinContent(n+1)
       e2 = h2.GetBinError(n+1)
       chisq += (y1 - y2)**2 / (e1**2 + e2**2 + 1e-99)
       ndof += 1
    hchisq.Fill(chisq)
    if len(h1.GetYaxis().GetTitle()) > 0:
      title = h1.GetYaxis().GetTitle()
    else:
      title = h1.GetTitle()
      h1.GetYaxis().SetTitle(title)
    h1.SetTitle(f"{title}, tcut={tcut}, #chi^{{2}} = {chisq:.1f} / {ndof}")
    h1.SetStats(0)
    h1.Draw()
    h2.Draw("same")
    c1 = ROOT.gROOT.FindObject("c1")
    c1.Update()
    if interactive:
      print(f"moment {m}: chisquare = {chisq}, ndof={ndof}")
      if prompt != 'g':
        prompt = input("<enter> to continue, " +
                       "g to continue without prompt, " +
                       "q to quit: ")
        if prompt == 'q':
          break
  return hchisq,hchisq2,hchisq3

def histogram_moments_correlations(support_moments=0):
  """
  Compute a statistical measure of the degree of correlation between the
  acceptance-corrected sample moments in the standard spherical basis
  (support_moments=0) or the support moments basis (support_moments=1).
  Return value is a histogram that reports the answer as a distribution
  of the product of the diagonal elements of the covariance matrix C[i,i]
  times the corresponding diagonal element of its inverse (1/C)[i,i]. 
  For an uncorrelated set of moments, all such products are 1, otherwise
  they are less than 1. Deviation from unity is a measure of the degree
  to which the moments can be treated as more or less uncorrelated.
  """
  hcorr = {}
  for tbin,mbin in standard_kinematic_bins():
    datadir = f"{workdir}/etapi0_moments_{mbin[0]},{mbin[1]}_{tbin[0]},{tbin[1]}"
    f5sample = h5py.File(datadir + "/Msample.h5")
    cov = f5sample["sample_covariance"][:]
    Nmoments = cov.shape[0]
    name = f"hcorr_{mbin[0]},{mbin[1]}_{tbin[0]},{tbin[1]}"
    if support_moments:
      title = f"net support moment correlation, m_{{X}}={mbin}, |t|={tbin}"
    else:
      title = f"net spherical moment correlation, m_{{X}}={mbin}, |t|={tbin}"
    hcorr[name] = ROOT.TH1D(name, title, Nmoments, 0, Nmoments)
    hcorr[name].GetXaxis().SetTitle("moments index")
    hcorr[name].GetYaxis().SetTitle("correlation parameter")
    hcorr[name].SetStats(0)
    if support_moments:
      f5saved = h5py.File(datadir + "/Msaved.h5")
      M = f5saved['Moments'][:]
      W,E,Vt = np.linalg.svd(M)
      cov = W.T @ cov @ W
    covinv = np.linalg.inv(cov)
    for i in range(Nmoments):
      hcorr[name].SetBinContent(i+1, covinv[i,i] * cov[i,i])
  return hcorr

def make_Kmatrix(Lmax=3):
  """
  Computes and returns the matrix K that gives the linear combination
  of density matrix elements corresponding to each moment. The density
  matrix is represented by (Lmax+1)**4 real values (see the r vector
  definition at the head of rank1_density_from_moments), which the
  K matrix maps onto the (2*Lmax+1)**2 moments as a linear transform.
  """
  Nmoments = (2*Lmax + 1)**2
  Nquantum = (Lmax + 1)**2
  K = np.zeros([Nmoments, Nquantum**2])
  a = -1
  for La in range(2*Lmax+1):
    for Ma in range(-La, La+1):
      a += 1
      i = -1
      for L1 in range(Lmax+1):
        for M1 in range(-L1, L1+1):
          i += 1
          j = -1
          for L2 in range(Lmax+1):
            for M2 in range(-L2, L2+1):
              j += 1
              g = 0
              if i == j:
                g = (ROOT.Math.wigner_3j(2*L1, 2*L2, 2*La,
                                         2*M1, -2*M2, -2*Ma)
                   * ROOT.Math.wigner_3j(2*L1, 2*L2, 2*La,
                                         0,    0,    0)
                   * (-1)**M2)
              elif Ma == 0 and i < j:
                  g = (ROOT.Math.wigner_3j(2*L1, 2*L2, 2*La,
                                           2*M1, -2*M2, -2*Ma)
                     * ROOT.Math.wigner_3j(2*L1, 2*L2, 2*La,
                                           0,    0,    0)
                     * (-1)**M2
                     + ROOT.Math.wigner_3j(2*L2, 2*L1, 2*La,
                                           2*M2, -2*M1, -2*Ma)
                     * ROOT.Math.wigner_3j(2*L1, 2*L2, 2*La,
                                           0,    0,    0)
                     * (-1)**M1) / 2**0.5
              elif Ma > 0 and i < j:
                  g = (ROOT.Math.wigner_3j(2*L1, 2*L2, 2*La,
                                           2*M1, -2*M2, -2*Ma)
                     * ROOT.Math.wigner_3j(2*L1, 2*L2, 2*La,
                                           0,    0,    0)
                     * (-1)**M2
                     + ROOT.Math.wigner_3j(2*L2, 2*L1, 2*La,
                                           2*M2, -2*M1, -2*Ma)
                     * ROOT.Math.wigner_3j(2*L1, 2*L2, 2*La,
                                           0,    0,    0)
                     * (-1)**M1)
              elif Ma < 0 and i > j:
                  g = (ROOT.Math.wigner_3j(2*L1, 2*L2, 2*La,
                                           2*M1, -2*M2, 2*Ma)
                     * ROOT.Math.wigner_3j(2*L1, 2*L2, 2*La,
                                           0,    0,    0)
                     * (-1)**M2
                     - ROOT.Math.wigner_3j(2*L2, 2*L1, 2*La,
                                           2*M2, -2*M1, 2*Ma)
                     * ROOT.Math.wigner_3j(2*L1, 2*L2, 2*La,
                                           0,    0,    0)
                     * (-1)**M1) * (-1)
              K[a, Nquantum*i+j] = g * ((2*L1+1) * (2*L2+1) * (2*La+1)
                                        / (4 * np.pi))**0.5
  return K

def make_Bmatrix(Lmax=3):
  """
  Computes and returns the matrix B that enforces the constraint that
  the density matrix is rank 1, which is written in terms of the rvector
  representation of the density matrix as r @ B @ r = 0.
  """
  Nquantum = (Lmax + 1)**2
  f = np.zeros([Nquantum**2], dtype=float)
  for i in range(Nquantum):
    f[i * (Nquantum + 1)] = 1
  B = (np.kron(f.reshape(Nquantum**2,1), f.reshape(1,Nquantum**2))
       - np.diag(np.ones([Nquantum**2], dtype=float)))
  return B

def get_model1_rvector(massEtaPi0=0, abst=0, rmatrix=[]):
  """
  Returns the r vector representation of the density matrix defined
  by model 1 at kinematics massEtaPi0,abst. The r vector is a compact
  representation of the n x n Hermetian density matrix where n=Lmax+1
  in terms of n**4 real components defined as follows. Basis quantum
  numbers (L,M) are listed in increasing order of allowed values, as
  L=0..Lmax, M=-L..L first incrementing M for fixed L, then increasing
  L, with list index LM = 0,1,2,...,n^2-1. Then, pairs of elements from
  this list are formed (LM,L'M') counting first over L'M', then over
  LM, with final list elements ((L,M),(L',M'))[i] uniquely labeling
  each component of the density matrix, i = 0,1,2,...n**4-1. The r[i]
  vector elements are defined as
              /
              | sqrt(2) |R| Re{ rho[LM(i),L'M'(i)] }, L'M' > LM
      r[i] = <  |R| rho[LM(i),L'M'(i)], L'M' == LM
              | sqrt(2) |R| Im{ rho[L'M'(i),LM(i)] }, L'M' < LM
              \
  where |R| is the trace of the rho matrix. If rmatrix is given then
  massEtaPi0,abst are ignored, and the information is taken from the
  density matrix specified by rmatrix.
  """
  if len(rmatrix) > 0:
    N = rmatrix.shape[0]
    rvector = np.empty([N**2], dtype=float)
    for i in range(N):
      rvector[i*N + i] = np.real(rmatrix[i,i])
      for j in range(i+1, N):
        rvector[i*N + j] = np.real(rmatrix[i,j]) * 2**0.5
        rvector[j*N + i] = np.imag(rmatrix[i,j]) * 2**0.5
    return rvector
  try:
    model1 = ROOT.trial_model1()
  except:
    ROOT.gROOT.ProcessLine(".L src/trial_model.C+O")
    model1 = ROOT.trial_model1()
  Lmax = int(model1.amplitude_Lmax())
  rvector = np.empty([(Lmax + 1)**4], dtype=float)
  sqrt2 = 2**0.5
  a = -1
  i = -1
  for L1 in range(Lmax+1):
    for M1 in range(-L1,L1+1):
      i += 1
      amp1 = model1.amplitude(L1, M1, massEtaPi0, abst)
      j = -1
      for L2 in range(Lmax+1):
        for M2 in range(-L2,L2+1):
          j += 1
          amp2 = model1.amplitude(L2, M2, massEtaPi0, abst)
          a += 1
          if i < j:
            rvector[a] = sqrt2 * (amp1[0] * amp2[0] + amp1[1] * amp2[1])
          elif i == j:
            rvector[a] = amp1[0] * amp2[0] + amp1[1] * amp2[1]
          else:
            rvector[a] = sqrt2 * (amp1[1] * amp2[0] - amp1[0] * amp2[1])
  return rvector

def get_model1_rmatrix(massEtaPi0=0, abst=0, rvector=[]):
  """
  Returns the complex density matrix defined by model 1 at kinematics
  massEtaPi0,abst, either by reading it from the rvector or by computing
  it directly from the quantum amplitudes if rvector=[]. If rvector is
  given then arguments massEtaPi0,abst are ignored. The rmatrix returned
  is converted to a density matrix by dividing by its trace.
  """
  try:
    model1 = ROOT.trial_model1()
  except:
    ROOT.gROOT.ProcessLine(".L src/trial_model.C+O")
    model1 = ROOT.trial_model1()
  Lmax = int(model1.amplitude_Lmax())
  nbasis = (Lmax + 1)**2
  sqrt2 = 2**0.5
  if len(rvector) > 0:
    rmatrix = np.empty([nbasis, nbasis], dtype=complex)
    for i in range(nbasis):
      rmatrix[i,i] = rvector[i * (nbasis + 1)]
      for j in range(i+1, nbasis):
        rmatrix[i,j] = (rvector[i * nbasis + j] +
                        rvector[j * nbasis + i] * 1j) / sqrt2
        rmatrix[j,i] = (rvector[i * nbasis + j] -
                        rvector[j * nbasis + i] * 1j) / sqrt2
  else:
    amp = np.empty([nbasis], dtype=complex)
    i = -1
    for L in range(Lmax+1):
      for M in range(-L, L+1):
        i += 1
        a = model1.amplitude(L, M, massEtaPi0, abst)
        amp[i] = a[0] + a[1] * 1j
    rmatrix = np.array([[amp[i] * np.conjugate(amp[j]) for i in range(nbasis)]
                                                       for j in range(nbasis)])
  return rmatrix

def get_model1_moments(massEtaPi0, abst):
  """
  Returns the moments vector defined by model 1 at kinematics massEtaPi0,abst. 
  """
  try:
    model1 = ROOT.trial_model1()
  except:
    ROOT.gROOT.ProcessLine(".L src/trial_model.C+O")
    model1 = ROOT.trial_model1()
  Lmax = int(model1.amplitude_Lmax())
  moments = np.empty([(2*Lmax+1)**2], dtype=float)
  i = -1
  for L in range(2*Lmax+1):
    for M in range(-L,L+1):
      i += 1
      moments[i] = model1.real_moment(L, M, massEtaPi0, abst)
  return moments

srcdir = "/".join(__file__.split('/')[:-1]) + "/src"
ROOT.gROOT.ProcessLine(f".L {srcdir}/fill_M_histograms.C+O")
