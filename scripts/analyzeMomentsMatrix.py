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

import buildMomentsMatrix as bmm

def usage():
  print("usage: >>> import analyzeMomentsMatrix.py as ana")
  print("  >>> ana.open(\"filename.h5\")")
  print("  >>> w,v = ana.do_svd()")
  print("  >>> h = ana.histogram_M_diagonal()")
  print("  >>> h = ana.histogram_M_diagonal_GJ0()")
  print("  >>> h = ana.histogram_M_diagonal_GJ0_Eta0()")
  print("  >>> h = ana.histogram_M_offdiagonal()")
  print("  >>> h = ana.histogram_eigenvalues(name, title, w)")
  print("  >>> h = ana.histogram_moments(name, title, moments, errors)")

if len(sys.argv) < 2 or sys.argv[1][0] == '-':
  usage()

h5inputs = {}

def open(h5input):
  global f5
  global M
  global Mvar
  f5 = h5py.File(h5input, 'r')
  M = f5['Moments']
  Mvar = f5['Moments_variance']
  h5inputs[0] = f5
  try:
    accepted_subset = f5['accepted_subset'][()]
  except:
    accepted_subset = 0
  return accepted_subset

def do_svd(h5input=0, lower=True):
  if not h5input in h5inputs:
    try:
      h5inputs[0] = 1
    except:
      print("unable to open input moments file", h5input)
      return 0
  return scipy.linalg.eigh(M, lower=lower)

def histogram_M_diagonal(h5input=0, MCsample_size=1e7):
  if not h5input in h5inputs:
    try:
      h5inputs[open(h5input)] = 1
    except:
      print("unable to open input moments file", h5input)
      return 0
  mdim = M.shape[0]
  hMdiag = ROOT.TH1D("hMdiag", "M diagonal elements", mdim, 0, mdim)
  hMdiag.GetXaxis().SetTitle("moment index")
  hMdiag.GetYaxis().SetTitle("M[i,i]")
  hMdiag.SetMinimum(0)
  normfact = (4 * math.pi)**3 / MCsample_size
  smallest = [1e9,0]
  largest = [0,0]
  smallerr = [1e9,0]
  largerr = [0,0]
  for m in range(mdim):
    hMdiag.SetBinContent(m+1, M[m,m] * normfact)
    hMdiag.SetBinError(m+1, Mvar[m,m]**0.5 * normfact)
    if M[m,m] < smallest[0]:
      smallest = [M[m,m], Mvar[m,m]]
    elif M[m,m] > largest[0]:
      largest = [M[m,m], Mvar[m,m]]
    if Mvar[m,m]**0.5 / M[m,m] > largerr[0]:
      largerr = [Mvar[m,m]**0.5 / M[m,m], 0]
    elif Mvar[m,m]**0.5 / M[m,m] < smallerr[0]:
      smallerr = [Mvar[m,m]**0.5 / M[m,m], 0]
  print("smallest diagonal element is ", end='')
  print(smallest[0] * normfact, "+/-", smallest[1]**0.5 * normfact, end='')
  print(", largest is ", end='')
  print(largest[0] * normfact, "+/-", largest[1]**0.5 * normfact)
  print("smallest fractional error is ", end='')
  print(smallerr[0], "+/-", smallerr[1], end='')
  print(", largest is ", end='')
  print(largerr[0], "+/-", largerr[1])
  return hMdiag

def histogram_M_diagonal_GJ0(h5input=0, MCsample_size=1e7):
  if not h5input in h5inputs:
    try:
      h5inputs[open(h5input)] = 1
    except:
      print("unable to open input moments file", h5input)
      return 0
  mdim = M.shape[0]
  hMdiag = ROOT.TH1D("hMdiagGJ0", "M diagonal elements with L(GJ)=0", 768, 0, 768)
  hMdiag.GetXaxis().SetTitle("moment index")
  hMdiag.GetYaxis().SetTitle("M[i,i]")
  hMdiag.SetMinimum(0)
  normfact = (4 * math.pi)**3 / MCsample_size
  for m in range(mdim):
    if (m % 49) == 0:
      hMdiag.SetBinContent(int(m/49)+1, M[m,m] * normfact)
      hMdiag.SetBinError(int(m/49)+1, Mvar[m,m]**0.5 * normfact)
  return hMdiag

def histogram_M_diagonal_GJ0_Eta0(h5input=0, MCsample_size=1e7):
  if not h5input in h5inputs:
    try:
      h5inputs[open(h5input)] = 1
    except:
      print("unable to open input moments file", h5input)
      return 0
  mdim = M.shape[0]
  hMdiag = ROOT.TH1D("hMdiagGJ0Eta0", "M diagonal elements with L(GJ)=0, L(Eta)=0", 28, 0, 28)
  hMdiag.GetXaxis().SetTitle("moment index")
  hMdiag.GetYaxis().SetTitle("M[i,i]")
  hMdiag.SetMinimum(0)
  normfact = (4 * math.pi)**3 / MCsample_size
  for m in range(mdim):
    if (m % 1372) == 0:
      hMdiag.SetBinContent(int(m/1372)+1, M[m,m] * normfact)
      hMdiag.SetBinError(int(m/1372)+1, Mvar[m,m]**0.5 * normfact)
  return hMdiag

def histogram_M_offdiagonal(h5input=0):
  if not h5input in h5inputs:
    try:
      h5inputs[open(h5input)] = 1
    except:
      print("unable to open input moments file", h5input)
      return 0
  Marray = np.asarray(M)
  Mvararray = np.asarray(Mvar)
  mdim = M.shape[0]
  nthreads = 10
  hMabove = []
  hMbelow = []
  hEabove = []
  hEbelow = []
  for i in range(nthreads):
    hMabovediag = ROOT.TH1D(f"hMabovediag{i}", "M above-diagonal elements normalized to diagonal", 200,-0.1, 0.1)
    hMbelowdiag = ROOT.TH1D(f"hMbelowdiag{i}", "M below-diagonal elements normalized to diagonal", 200,-0.1, 0.1)
    hEabovediag = ROOT.TH1D(f"hEabovediag{i}", "M above-diagonal error normalized to diagonal", 200,-0.1, 0.1)
    hEbelowdiag = ROOT.TH1D(f"hEbelowdiag{i}", "M below-diagonal error normalized to diagonal", 200,-0.1, 0.1)
    hMabovediag.GetXaxis().SetTitle("normalized M matrix element")
    hMabovediag.GetYaxis().SetTitle("moments")
    hMbelowdiag.GetXaxis().SetTitle("normalized M matrix element")
    hMbelowdiag.GetYaxis().SetTitle("moments")
    hEabovediag.GetXaxis().SetTitle("normalized M matrix element error")
    hEabovediag.GetYaxis().SetTitle("moments")
    hEbelowdiag.GetXaxis().SetTitle("normalized M matrix element error")
    hEbelowdiag.GetYaxis().SetTitle("moments")
    hMabove.append(hMabovediag)
    hMbelow.append(hMbelowdiag)
    hEabove.append(hEabovediag)
    hEbelow.append(hEbelowdiag)
  ROOT.gROOT.ProcessLine(".L fill_M_histograms.C+O")
  ROOT.fill_M_histograms.__release_gil__ = False
  ithread = 0
  threads = []
  for m in range(mdim):
    t = threading.Thread(target=ROOT.fill_M_histograms, args=(Marray, Mvararray, m, mdim, hMabove[ithread], hMbelow[ithread], hEabove[ithread], hEbelow[ithread]))
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
  hEabovediag = hEabove[0]
  hEbelowdiag = hEbelow[0]
  for i in range(1,nthreads):
    hMabovediag.Add(hMabove[i])
    hMbelowdiag.Add(hMbelow[i])
    hEabovediag.Add(hEabove[i])
    hEbelowdiag.Add(hEbelow[i])
  return (hMabovediag,hMbelowdiag,hEabovediag,hEbelowdiag)

def histogram_offdiagonal(name, M, Mvar):
  hMabovediag = ROOT.TH1D(f"hMabovediag{name}", "M above-diagonal elements normalized to diagonal", 10000,-0.1, 0.1)
  hMbelowdiag = ROOT.TH1D(f"hMbelowdiag{name}", "M below-diagonal elements normalized to diagonal", 10000,-0.1, 0.1)
  hEabovediag = ROOT.TH1D(f"hEabovediag{name}", "M above-diagonal errors normalized to diagonal", 10000,-0.1, 0.1)
  hEbelowdiag = ROOT.TH1D(f"hEbelowdiag{name}", "M below-diagonal errors normalized to diagonal", 10000,-0.1, 0.1)
  hMabovediag.GetXaxis().SetTitle("normalized M matrix element")
  hMabovediag.GetYaxis().SetTitle("moments")
  hMbelowdiag.GetXaxis().SetTitle("normalized M matrix element")
  hMbelowdiag.GetYaxis().SetTitle("moments")
  hEabovediag.GetXaxis().SetTitle("normalized M matrix errors")
  hEabovediag.GetYaxis().SetTitle("moments")
  hEbelowdiag.GetXaxis().SetTitle("normalized M matrix errors")
  hEbelowdiag.GetYaxis().SetTitle("moments")
  mdim = M.shape[0]
  for i in range(mdim):
    for j in range(i):
      hMbelowdiag.Fill(M[i,j] / (M[i,i] * M[j,j])**0.5)
      hEbelowdiag.Fill((Mvar[i,j] / (M[i,i] * M[j,j]))**0.5)
    for j in range(i+1,mdim):
      hMabovediag.Fill(M[i,j] / (M[i,i] * M[j,j])**0.5)
      hEabovediag.Fill((Mvar[i,j] / (M[i,i] * M[j,j]))**0.5)
  return (hMabovediag, hMbelowdiag, hEabovediag, hEbelowdiag)

def histogram_eigenvalues(name, title, w, normfactor=1):
  ndim = len(w)
  h = ROOT.TH1D(name, title, ndim, 0, ndim)
  for i in range(ndim):
    h.SetBinContent(ndim - i, w[i] * normfactor)
  h.GetXaxis().SetTitle("eigenvalue index")
  h.GetYaxis().SetTitle("acceptance eigenvalue")
  return h

def histogram_moments(name, title, moments, errors):
  ndim = len(moments)
  h = ROOT.TH1D(name, title, ndim, 0, ndim)
  for i in range(ndim):
    h.SetBinContent(i+1, moments[i])
    h.SetBinError(i+1, errors[i])
  h.GetXaxis().SetTitle("moment index")
  h.GetYaxis().SetTitle("sample moment")
  return h

def list_boxes():
  tmlist = {'all': "etapi0_moments_all.h5"}
  for tlim in ("t(0.0,0.2)", "t(0.2,0.5)", "t(0.5,1.0)"):
    for mlim in ("m(1.0,1.2)", "m(1.2,1.4)", "m(1.4,1.6)", "m(1.6,1.8)", "m(1.8,2.0)"):
      tmlist[tlim+mlim] = f"etapi0_moments_{mlim}_{tlim}.h5"
  return tmlist

def hinhout(massEtaPi0_limits=(0,99), abst_limits=(0,99), model=1, sample_subset=range(10), maxweight=0):
  try:
    f5 = h5py.File("Msaved.h5")
    M = f5['Moments']
    Mvar = f5['Moments_variance']
    Nacc = f5['accepted_subset'][()]
    Ngen = f5['generated_subset'][()]
  except:
    try:
      f5.close()
    except:
      pass
    print("hinhout regenerating Msaved.h5")
    acc_union = []
    gen_union = []
    for i in sample_subset:
      bmm.open(f"../etapi0_moments_x10_{i}.root:etapi0_moments")
      acc_events = bmm.select_events(massEtaPi0_limits=massEtaPi0_limits, abst_limits=abst_limits)
      M_,Mvar_ = bmm.buildMomentsMatrix_sequential(acc_events, mPi0=1, mEta=1)
      acc_union += acc_events
      try:
        M += M_
        Mvar += Mvar_
      except:
        M = M_
        Mvar = Mvar_
      bmm.open(f"../generated_moments_x10_{i}.root:etapi0_moments")
      gen_events = bmm.select_events(massEtaPi0_limits=massEtaPi0_limits, abst_limits=abst_limits)
      gen_union += gen_events
      """ disable this, unless you need to check the orthonormality of the moments
      M_,Mvar_ = bmm.buildMomentsMatrix_sequential(gen_events, mPi0=1, mEta=1)
      try:
        Mgen += M_
        Mgenvar += Mvar_
      except:
        Mgen = M_
        Mgenvar = Mvar_
    bmm.save_output(Mgen, Mgenvar, gen_union, gen_union, "Mperfect.h5")
      """
    bmm.save_output(M, Mvar, acc_union, gen_union, "Msaved.h5")
    Nacc = len(acc_union)
    Ngen = len(gen_union)

  M *= (4 * math.pi)**3 / Ngen
  Mvar *= ((4 * math.pi)**3 / Ngen)**2

  try:
    f5 = h5py.File("Msample.h5")
    samplemom = f5['sample_moment']
    model1mom = f5['model1_moment']
    samplevar = f5['sample_variance']
    model1var = f5['model1_variance']
  except:
    try:
      f5.close()
    except:
      pass
    print("hinhout regenerating Msample.h5")
    acc_union = []
    gen_union = []
    for i in sample_subset:
      bmm.open(f"../etapi0_moments_x10_{i}.root:etapi0_moments")
      acc_events = bmm.select_events(massEtaPi0_limits=massEtaPi0_limits, abst_limits=abst_limits, model=1, maxweight=maxweight)
      samplemom_, samplevar_ = bmm.compute_moments(acc_events, mPi0=1, mEta=1, use_generated=1)
      try:
        samplemom += samplemom_
        samplevar += samplevar_
      except:
        samplemom = samplemom_
        samplevar = samplevar_
      acc_union += acc_events
      bmm.open(f"../generated_moments_x10_{i}.root:etapi0_moments")
      gen_events = bmm.select_events(massEtaPi0_limits=massEtaPi0_limits, abst_limits=abst_limits, model=1, maxweight=maxweight)
      model1mom_, model1var_ = bmm.compute_moments(gen_events, mPi0=1, mEta=1)
      try:
        model1mom += model1mom_
        model1var += model1var_
      except:
        model1mom = model1mom_
        model1var = model1var_
      gen_union += gen_events
    f5 = h5py.File("Msample.h5", "w")
    f5.create_dataset("sample_moment", data=samplemom)
    f5.create_dataset("model1_moment", data=model1mom)
    f5.create_dataset("sample_variance", data=samplevar)
    f5.create_dataset("model1_variance", data=model1var)
    f5.close()

  Minv = 2 * np.linalg.inv(M + np.transpose(M))
  correctmom = Minv @ samplemom
  correctvar = np.diag(Minv @ np.diag(samplevar) @ np.transpose(Minv))
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
    hsamp.SetBinError(i+1, samplevar[i]**0.5)
    hcorr.SetBinContent(i+1, correctmom[i])
    hcorr.SetBinError(i+1, correctvar[i]**0.5)
    hmodel.SetBinContent(i+1, model1mom[i])
    hmodel.SetBinError(i+1, model1var[i]**0.5)
    hdiff.SetBinContent(i+1, correctmom[i])
    hdiff.SetBinError(i+1, correctvar[i]**0.5)
  hdiff.Add(hmodel, -1)
  return hsamp,hcorr,hmodel,hdiff
