#!/usr/bin/env python3
#
# analyzeMomentsMatrix.py - general utilities for examining and histogramming
#                           a Monte Carlo acceptance moments matrix.
#
# author: richard.t.jones at uconn.edu
# version: june 10, 2023

MCsample_size = 1e7

import numpy as np
import h5py
import sys
import ROOT
import math
import threading
import scipy.linalg

def usage():
  print("usage: >>> import analyzeMomentsMatrix.py as ana")
  print("  >>> ana.open(\"filename.h5\")")
  print("  >>> w,v = ana.do_svd()")
  print("  >>> h = ana.histogram_M_diagonal()")
  print("  >>> h = ana.histogram_M_diagonal_GJ0()")
  print("  >>> h = ana.histogram_M_diagonal_GJ0_Eta0()")
  print("  >>> h = ana.histogram_M_offdiagonal()")
  print("  >>> h = ana.histogram_eigenvalues(name, title, w)")

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

def histogram_M_diagonal(h5input=0):
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

def histogram_M_diagonal_GJ0(h5input=0):
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

def histogram_M_diagonal_GJ0_Eta0(h5input=0):
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

def list_boxes():
  tmlist = {'all': "etapi0_moments_all.h5"}
  for tlim in ("t(0.0,0.2)", "t(0.2,0.5)", "t(0.5,1.0)"):
    for mlim in ("m(1.0,1.2)", "m(1.2,1.4)", "m(1.4,1.6)", "m(1.6,1.8)", "m(1.8,2.0)"):
      tmlist[tlim+mlim] = f"etapi0_moments_{mlim}_{tlim}.h5"
  return tmlist

def find_all_eigenvalues(half='lower'):
  """
  Argument half can be 'lower', 'upper', or 'sym'.
  """
  boxes = list_boxes()
  h = {}
  w = {}
  v = {}
  M = {}
  Mvar = {}
  for box in boxes:
    h5box = h5py.File(boxes[box], 'r')
    print(box, h5box['generated_subset'][()])
    normfactor = (4 * math.pi)**3 / h5box['generated_subset'][()]
    M[box] = h5box['Moments'][:] * normfactor
    Mvar[box] = h5box['Moments_variance'][:] * normfactor
    if half == 'lower':
       svd = scipy.linalg.eigh(M[box], lower=True)
    elif half == 'upper':
       svd = scipy.linalg.eigh(M[box], lower=False)
    else:
       M[box] = (M[box] + np.transpose(M[box])) / 2
       svd = scipy.linalg.eigh(M[box])
    w[box] = svd[0]
    v[box] = svd[1]
    h[box] = histogram_eigenvalues('heig_' + box, "eigenvalues for bin " + box,  w[box])
  return h,M,Mvar,w,v
