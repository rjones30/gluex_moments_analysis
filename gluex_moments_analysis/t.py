import gluex_moments_analysis.analyzeMomentsMatrix as ana
import numpy as np
import random

k=ana.make_Kmatrix()
b=ana.make_Bmatrix()

def test_Kmatrix(nrandom=1000):
  for i in range(nrandom):
    abst = random.uniform(0,2)
    mX = random.uniform(0.6,2.4)
    h = ana.get_model1_moments(mX, abst)
    r = ana.get_model1_rvector(mX, abst)
    diff = np.linalg.norm(h - k @ r)
    print(f"{i}: {diff}")
    if diff > 1e-12:
      input("ok?")

def test_rmatrix(nrandom=1000):
  for i in range(nrandom):
    abst = random.uniform(0,2)
    mX = random.uniform(0.6,2.4)
    rvector = ana.get_model1_rvector(mX, abst)
    rho = ana.get_model1_rmatrix(mX, abst)
    rrho = ana.get_model1_rmatrix(rvector=rvector)
    diff = np.linalg.norm(rrho - rrho)
    rbr = r @ b @ r
    tracediff = np.linalg.norm(np.trace(rho)**2 - np.trace(rho @ rho))
    print(f"{i}: {diff} {rbr} {tracediff}")
    if diff > 1e-12 or abs(rbr) > 1e-12 or tracediff > 1e-10:
      input("ok?")
