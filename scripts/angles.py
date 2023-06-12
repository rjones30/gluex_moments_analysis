#!/usr/bin/env python3

import hddm_s
import ROOT
import sys
import numpy as np

hm12 = ROOT.TH1D("hm12", "mass of gammas 1,2", 1000, 0, 2)
hm12.GetXaxis().SetTitle("mass of gammas 1,2 (GeV/c^2)")
hm12.GetYaxis().SetTitle("events")
hm34 = ROOT.TH1D("hm34", "mass of gammas 3,4", 1000, 0, 2)
hm34.GetXaxis().SetTitle("mass of gammas 3,4 (GeV/c^2)")
hm34.GetYaxis().SetTitle("events")
hcost12 = ROOT.TH1D("hcost12", "costhetaCM of gammas 1,2", 100, -1, 1)
hcost12.GetXaxis().SetTitle("cos(thetaCM) of gammas 1")
hcost12.GetYaxis().SetTitle("events")
hcost34 = ROOT.TH1D("hcost34", "costhetaCM of gammas 3,4", 100, -1, 1)
hcost34.GetXaxis().SetTitle("cos(thetaCM) of gammas 3")
hcost34.GetYaxis().SetTitle("events")
hcostGJ = ROOT.TH1D("hcostGJ", "costhetaGJ of eta", 100, -1, 1)
hcostGJ.GetXaxis().SetTitle("cos(thetaGJ)")
hcostGJ.GetYaxis().SetTitle("events")
hphi12 = ROOT.TH1D("hphi12", "phiCM of gammas 1,2", 100, -180, 180)
hphi12.GetXaxis().SetTitle("phiCM of gammas 1 (deg)")
hphi12.GetYaxis().SetTitle("events")
hphi34 = ROOT.TH1D("hphi34", "phiCM of gammas 3,4", 100, -180, 180)
hphi34.GetXaxis().SetTitle("phiCM of gammas 3 (deg)")
hphi34.GetYaxis().SetTitle("events")
hphiGJ = ROOT.TH1D("hphiGJ", "phiGJ of eta", 100, -180, 180)
hphiGJ.GetXaxis().SetTitle("phiGJ (deg)")
hphiGJ.GetYaxis().SetTitle("events")
hmandt = ROOT.TH1D("hmandt", "Mandelstam -t", 100, 0, 3)
hmandt.GetXaxis().SetTitle("|t| (GeV^2)")
hmandt.GetYaxis().SetTitle("events")
hphiR = ROOT.TH1D("hphiR", "phi of recoil proton", 100, -180, 180)
hphiR.GetXaxis().SetTitle("phiRecoil (deg)")
hphiR.GetYaxis().SetTitle("events")

for fin in sys.argv[1:]:
   p = [0]*7
   for rec in hddm_s.istream(fin):
      moms = rec.getMomenta()
      for i in range(7):
         p[i] = ROOT.TLorentzVector(moms[i].px, moms[i].py, moms[i].pz, moms[i].E)
      mandt = (p[6] - p[1]).Mag2()
      phiR = np.arctan2(p[6][1], p[6][0])
      for i in range(7):
         p[i].RotateZ(np.pi - phiR)
      pfwd = p[2] + p[3] + p[4] + p[5]
      for i in range(7):
         p[i].Boost(-pfwd[0] / pfwd[3], -pfwd[1] / pfwd[3], -pfwd[2] / pfwd[3])
      thetaR = np.arctan2(p[0][0], p[0][2])
      for i in range(7):
         p[i].RotateY(thetaR)
      p12 = p[2] + p[3]
      thetaGJ = p12.Theta()
      phiGJ = p12.Phi()
      for i in range(2,4):
         p[i].RotateZ(-phiGJ)
         p[i].RotateY(thetaGJ)
      p12 = p[2] + p[3]
      for i in range(2,4):
         p[i].Boost(-p12[0] / p12[3], -p12[1] / p12[3], -p12[2] / p12[3])
      for i in range(4,6):
         p[i].RotateZ(-phiGJ)
         p[i].RotateY(np.pi - thetaGJ)
      p34 = p[4] + p[5]
      for i in range(4,6):
         p[i].Boost(-p34[0] / p34[3], -p34[1] / p34[3], -p34[2] / p34[3])
      theta12 = p[2].Theta()
      phi12 = p[2].Phi()
      theta34 = p[4].Theta()
      phi34 = p[4].Phi()
      hm12.Fill(p12.Mag())
      hm34.Fill(p34.Mag())
      hcost12.Fill(np.cos(theta12))
      hcost34.Fill(np.cos(theta34))
      hcostGJ.Fill(np.cos(thetaGJ))
      hphi12.Fill(phi12 * 180/np.pi)
      hphi34.Fill(phi34 * 180/np.pi)
      hphiGJ.Fill(phiGJ * 180/np.pi)
      hmandt.Fill(-mandt)
      hphiR.Fill(phiR * 180/np.pi)

fout = ROOT.TFile("angles.root", "recreate")
hm12.Write()
hm34.Write()
hcost12.Write()
hcost34.Write()
hcostGJ.Write()
hphi12.Write()
hphi34.Write()
hphiGJ.Write()
hmandt.Write()
hphiR.Write()
