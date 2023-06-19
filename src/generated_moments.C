#include <HDDM/hddm_s.hpp>
#include <iostream>
#include <fstream>
#include <TFile.h>
#include <TTree.h>
#include <TLorentzVector.h>
#include <Math/SpecFuncMathMore.h>

#include "generated_moments.h"

std::string default_outfilename("generated_moments.root");
char *outfilename = (char*)default_outfilename.c_str();

void usage()
{
   std::cout << "usage: generated_moments [options] <infile1.hddm> [...]" << std::endl
             << " where options may include any of the following." << std::endl
             << "  -o <outfile.root> - write output tree to <outfile.root>"
             << ", default is " << outfilename << std::endl;
   exit(1);
  
}

double angular_moment(int L, int M, double theta, double phi)
{
    // These are the real-valued spherical harmonics, Condon and Shortley convention

#define SQRT2 1.4142135623730951

    int Mabs = abs(M);
    double dlm = ROOT::Math::sph_legendre(L, Mabs, theta);
    if (M < 0)
        return SQRT2 * dlm * sin(Mabs * phi) * ((Mabs % 2)? -1 : 1);
    else if (M > 0)
        return SQRT2 * dlm * cos(Mabs * phi) * ((Mabs % 2)? -1 : 1);
    else
        return dlm;
}

int main(int argc, char *argv[])
{
   int iarg(1);
   while (iarg < argc && argv[iarg][0] == '-') {
      if (strstr(argv[iarg], "-o") == 0)
         usage();
      else if (++iarg < argc) {
         outfilename = argv[iarg];
      }
      else {
         usage();
      }
      ++iarg;
   }
   if (iarg >= argc)
      usage();

   TFile ftree(outfilename, "recreate");
   TTree *tree = new TTree("etapi0_moments", "eta pi0 moments");

   tree->Branch("runNo", &runNo, "runNo/L");
   tree->Branch("eventNo", &eventNo, "eventNo/L");
   tree->Branch("weight", &weight, "weight/D");
   tree->Branch("sqrts", &sqrts, "sqrts/D");
   tree->Branch("abst", &abst, "abst/D");
   tree->Branch("massEtaPi0", &massEtaPi0, "massEtaPi0/D");
   tree->Branch("massEta", &massEta, "massEta/D");
   tree->Branch("massPi0", &massPi0, "massPi0/D");
   tree->Branch("phiR", &phiR, "phiR/D");
   tree->Branch("thetaGJ", &thetaGJ, "thetaGJ/D");
   tree->Branch("phiGJ", &phiGJ, "phiGJ/D");
   tree->Branch("thetaEta", &thetaEta, "thetaEta/D");
   tree->Branch("phiEta", &phiEta, "phiEta/D");
   tree->Branch("thetaPi0", &thetaPi0, "thetaPi0/D");
   tree->Branch("phiPi0", &phiPi0, "phiPi0/D");
   tree->Branch("momentsGJ", &momentsGJ, "momentsGJ/I[0,169]");
   tree->Branch("momentsEta", &momentsEta, "momentsEta/I[0,169]");
   tree->Branch("momentsPi0", &momentsPi0, "momentsPi0/I[0,169]");
   tree->Branch("YmomGJ", YmomGJ, "YmomGJ[momentsGJ]/D");
   tree->Branch("YmomEta", YmomEta, "YmomEta[momentsEta]/D");
   tree->Branch("YmomPi0", YmomPi0, "YmomPi0[momentsPi0]/D");

   for (;iarg < argc; ++iarg) {
      std::ifstream infile(argv[iarg]);
      hddm_s::istream inhddm(infile);
      hddm_s::HDDM record;
      while (inhddm >> record) {
         hddm_s::PhysicsEventList pes = record.getPhysicsEvents();
         runNo = pes(0).getRunNo();
         eventNo = pes(0).getEventNo();
         hddm_s::BeamList beams = record.getBeams();
         hddm_s::TargetList targets = record.getTargets();
         hddm_s::ProductList products = record.getProducts();
         hddm_s::Momentum mom_beam = beams(0).getMomenta()(0);
         hddm_s::Momentum mom_targ = targets(0).getMomenta()(0);
         TLorentzVector pbeam(mom_beam.getPx(), mom_beam.getPy(), 
                              mom_beam.getPz(), mom_beam.getE());
         TLorentzVector ptarg(mom_targ.getPx(), mom_targ.getPy(),
                              mom_targ.getPz(), mom_targ.getE());
         TLorentzVector *pp[7] = {&pbeam, &ptarg};
         for (int ip=0; ip < 5; ++ip) {
            hddm_s::Momentum mom_prod = products(ip).getMomenta()(0);
            pp[ip+2] = new TLorentzVector(mom_prod.getPx(), mom_prod.getPy(),
                                          mom_prod.getPz(), mom_prod.getE());
         }
         sqrts = (*pp[0] + *pp[1]).Mag();
         abst = -(*pp[1] - *pp[6]).Mag2();
         massEtaPi0 = (*pp[2] + *pp[3] + *pp[4] + *pp[5]).Mag();
         massEta = (*pp[2] + *pp[3]).Mag();
         massPi0 = (*pp[4] + *pp[5]).Mag();
         phiR = atan2((*pp[6])[1], (*pp[6])[0]);
         for (int i=0; i < 7; ++i)
            pp[i]->RotateZ(M_PI - phiR);
         TLorentzVector betaGJ{*pp[2] + *pp[3] + *pp[4] + *pp[5]};
         betaGJ *= 1 / betaGJ(3);
         for (int i=0; i < 7; ++i)
            pp[i]->Boost(-betaGJ(0), -betaGJ(1), -betaGJ(2));
         double thetaR = atan2((*pp[0])[0], (*pp[0])[2]);
         for (int i=0; i < 7; ++i)
            pp[i]->RotateY(thetaR);
         TLorentzVector pEtaGJ{*pp[2] + *pp[3]};
         thetaGJ = pEtaGJ.Theta();
         phiGJ = pEtaGJ.Phi();
         for (int i=2; i < 4; ++i) {
            pp[i]->RotateZ(-phiGJ);
            pp[i]->RotateY(thetaGJ);
         }
         TLorentzVector beta12{*pp[2] + *pp[3]};
         beta12 *= 1 / beta12(3);
         for (int i=2; i < 4; ++i)
            pp[i]->Boost(-beta12(0), -beta12(1), -beta12(2));
         for (int i=4; i < 6; ++i) {
            pp[i]->RotateZ(-phiGJ);
            pp[i]->RotateY(M_PI - thetaGJ);
         }
         TLorentzVector beta34{*pp[4] + *pp[5]};
         beta34 *= 1 / beta34(3);
         for (int i=4; i < 6; ++i)
            pp[i]->Boost(-beta34(0), -beta34(1), -beta34(2));
         thetaEta = pp[2]->Theta();
         phiEta = pp[2]->Phi();
         thetaPi0 = pp[4]->Theta();
         phiPi0 = pp[4]->Phi();
         tree->Fill();
      }
   }

   tree->Write();
}
