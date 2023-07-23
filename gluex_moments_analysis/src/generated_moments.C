#include <HDDM/hddm_s.hpp>
#include <iostream>
#include <fstream>
#include <TFile.h>
#include <TTree.h>
#include <TLorentzVector.h>
#include <Math/SpecFuncMathMore.h>
#include <climits>

#include "DSelector_etapi0_moments.h"
#include "generated_moments.h"

std::string default_outfilename("generated_moments.root");
char *outfilename = (char*)default_outfilename.c_str();

DSelector_etapi0_moments models;

void usage()
{
   std::cout << "usage: generated_moments [options] <infile1.hddm> [...]" << std::endl
             << " where options may include any of the following." << std::endl
             << "  -n <noutput> - stop after writing <noutput> events [all]" << std::endl
             << "  -s <nskip> - skip <nskip> events at start of input [0]" << std::endl
             << "  -o <outfile.root> - write output tree to <outfile.root>"
             << ", default is " << outfilename << std::endl;
   exit(1);
  
}

int main(int argc, char *argv[])
{
   long int nskip(0);
   long int maxnout(LLONG_MAX);

   int iarg(1);
   while (iarg < argc && argv[iarg][0] == '-') {
      if (strstr(argv[iarg], "-o") == argv[iarg]) {
         outfilename = argv[++iarg];
      }
      else if (strstr(argv[iarg], "-n") == argv[iarg]) {
         maxnout = std::atoi(argv[++iarg]);
      }
      else if (strstr(argv[iarg], "-s") == argv[iarg]) {
         nskip = std::atoi(argv[++iarg]);
      }
      else {
         usage();
      }
      ++iarg;
   }
   if (iarg >= argc)
      usage();

   std::cout << "writing at most " << maxnout << " events to " << outfilename
             << " starting with input event " << nskip << std::endl;
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
   tree->Branch("model1moment", model1moment, "model1moment[momentsGJ]/D");

   int nout(0);
   for (;iarg < argc; ++iarg) {
      std::ifstream infile(argv[iarg]);
      hddm_s::istream inhddm(infile);
      if (nskip > 0) {
         inhddm.skip(nskip);
	 nskip = 0;
      }
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

         models.model1_use_generated(true);
         models.massEtaPi0_ = massEtaPi0;
         models.abst_ = abst;
         momentsGJ = 0;
         for (int L=0; L < 13; L += 1) {
            for (int M=-L; M <= L; M += 1) {
               YmomGJ[momentsGJ] = models.angular_moment(L, M, thetaGJ, phiGJ);
               model1moment[momentsGJ] = models.model1_moment(L, M);
               ++momentsGJ;
            }
         }
         momentsEta = 0;
         for (int L=0; L < 13; L += 2) {
            for (int M=-L; M <= L; M += 1) {
               YmomEta[momentsEta] = models.angular_moment(L, M, thetaEta, phiEta);
               ++momentsEta;
            }
         }
         momentsPi0 = 0;
         for (int L=0; L < 13; L += 2) {
            for (int M=-L; M <= L; M += 1) {
               YmomPi0[momentsPi0] = models.angular_moment(L, M, thetaPi0, phiPi0);
               ++momentsPi0;
            }
         }
         tree->Fill();
         if (++nout == maxnout)
            break;
      }
   }

   tree->Write();
}
