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

double model1_Lmax()
{
   return 3;
}

std::vector<double> model1_amplitude(int L, int M)
{
   // Returns the complex amplitude specified in model 1 for partial wave L,M
   // in the Gottfried Jackson frame, represented as a two-component vector with
   // real part first, imaginary part second.

#define SQR(x) ((x)*(x))

   double mX(massEtaPi0);
   double t(abst);

   if (L > model1_Lmax()) {
      return std::vector<double>{0,0};
   }
   else if (L == 0 && M == 0) {
      double mag = 11.5 * exp(-0.5 * SQR((mX - 1.0) / 0.15)) * cos(3 * t);
      return std::vector<double>{mag, mag * (mX - 1.0)};
   }
   else if (L == 1 && M == 1) {
      double mag = 0.5 * exp(-0.5 * SQR((mX - 1.8) / 0.500)) * sin(2 * t);
      return std::vector<double>{mag, mag * (mX - 1.8)};
   }
   else if (L == 1 && M == 0) {
      double mag = 0.59 * exp(-0.5 * SQR((mX - 1.7) / 0.500));
      return std::vector<double>{mag, mag * (mX - 1.7)};
   }
   else if (L == 1 && M == -1) {
      double mag = -4.5 * exp(-0.5 * SQR((mX - 1.6) / 0.500)) * cos(2 * t);
      return std::vector<double>{mag, mag * (mX - 1.6)};
   }
   else if (L == 2 && M == 2) {
      double mag = 8.5 * exp(-0.5 * SQR((mX - 1.4) / 0.200)) * sin(2 * t) / (2*t + 1e-99);
      return std::vector<double>{mag, mag * (mX - 1.4)};
   }
   else if (L == 2 && M == 1) {
      double mag = 1.5 * exp(-0.5 * SQR((mX - 1.2) / 0.300)) * cos(2 * t);
      return std::vector<double>{mag, mag * (mX - 1.2)};
   }
   else if (L == 2 && M == 0) {
      double mag = -0.5 * exp(-0.5 * SQR((mX - 1.1) / 0.200)) * sin(1 * t) / (1*t + 1e-99);
      return std::vector<double>{mag, mag * (mX - 1.1)};
   }
   else if (L == 2 && M == -1) {
      double mag = -3.5 * exp(-0.5 * SQR((mX - 1.3) / 0.300)) * cos(2 * t);
      return std::vector<double>{mag, mag * (mX - 1.3)};
   }
   else if (L == 2 && M == -2) {
      double mag = -17.5 * exp(-0.5 * SQR((mX - 1.5) / 0.200)) * sin(2 * t) / (2*t + 1e-99);
      return std::vector<double>{mag, mag * (mX - 1.5)};
   }
   else if (L == 3 && M == 3) {
      double mag = 0.85 * exp(-0.5 * SQR((mX - 1.48) / 0.200)) * sin(2 * t) / (2*t + 1e-99);
      return std::vector<double>{mag, mag * (mX - 1.48)};
   }
   else if (L == 3 && M == 2) {
      double mag = 0.7 * exp(-0.5 * SQR((mX - 1.23) / 0.300)) * cos(2 * t);
      return std::vector<double>{mag, mag * (mX - 1.23)};
   }
   else if (L == 3 && M == 1) {
      double mag = -0.15 * exp(-0.5 * SQR((mX - 1.18) / 0.200)) * sin(1 * t) / (1*t + 1e-99);
      return std::vector<double>{mag, mag * (mX - 1.18)};
   }
   else if (L == 3 && M == 0) {
      double mag = -9.5 * exp(-0.5 * SQR((mX - 1.32) / 0.300)) * cos(2 * t);
      return std::vector<double>{mag, mag * (mX - 1.32)};
   }
   else if (L == 3 && M == -1) {
      double mag = -0.95 * exp(-0.5 * SQR((mX - 1.56) / 0.200)) * sin(2 * t) / (2*t + 1e-99);
      return std::vector<double>{mag, mag * (mX - 1.56)};
   }
   else if (L == 3 && M == -2) {
      double mag = 0.65 * exp(-0.5 * SQR((mX - 1.26) / 0.200)) * sin(2 * t) / (2*t + 1e-99);
      return std::vector<double>{mag, mag * (mX - 1.26)};
   }
   else if (L == 3 && M == -3) {
      double mag = -4.95 * exp(-0.5 * SQR((mX - 1.66) / 0.200)) * sin(2 * t) / (2*t + 1e-99);
      return std::vector<double>{mag, mag * (mX - 1.66)};
   }
   else {
      return std::vector<double>{0,0};
   }
}

double model1_moment(int L, int M)
{
   // Computes the model 1 moment corresponding to real spherical harmonic index L,M.
   // Model 1 here refers to a set of functions c[L,M] of arguments (mEtaPi0, abst)
   // that are used to construct a mock data sample for testing a moments extraction
   // procedure. Events moments computed with weigh model1_moment(L,M) emulate a sample
   // whose parent kinematic density function is that produced by the genr8 generator,
   // uniform in CM angles, multiplied by 
   //        __                                                2
   //     |  \                           L                   |
   //     |  /_   c[L,M](mEtaPi0, abst) Y (theta_GJ, phi_GJ) |
   //     |  L,M                         M                   |
   //
   // where the Y{LM} are the ordinary complex spherical harmonics. The evaluation of
   // the functions c[L,M] is delegated to method model1_amplitude(L,M).

   int M0(abs(M));
   std::vector<double> moment{0,0};
   for (int L1=0; L1 <= model1_Lmax(); ++L1) {
      for (int M1=-L1; M1 <= L1; ++M1) {
         for (int L2=abs(L-L1); L2 <= model1_Lmax(); ++L2) {
            if (L2 <= L+L1 && abs(M0-M1) <= L2) {
               int M2 = M0-M1;
               std::vector<double> amp1 = model1_amplitude(L1, M1);
               std::vector<double> amp2 = model1_amplitude(L2, -M2);
               double a = ROOT::Math::wigner_3j(2*L1, 2*L2, 2*L, 2*M1, 2*M2, -2*M0)
                        * ROOT::Math::wigner_3j(2*L1, 2*L2, 2*L, 0, 0, 0)
                        * ((M2 % 2)? -1 : 1)
                        * sqrt(2*L1+1) * sqrt(2*L2+1) * sqrt(2*L+1)
                        / sqrt(4*M_PI);
               moment[0] += a * (amp1[0] * amp2[0] + amp1[1] * amp2[1]);
               moment[1] += a * (amp1[1] * amp2[0] - amp1[0] * amp2[1]);
            }
         }
      }
   }
   if (M == 0)
      return moment[0];
   else if (M > 0)
      return moment[0] * SQRT2;
   else
      return -moment[1] * SQRT2;
}

double model1_density(int source)
{
   // Returns the model 1 density at Gottfried-Jackson angles thetaGJ,phiGJ
   // computed either from the complex amplitude sum (source==0) or from
   // a sum over computed moments (source==1).

   double theta(thetaGJ);
   double phi(phiGJ);

   if (source == 0) {
      double areal(0), aimag(0);
      for (int L=0; L <= model1_Lmax(); ++L) {
         for (int M=-L; M <= L; ++M) {
            std::vector<double> clm = model1_amplitude(L, M);
            double dlm = ROOT::Math::sph_legendre(L, abs(M), thetaGJ);
            if (M < 0) {
               dlm *= ((abs(M) % 2)? -1 : 1);
            }
            areal += dlm * (clm[0] * cos(M * phi) - clm[1] * sin(M * phi));
            aimag += dlm * (clm[0] * sin(M * phi) + clm[1] * cos(M * phi));
         }
      }
      return areal*areal + aimag*aimag;
   }
   else {
      double pdf = 0;
      for (int L=0; L <= 2*model1_Lmax(); ++L) {
         for (int M=-L; M <= L; ++M) {
            pdf += angular_moment(L, M, theta, phi) * model1_moment(L, M);
         }
      }
      return pdf;
   }
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
   tree->Branch("model1moment", model1moment, "model1moment[momentsGJ]/D");

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

         momentsGJ = 0;
         for (int L=0; L < 13; L += 1) {
            for (int M=-L; M <= L; M += 1) {
               YmomGJ[momentsGJ] = angular_moment(L, M, thetaGJ, phiGJ);
               model1moment[momentsGJ] = model1_moment(L, M);
               ++momentsGJ;
            }
         }
         momentsEta = 0;
         for (int L=0; L < 13; L += 2) {
            for (int M=-L; M <= L; M += 1) {
               YmomEta[momentsEta] = angular_moment(L, M, thetaEta, phiEta);
               ++momentsEta;
            }
         }
         momentsPi0 = 0;
         for (int L=0; L < 13; L += 2) {
            for (int M=-L; M <= L; M += 1) {
               YmomPi0[momentsPi0] = angular_moment(L, M, thetaPi0, phiPi0);
               ++momentsPi0;
            }
         }
         tree->Fill();
      }
   }

   tree->Write();
}
