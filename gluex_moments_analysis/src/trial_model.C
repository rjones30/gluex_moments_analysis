//
// trial_model.h - header for class trial_model 
//
// author: richard.t.jones at uconn.edu
// version: july 23, 2023
//
// Defines the partial-wave model for testing analysis of 
// experimental angular distributions in the eta,pi0 system.

#include "trial_model.h"

double trial_model::angular_moment(int L, int M, double theta, double phi)
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

double trial_model::real_moment(int L, int M, double mX, double abst) const
{
   // Computes the model 1 moment corresponding to real spherical harmonic index L,M.
   // Model 1 here refers to a set of functions c[L,M] of arguments (mEtaPi0, abst)
   // that are used to construct a mock data sample for testing a moments extraction
   // procedure. Events moments computed with weigh model1_moment(L,M) emulate a sample
   // whose parent kinematic density function is that produced by the genr8 generator,
   // uniform in CM angles, multiplied by 
   //       __                                               2
   //    |  \                           L                   |
   //    |  /_   c[L,M](mEtaPi0, abst) Y (theta_GJ, phi_GJ) |
   //    |  L,M                         M                   |
   //
   // where the Y{LM} are the ordinary complex spherical harmonics. The evaluation of
   // the functions c[L,M] is delegated to method model1_amplitude(L,M).

   int M0(abs(M));
   std::vector<double> moment{0,0};
   for (int L1=0; L1 <= amplitude_Lmax(); ++L1) {
     for (int M1=-L1; M1 <= L1; ++M1) {
       for (int L2=abs(L-L1); L2 <= amplitude_Lmax(); ++L2) {
         if (L2 <= L+L1 && abs(M0-M1) <= L2) {
            int M2 = M0-M1;
            std::vector<double> amp1 = amplitude(L1, M1, mX, abst);
            std::vector<double> amp2 = amplitude(L2, -M2, mX, abst);
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

double trial_model1::amplitude_Lmax() const
{
   return 3;
}

std::vector<double> trial_model1::amplitude(int L, int M, double mX, double abst) const
{
   // Returns the complex amplitude specified in model 1 for partial wave L,M
   // in the Gottfried Jackson frame, represented as a two-component vector with
   // real part first, imaginary part second.

#define SQR(x) ((x)*(x))

   if (L > amplitude_Lmax()) {
     return std::vector<double>{0,0};
   }
   else if (L == 0 && M == 0) {
     double mag = 11.5 * exp(-0.5 * SQR((mX - 1.0) / 0.15)) * cos(3 * abst);
     return std::vector<double>{mag, mag * (mX - 1.0)};
   }
   else if (L == 1 && M == 1) {
     double mag = 0.5 * exp(-0.5 * SQR((mX - 1.8) / 0.500)) * sin(2 * abst);
     return std::vector<double>{mag, mag * (mX - 1.8)};
   }
   else if (L == 1 && M == 0) {
     double mag = 0.59 * exp(-0.5 * SQR((mX - 1.7) / 0.500));
     return std::vector<double>{mag, mag * (mX - 1.7)};
   }
   else if (L == 1 && M == -1) {
     double mag = -4.5 * exp(-0.5 * SQR((mX - 1.6) / 0.500)) * cos(2 * abst);
     return std::vector<double>{mag, mag * (mX - 1.6)};
   }
   else if (L == 2 && M == 2) {
     double mag = 8.5 * exp(-0.5 * SQR((mX - 1.4) / 0.200)) * sin(2 * abst) / (2*abst + 1e-99);
     return std::vector<double>{mag, mag * (mX - 1.4)};
   }
   else if (L == 2 && M == 1) {
     double mag = 1.5 * exp(-0.5 * SQR((mX - 1.2) / 0.300)) * cos(2 * abst);
     return std::vector<double>{mag, mag * (mX - 1.2)};
   }
   else if (L == 2 && M == 0) {
     double mag = -0.5 * exp(-0.5 * SQR((mX - 1.1) / 0.200)) * sin(1 * abst) / (1*abst + 1e-99);
     return std::vector<double>{mag, mag * (mX - 1.1)};
   }
   else if (L == 2 && M == -1) {
     double mag = -3.5 * exp(-0.5 * SQR((mX - 1.3) / 0.300)) * cos(2 * abst);
     return std::vector<double>{mag, mag * (mX - 1.3)};
   }
   else if (L == 2 && M == -2) {
     double mag = -17.5 * exp(-0.5 * SQR((mX - 1.5) / 0.200)) * sin(2 * abst) / (2*abst + 1e-99);
     return std::vector<double>{mag, mag * (mX - 1.5)};
   }
   else if (L == 3 && M == 3) {
     double mag = 0.85 * exp(-0.5 * SQR((mX - 1.48) / 0.200)) * sin(2 * abst) / (2*abst + 1e-99);
     return std::vector<double>{mag, mag * (mX - 1.48)};
   }
   else if (L == 3 && M == 2) {
     double mag = 0.7 * exp(-0.5 * SQR((mX - 1.23) / 0.300)) * cos(2 * abst);
     return std::vector<double>{mag, mag * (mX - 1.23)};
   }
   else if (L == 3 && M == 1) {
     double mag = -0.15 * exp(-0.5 * SQR((mX - 1.18) / 0.200)) * sin(1 * abst) / (1*abst + 1e-99);
     return std::vector<double>{mag, mag * (mX - 1.18)};
   }
   else if (L == 3 && M == 0) {
     double mag = -9.5 * exp(-0.5 * SQR((mX - 1.32) / 0.300)) * cos(2 * abst);
     return std::vector<double>{mag, mag * (mX - 1.32)};
   }
   else if (L == 3 && M == -1) {
     double mag = -0.95 * exp(-0.5 * SQR((mX - 1.56) / 0.200)) * sin(2 * abst) / (2*abst + 1e-99);
     return std::vector<double>{mag, mag * (mX - 1.56)};
   }
   else if (L == 3 && M == -2) {
     double mag = 0.65 * exp(-0.5 * SQR((mX - 1.26) / 0.200)) * sin(2 * abst) / (2*abst + 1e-99);
     return std::vector<double>{mag, mag * (mX - 1.26)};
   }
   else if (L == 3 && M == -3) {
     double mag = -4.95 * exp(-0.5 * SQR((mX - 1.66) / 0.200)) * sin(2 * abst) / (2*abst + 1e-99);
     return std::vector<double>{mag, mag * (mX - 1.66)};
   }
   else {
     return std::vector<double>{0,0};
   }
}

double trial_model::angular_density(double theta, double phi, 
                                    double mX, double abst,
                                    int from_amplitudes) const
{
   // Returns the model density at Gottfried-Jackson angles theta,phi
   // computed either from the complex amplitude sum (from_amplitudes==1)
   // or from the sum over angular moments (from_amplitudes==0).

   if (from_amplitudes) {
     double areal(0), aimag(0);
     for (int L=0; L <= amplitude_Lmax(); ++L) {
       for (int M=-L; M <= L; ++M) {
         std::vector<double> clm = amplitude(L, M, mX, abst);
         double dlm = ROOT::Math::sph_legendre(L, abs(M), theta);
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
     for (int L=0; L <= 2*amplitude_Lmax(); ++L) {
       for (int M=-L; M <= L; ++M) {
         pdf += angular_moment(L, M, theta, phi) * real_moment(L, M, mX, abst);
       }
     }
     return pdf;
   }
}
