//
// trial_model.h - header for class trial_model 
//
// author: richard.t.jones at uconn.edu
// version: july 23, 2023
//
// Defines the partial-wave model for testing analysis of 
// experimental angular distributions in the eta,pi0 system.

#ifndef TRIAL_MODEL_H
#define TRIAL_MODEL_H

#include <vector>
#include <Math/SpecFuncMathMore.h>

class trial_model {
 public:
   virtual ~trial_model() {}

   static double angular_moment(int L, int M, double theta, double phi);

   virtual double real_moment(int L, int M, double mX, double abst) const;
   virtual double angular_density(double theta, double phi, double mX, double abst,
                                  int from_amplitudes=0) const;

   virtual double amplitude_Lmax() const = 0;
   virtual std::vector<double> amplitude(int L, int M, double mX, double abst) const = 0;
};

class trial_model1 : public trial_model {
 public:
   virtual ~trial_model1() {}

   virtual double amplitude_Lmax() const;
   virtual std::vector<double> amplitude(int L, int M, double mX, double abst) const;
};

#endif
