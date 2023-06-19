#ifndef DSelector_etapi0_moments_h
#define DSelector_etapi0_moments_h

#include <iostream>

#include "DSelector/DSelector.h"
#include "DSelector/DHistogramActions.h"
#include "DSelector/DCutActions.h"

#include "TH1I.h"
#include "TH2I.h"

class DSelector_etapi0_moments : public DSelector
{
	public:

		DSelector_etapi0_moments(TTree* locTree = NULL) : DSelector(locTree){}
		virtual ~DSelector_etapi0_moments(){}

		void Init(TTree *tree);
		Bool_t Process(Long64_t entry);

	private:

		void Get_ComboWrappers(void);
		void Finalize(void);

		// BEAM POLARIZATION INFORMATION
		UInt_t dPreviousRunNumber;
		bool dIsPolarizedFlag; //else is AMO
		bool dIsPARAFlag; //else is PERP or AMO

		bool dIsMC;

		// ANALYZE CUT ACTIONS
		// // Automatically makes mass histograms where one cut is missing
		DHistogramAction_AnalyzeCutActions* dAnalyzeCutActions;

		//CREATE REACTION-SPECIFIC PARTICLE ARRAYS

		//Step 0
		DParticleComboStep* dStep0Wrapper;
		DBeamParticle* dComboBeamWrapper;
		DChargedTrackHypothesis* dProtonWrapper;

		//Step 1
		DParticleComboStep* dStep1Wrapper;
		DNeutralParticleHypothesis* dPhoton1Wrapper;
		DNeutralParticleHypothesis* dPhoton2Wrapper;

		//Step 2
		DParticleComboStep* dStep2Wrapper;
		DNeutralParticleHypothesis* dPhoton3Wrapper;
		DNeutralParticleHypothesis* dPhoton4Wrapper;

		// DEFINE YOUR HISTOGRAMS HERE
		// EXAMPLES:
		TH1I* dHist_MissingMassSquared;
		TH1I* dHist_BeamEnergy;

		unsigned long int runNo;
		unsigned long int eventNo;
		double weight;
		double sqrts_;
		double abst_;
		double massEtaPi0_;
		double massEta_;
		double massPi0_;
		double sqrts;
		double abst;
		double massEtaPi0;
		double massEta;
		double massPi0;
		double phiR_;
		double thetaGJ_;
		double phiGJ_;
		double thetaEta_;
		double phiEta_;
		double thetaPi0_;
		double phiPi0_;
		double thetaGJ;
		double phiR;
		double phiGJ;
		double thetaEta;
		double phiEta;
		double thetaPi0;
		double phiPi0;
		unsigned int momentsGJ;
		unsigned int momentsEta;
		unsigned int momentsPi0;
		double YmomGJ[169]; // L=0...12
		double YmomEta[91]; // L=0,2,4,6,8,10,12
		double YmomPi0[91]; // L=0,2,4,6,8,10,12
		double YmomGJ_[169]; // L=0...12
		double YmomEta_[91]; // L=0,2,4,6,8,10,12
		double YmomPi0_[91]; // L=0,2,4,6,8,10,12

		double angular_moment(int L, int M, double theta, double phi) const;
		double model1_moment(int L, int M) const;
		std::vector<double> model1_amplitude(int L, int M) const;
		double model1_density(int source=0) const;
		void model1_use_generated(int yesno=1);
		double model1_Lmax() const;

 private:
    int model1_generated_kinematics;

	ClassDef(DSelector_etapi0_moments, 0);
};

void DSelector_etapi0_moments::Get_ComboWrappers(void)
{
	//Step 0
	dStep0Wrapper = dComboWrapper->Get_ParticleComboStep(0);
	dComboBeamWrapper = static_cast<DBeamParticle*>(dStep0Wrapper->Get_InitialParticle());
	dProtonWrapper = static_cast<DChargedTrackHypothesis*>(dStep0Wrapper->Get_FinalParticle(2));

	//Step 1
	dStep1Wrapper = dComboWrapper->Get_ParticleComboStep(1);
	dPhoton1Wrapper = static_cast<DNeutralParticleHypothesis*>(dStep1Wrapper->Get_FinalParticle(0));
	dPhoton2Wrapper = static_cast<DNeutralParticleHypothesis*>(dStep1Wrapper->Get_FinalParticle(1));

	//Step 2
	dStep2Wrapper = dComboWrapper->Get_ParticleComboStep(2);
	dPhoton3Wrapper = static_cast<DNeutralParticleHypothesis*>(dStep2Wrapper->Get_FinalParticle(0));
	dPhoton4Wrapper = static_cast<DNeutralParticleHypothesis*>(dStep2Wrapper->Get_FinalParticle(1));
}

#endif // DSelector_etapi0_moments_h
