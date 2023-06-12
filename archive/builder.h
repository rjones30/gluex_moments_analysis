//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Tue Jun  6 11:50:38 2023 by ROOT version 6.22/06
// from TTree etapi0_moments/etapi0_moments
// found on file: root://cn440.storrs.hpc.uconn.edu/Gluex/resilient/simulation/moments-6-2023/etapi0_moments.root
//////////////////////////////////////////////////////////

#ifndef builder_h
#define builder_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <TSelector.h>
#include <TTreeReader.h>
#include <TTreeReaderValue.h>
#include <TTreeReaderArray.h>

// Headers needed by this particular selector
#include "TLorentzVector.h"



class builder : public TSelector {
public :
   TTreeReader     fReader;  //!the tree reader
   TTree          *fChain = 0;   //!pointer to the analyzed TTree or TChain
   TFile          *fTreeFile = 0;
   unsigned int   fNevents = 0;

   // Readers to access the data (delete the ones you do not need).
   //TTreeReaderValue<UInt_t> run = {fReader, "run"};
   //TTreeReaderValue<ULong64_t> event = {fReader, "event"};
   //TTreeReaderValue<Float_t> weight = {fReader, "weight"};
   //TTreeReaderValue<ULong64_t> numtruepid_final = {fReader, "numtruepid_final"};
   //TTreeReaderValue<ULong64_t> truepids_decay = {fReader, "truepids_decay"};
   //TTreeReaderValue<Bool_t> is_truetop = {fReader, "is_truetop"};
   //TTreeReaderValue<Bool_t> is_truecombo = {fReader, "is_truecombo"};
   //TTreeReaderValue<Bool_t> is_bdtcombo = {fReader, "is_bdtcombo"};
   //TTreeReaderValue<Bool_t> rftime = {fReader, "rftime"};
   //TTreeReaderValue<Float_t> kin_chisq = {fReader, "kin_chisq"};
   //TTreeReaderValue<UInt_t> kin_ndf = {fReader, "kin_ndf"};
   //TTreeReaderValue<UInt_t> beam_beamid = {fReader, "beam_beamid"};
   //TTreeReaderValue<Bool_t> beam_isgen = {fReader, "beam_isgen"};
   //TTreeReaderValue<TLorentzVector> beam_x4_meas = {fReader, "beam_x4_meas"};
   //TTreeReaderValue<TLorentzVector> beam_p4_meas = {fReader, "beam_p4_meas"};
   //TTreeReaderValue<TLorentzVector> beam_x4_kin = {fReader, "beam_x4_kin"};
   //TTreeReaderValue<TLorentzVector> beam_p4_kin = {fReader, "beam_p4_kin"};
   //TTreeReaderValue<TLorentzVector> beam_x4_true = {fReader, "beam_x4_true"};
   //TTreeReaderValue<TLorentzVector> beam_p4_true = {fReader, "beam_p4_true"};
   //TTreeReaderValue<UInt_t> p_trkid = {fReader, "p_trkid"};
   //TTreeReaderValue<TLorentzVector> p_x4_meas = {fReader, "p_x4_meas"};
   //TTreeReaderValue<TLorentzVector> p_p4_meas = {fReader, "p_p4_meas"};
   //TTreeReaderValue<TLorentzVector> p_x4_kin = {fReader, "p_x4_kin"};
   //TTreeReaderValue<TLorentzVector> p_p4_kin = {fReader, "p_p4_kin"};
   //TTreeReaderValue<Float_t> p_true_fom = {fReader, "p_true_fom"};
   //TTreeReaderValue<TLorentzVector> p_x4_true = {fReader, "p_x4_true"};
   //TTreeReaderValue<TLorentzVector> p_p4_true = {fReader, "p_p4_true"};
   //TTreeReaderValue<Float_t> p_pid_fom = {fReader, "p_pid_fom"};
   //TTreeReaderValue<Float_t> p_beta_time = {fReader, "p_beta_time"};
   //TTreeReaderValue<Float_t> p_chisq_time = {fReader, "p_chisq_time"};
   //TTreeReaderValue<UInt_t> p_ndf_time = {fReader, "p_ndf_time"};
   //TTreeReaderValue<UInt_t> p_ndf_trk = {fReader, "p_ndf_trk"};
   //TTreeReaderValue<Float_t> p_chisq_trk = {fReader, "p_chisq_trk"};
   //TTreeReaderValue<UInt_t> p_ndf_dedx = {fReader, "p_ndf_dedx"};
   //TTreeReaderValue<Float_t> p_chisq_dedx = {fReader, "p_chisq_dedx"};
   //TTreeReaderValue<Float_t> p_dedx_cdc = {fReader, "p_dedx_cdc"};
   //TTreeReaderValue<Float_t> p_dedx_fdc = {fReader, "p_dedx_fdc"};
   //TTreeReaderValue<Float_t> p_dedx_tof = {fReader, "p_dedx_tof"};
   //TTreeReaderValue<Float_t> p_dedx_st = {fReader, "p_dedx_st"};
   //TTreeReaderValue<Float_t> p_ebcal = {fReader, "p_ebcal"};
   //TTreeReaderValue<Float_t> p_eprebcal = {fReader, "p_eprebcal"};
   //TTreeReaderValue<Float_t> p_efcal = {fReader, "p_efcal"};
   //TTreeReaderValue<Float_t> p_bcal_delphi = {fReader, "p_bcal_delphi"};
   //TTreeReaderValue<Float_t> p_bcal_delz = {fReader, "p_bcal_delz"};
   //TTreeReaderValue<Float_t> p_fcal_doca = {fReader, "p_fcal_doca"};
   //TTreeReaderValue<UInt_t> g1_showid = {fReader, "g1_showid"};
   //TTreeReaderValue<TLorentzVector> g1_x4_meas = {fReader, "g1_x4_meas"};
   //TTreeReaderValue<TLorentzVector> g1_p4_meas = {fReader, "g1_p4_meas"};
   //TTreeReaderValue<TLorentzVector> g1_x4_kin = {fReader, "g1_x4_kin"};
   //TTreeReaderValue<TLorentzVector> g1_p4_kin = {fReader, "g1_p4_kin"};
   //TTreeReaderValue<Float_t> g1_true_fom = {fReader, "g1_true_fom"};
   //TTreeReaderValue<TLorentzVector> g1_x4_true = {fReader, "g1_x4_true"};
   //TTreeReaderValue<TLorentzVector> g1_p4_true = {fReader, "g1_p4_true"};
   //TTreeReaderValue<Float_t> g1_beta_time = {fReader, "g1_beta_time"};
   //TTreeReaderValue<Float_t> g1_chisq_time = {fReader, "g1_chisq_time"};
   //TTreeReaderValue<UInt_t> g1_ndf_time = {fReader, "g1_ndf_time"};
   //TTreeReaderValue<Float_t> g1_ebcal = {fReader, "g1_ebcal"};
   //TTreeReaderValue<Float_t> g1_eprebcal = {fReader, "g1_eprebcal"};
   //TTreeReaderValue<Float_t> g1_efcal = {fReader, "g1_efcal"};
   //TTreeReaderValue<Float_t> g1_bcal_delphi = {fReader, "g1_bcal_delphi"};
   //TTreeReaderValue<Float_t> g1_bcal_delz = {fReader, "g1_bcal_delz"};
   //TTreeReaderValue<Float_t> g1_fcal_doca = {fReader, "g1_fcal_doca"};
   //TTreeReaderValue<UInt_t> g2_showid = {fReader, "g2_showid"};
   //TTreeReaderValue<TLorentzVector> g2_x4_meas = {fReader, "g2_x4_meas"};
   //TTreeReaderValue<TLorentzVector> g2_p4_meas = {fReader, "g2_p4_meas"};
   //TTreeReaderValue<TLorentzVector> g2_x4_kin = {fReader, "g2_x4_kin"};
   //TTreeReaderValue<TLorentzVector> g2_p4_kin = {fReader, "g2_p4_kin"};
   //TTreeReaderValue<Float_t> g2_true_fom = {fReader, "g2_true_fom"};
   //TTreeReaderValue<TLorentzVector> g2_x4_true = {fReader, "g2_x4_true"};
   //TTreeReaderValue<TLorentzVector> g2_p4_true = {fReader, "g2_p4_true"};
   //TTreeReaderValue<Float_t> g2_beta_time = {fReader, "g2_beta_time"};
   //TTreeReaderValue<Float_t> g2_chisq_time = {fReader, "g2_chisq_time"};
   //TTreeReaderValue<UInt_t> g2_ndf_time = {fReader, "g2_ndf_time"};
   //TTreeReaderValue<Float_t> g2_ebcal = {fReader, "g2_ebcal"};
   //TTreeReaderValue<Float_t> g2_eprebcal = {fReader, "g2_eprebcal"};
   //TTreeReaderValue<Float_t> g2_efcal = {fReader, "g2_efcal"};
   //TTreeReaderValue<Float_t> g2_bcal_delphi = {fReader, "g2_bcal_delphi"};
   //TTreeReaderValue<Float_t> g2_bcal_delz = {fReader, "g2_bcal_delz"};
   //TTreeReaderValue<Float_t> g2_fcal_doca = {fReader, "g2_fcal_doca"};
   //TTreeReaderValue<UInt_t> g3_showid = {fReader, "g3_showid"};
   //TTreeReaderValue<TLorentzVector> g3_x4_meas = {fReader, "g3_x4_meas"};
   //TTreeReaderValue<TLorentzVector> g3_p4_meas = {fReader, "g3_p4_meas"};
   //TTreeReaderValue<TLorentzVector> g3_x4_kin = {fReader, "g3_x4_kin"};
   //TTreeReaderValue<TLorentzVector> g3_p4_kin = {fReader, "g3_p4_kin"};
   //TTreeReaderValue<Float_t> g3_true_fom = {fReader, "g3_true_fom"};
   //TTreeReaderValue<TLorentzVector> g3_x4_true = {fReader, "g3_x4_true"};
   //TTreeReaderValue<TLorentzVector> g3_p4_true = {fReader, "g3_p4_true"};
   //TTreeReaderValue<Float_t> g3_beta_time = {fReader, "g3_beta_time"};
   //TTreeReaderValue<Float_t> g3_chisq_time = {fReader, "g3_chisq_time"};
   //TTreeReaderValue<UInt_t> g3_ndf_time = {fReader, "g3_ndf_time"};
   //TTreeReaderValue<Float_t> g3_ebcal = {fReader, "g3_ebcal"};
   //TTreeReaderValue<Float_t> g3_eprebcal = {fReader, "g3_eprebcal"};
   //TTreeReaderValue<Float_t> g3_efcal = {fReader, "g3_efcal"};
   //TTreeReaderValue<Float_t> g3_bcal_delphi = {fReader, "g3_bcal_delphi"};
   //TTreeReaderValue<Float_t> g3_bcal_delz = {fReader, "g3_bcal_delz"};
   //TTreeReaderValue<Float_t> g3_fcal_doca = {fReader, "g3_fcal_doca"};
   //TTreeReaderValue<UInt_t> g4_showid = {fReader, "g4_showid"};
   //TTreeReaderValue<TLorentzVector> g4_x4_meas = {fReader, "g4_x4_meas"};
   //TTreeReaderValue<TLorentzVector> g4_p4_meas = {fReader, "g4_p4_meas"};
   //TTreeReaderValue<TLorentzVector> g4_x4_kin = {fReader, "g4_x4_kin"};
   //TTreeReaderValue<TLorentzVector> g4_p4_kin = {fReader, "g4_p4_kin"};
   //TTreeReaderValue<Float_t> g4_true_fom = {fReader, "g4_true_fom"};
   //TTreeReaderValue<TLorentzVector> g4_x4_true = {fReader, "g4_x4_true"};
   //TTreeReaderValue<TLorentzVector> g4_p4_true = {fReader, "g4_p4_true"};
   //TTreeReaderValue<Float_t> g4_beta_time = {fReader, "g4_beta_time"};
   //TTreeReaderValue<Float_t> g4_chisq_time = {fReader, "g4_chisq_time"};
   //TTreeReaderValue<UInt_t> g4_ndf_time = {fReader, "g4_ndf_time"};
   //TTreeReaderValue<Float_t> g4_ebcal = {fReader, "g4_ebcal"};
   //TTreeReaderValue<Float_t> g4_eprebcal = {fReader, "g4_eprebcal"};
   //TTreeReaderValue<Float_t> g4_efcal = {fReader, "g4_efcal"};
   //TTreeReaderValue<Float_t> g4_bcal_delphi = {fReader, "g4_bcal_delphi"};
   //TTreeReaderValue<Float_t> g4_bcal_delz = {fReader, "g4_bcal_delz"};
   //TTreeReaderValue<Float_t> g4_fcal_doca = {fReader, "g4_fcal_doca"};
   //TTreeReaderValue<Long64_t> runNo = {fReader, "runNo"};
   //TTreeReaderValue<Long64_t> eventNo = {fReader, "eventNo"};
   //TTreeReaderValue<Double_t> weight_ = {fReader, "weight_"};
   //TTreeReaderValue<Double_t> sqrts_ = {fReader, "sqrts_"};
   //TTreeReaderValue<Double_t> abst_ = {fReader, "abst_"};
   //TTreeReaderValue<Double_t> massEta_ = {fReader, "massEta_"};
   //TTreeReaderValue<Double_t> massPi0_ = {fReader, "massPi0_"};
   //TTreeReaderValue<Double_t> sqrts = {fReader, "sqrts"};
   //TTreeReaderValue<Double_t> abst = {fReader, "abst"};
   //TTreeReaderValue<Double_t> massEta = {fReader, "massEta"};
   //TTreeReaderValue<Double_t> massPi0 = {fReader, "massPi0"};
   //TTreeReaderValue<Double_t> phiR_ = {fReader, "phiR_"};
   //TTreeReaderValue<Double_t> thetaGJ_ = {fReader, "thetaGJ_"};
   //TTreeReaderValue<Double_t> phiGJ_ = {fReader, "phiGJ_"};
   //TTreeReaderValue<Double_t> thetaEta_ = {fReader, "thetaEta_"};
   //TTreeReaderValue<Double_t> phiEta_ = {fReader, "phiEta_"};
   //TTreeReaderValue<Double_t> thetaPi0_ = {fReader, "thetaPi0_"};
   //TTreeReaderValue<Double_t> phiPi0_ = {fReader, "phiPi0_"};
   //TTreeReaderValue<Double_t> phiR = {fReader, "phiR"};
   //TTreeReaderValue<Double_t> thetaGJ = {fReader, "thetaGJ"};
   //TTreeReaderValue<Double_t> phiGJ = {fReader, "phiGJ"};
   //TTreeReaderValue<Double_t> thetaEta = {fReader, "thetaEta"};
   //TTreeReaderValue<Double_t> phiEta = {fReader, "phiEta"};
   //TTreeReaderValue<Double_t> thetaPi0 = {fReader, "thetaPi0"};
   //TTreeReaderValue<Double_t> phiPi0 = {fReader, "phiPi0"};
   //TTreeReaderValue<Int_t> momentsGJ = {fReader, "momentsGJ"};
   //TTreeReaderValue<Int_t> momentsEta = {fReader, "momentsEta"};
   //TTreeReaderValue<Int_t> momentsPi0 = {fReader, "momentsPi0"};
   //TTreeReaderArray<Double_t> YmomGJ = {fReader, "YmomGJ"};
   //TTreeReaderArray<Double_t> YmomEta = {fReader, "YmomEta"};
   //TTreeReaderArray<Double_t> YmomPi0 = {fReader, "YmomPi0"};
   //TTreeReaderArray<Double_t> YmomGJ_ = {fReader, "YmomGJ_"};
   //TTreeReaderArray<Double_t> YmomEta_ = {fReader, "YmomEta_"};
   //TTreeReaderArray<Double_t> YmomPi0_ = {fReader, "YmomPi0_"};
 
   Long64_t runNo;
   Long64_t eventNo;
   Double_t weight_;
   Double_t sqrts_;
   Double_t abst_;
   Double_t massEta_;
   Double_t massPi0_;
   Double_t sqrts;
   Double_t abst;
   Double_t massEta;
   Double_t massPi0;
   Double_t phiR_;
   Double_t thetaGJ_;
   Double_t phiGJ_;
   Double_t thetaEta_;
   Double_t phiEta_;
   Double_t thetaPi0_;
   Double_t phiPi0_;
   Double_t phiR;
   Double_t thetaGJ;
   Double_t phiGJ;
   Double_t thetaEta;
   Double_t phiEta;
   Double_t thetaPi0;
   Double_t phiPi0;
   Int_t momentsGJ;
   Int_t momentsEta;
   Int_t momentsPi0;
   Double_t *YmomGJ;
   Double_t *YmomEta;
   Double_t *YmomPi0;
   Double_t *YmomGJ_;
   Double_t *YmomEta_;
   Double_t *YmomPi0_;

   builder(TTree * /*tree*/ =0) { }
   virtual ~builder() { }
   virtual Int_t   Version() const { return 2; }
   virtual void    Begin(TTree *tree);
   virtual void    SlaveBegin(TTree *tree);
   virtual void    Init(TTree *tree);
   virtual Bool_t  Notify();
   virtual Bool_t  Process(Long64_t entry);
   virtual Int_t   GetEntry(Long64_t entry, Int_t getall = 0) { return fChain ? fChain->GetTree()->GetEntry(entry, getall) : 0; }
   virtual void    SetOption(const char *option) { fOption = option; }
   virtual void    SetObject(TObject *obj) { fObject = obj; }
   virtual void    SetInputList(TList *input) { fInput = input; }
   virtual TList  *GetOutputList() const { return fOutput; }
   virtual void    SlaveTerminate();
   virtual void    Terminate();

   ClassDef(builder,0);

};

#endif

#ifdef builder_cxx
void builder::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the reader is initialized.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   //fReader.SetTree(tree);
}

Bool_t builder::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}


#endif // #ifdef builder_cxx
