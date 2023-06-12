#define builder_cxx
// The class definition in builder.h has been generated automatically
// by the ROOT utility TTree::MakeSelector(). This class is derived
// from the ROOT class TSelector. For more information on the TSelector
// framework see $ROOTSYS/README/README.SELECTOR or the ROOT User Manual.


// The following methods are defined in this file:
//    Begin():        called every time a loop on the tree starts,
//                    a convenient place to create your histograms.
//    SlaveBegin():   called after Begin(), when on PROOF called only on the
//                    slave servers.
//    Process():      called for each event, in this function you decide what
//                    to read and fill your histograms.
//    SlaveTerminate: called at the end of the loop on the tree, when on PROOF
//                    called only on the slave servers.
//    Terminate():    called at the end of the loop on the tree,
//                    a convenient place to draw/fit your histograms.
//
// To use this file, try the following session on your Tree T:
//
// root> T->Process("builder.C")
// root> T->Process("builder.C","some options")
// root> T->Process("builder.C+")
//

#define TREENAME "etapi0_moments"
#define TREEFILE "etapi0_moments.root"
#define XROOTDURL "root://cn440.storrs.hpc.uconn.edu/Gluex/resilient/"
#define TREEDIR "simulation/moments-6-2023/"

#include "builder.h"
#include <TH2.h>
#include <TStyle.h>

void builder::Begin(TTree * /*tree*/)
{
   // The Begin() function is called at the start of the query.
   // When running with PROOF Begin() is only called on the client.
   // The tree argument is deprecated (on PROOF 0 is passed).

   TString option = GetOption();
}

void builder::SlaveBegin(TTree * /*tree*/)
{
   // The SlaveBegin() function is called after the Begin() function.
   // When running with PROOF SlaveBegin() is called on each slave server.
   // The tree argument is deprecated (on PROOF 0 is passed).

   TString option = GetOption();

   fTreeFile = TFile::Open(XROOTDURL TREEDIR TREEFILE);
   fChain = (TTree*)fTreeFile->Get(TREENAME);
   fNevents = fChain->GetEntries();

   fChain->SetBranchAddress("momentsGJ", &momentsGJ);
   fChain->SetBranchAddress("momentsEta", &momentsEta);
   fChain->SetBranchAddress("momentsPi0", &momentsPi0);
   fChain->GetEntry(1);
   YmomGJ = new double [fNevents * momentsGJ];
   YmomGJ_ = new double [fNevents * momentsGJ];
   YmomEta = new double [fNevents * momentsEta];
   YmomEta_ = new double [fNevents * momentsEta];
   YmomPi0 = new double [fNevents * momentsPi0];
   YmomPi0_ = new double [fNevents * momentsPi0];
}

Bool_t builder::Process(Long64_t entry)
{
   // The Process() function is called for each entry in the tree (or possibly
   // keyed object in the case of PROOF) to be processed. The entry argument
   // specifies which entry in the currently loaded tree is to be processed.
   // When processing keyed objects with PROOF, the object is already loaded
   // and is available via the fObject pointer.
   //
   // This function should contain the \"body\" of the analysis. It can contain
   // simple or elaborate selection criteria, run algorithms on the data
   // of the event and typically fill histograms.
   //
   // The processing can be stopped by calling Abort().
   //
   // Use fStatus to set the return value of TTree::Process().
   //
   // The return value is currently not used.

   //fReader.SetLocalEntry(entry);

   fChain->SetBranchAddress("YmomGJ", &YmomGJ[entry * momentsGJ]);
   fChain->SetBranchAddress("YmomGJ_", &YmomGJ_[entry * momentsGJ]);
   fChain->SetBranchAddress("YmomEta", &YmomEta[entry * momentsEta]);
   fChain->SetBranchAddress("YmomEta_", &YmomEta_[entry * momentsEta]);
   fChain->SetBranchAddress("YmomPi0", &YmomPi0[entry * momentsPi0]);
   fChain->SetBranchAddress("YmomPi0_", &YmomPi0_[entry * momentsPi0]);
   fChain->GetEntry(entry);
   if (entry % 10000 == 0) {
      std::cout << "read " << entry << " events from input sample" << std::endl;
   }

   return kTRUE;
}

void builder::SlaveTerminate()
{
   // The SlaveTerminate() function is called after all entries or objects
   // have been processed. When running with PROOF SlaveTerminate() is called
   // on each slave server.

}

void builder::Terminate()
{
   // The Terminate() function is the last function to be called during
   // a query. It always runs on the client, it can be used to present
   // the results graphically or save the results to file.

}
