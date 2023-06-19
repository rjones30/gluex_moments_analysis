#define moments_plotter_cxx
// The class definition in moments_plotter.h has been generated automatically
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
// root> T->Process("moments_plotter.C")
// root> T->Process("moments_plotter.C","some options")
// root> T->Process("moments_plotter.C+")
//


#include "moments_plotter.h"
#include <TH2.h>
#include <TStyle.h>

#include <sstream>

void moments_plotter::Begin(TTree * /*tree*/)
{
   // The Begin() function is called at the start of the query.
   // When running with PROOF Begin() is only called on the client.
   // The tree argument is deprecated (on PROOF 0 is passed).

   TString option = GetOption();
}

void moments_plotter::SlaveBegin(TTree *tree)
{
   // The SlaveBegin() function is called after the Begin() function.
   // When running with PROOF SlaveBegin() is called on each slave server.
   // The tree argument is deprecated (on PROOF 0 is passed).

   TString option = GetOption();

   int imom = 0;
   for (int L=0; L <= 6; ++L) {
      for (int M=-L; M <= L; ++M, ++imom) {
         for (int it=0; it < 3; ++it) {
            std::stringstream name, title;
            name << "h" << it << "m" << imom;
            title << "moment " << imom << ", t bin " << it;
            hmoment[it][imom] = new TProfile(name.str().c_str(), title.str().c_str(), 250, 0, 2.5);
         }
      }
   }

   tree->GetBranch("massEtaPi0")->SetAddress(&massEtaPi0);
   //tree->GetBranch("massEtaPi0_")->SetAddress(&massEtaPi0_);
   tree->GetBranch("abst")->SetAddress(&abst);
   //tree->GetBranch("abst_")->SetAddress(&abst_);
   tree->GetBranch("momentsGJ")->SetAddress(&momentsGJ);
   tree->GetBranch("model1moment")->SetAddress(model1moment);
   //tree->GetBranch("model1moment_")->SetAddress(model1moment_);
   fChain = tree;
}

Bool_t moments_plotter::Process(Long64_t entry)
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

   //fReader.SetEntry(entry);
   fChain->GetEntry(entry);

   int it;
   if (fabs(abst - 0.2) < 0.01)
      it = 0;
   else if (fabs(abst - 0.5) < 0.01)
      it = 1;
   else if (fabs(abst - 0.9) < 0.01)
      it = 2;
   else
      return kFALSE;
   for (int im=0; im < 49; ++im)
      hmoment[it][im]->Fill(massEtaPi0, model1moment[im]);
   return kTRUE;
}

void moments_plotter::SlaveTerminate()
{
   // The SlaveTerminate() function is called after all entries or objects
   // have been processed. When running with PROOF SlaveTerminate() is called
   // on each slave server.

   TFile fout("moments_plots.root", "recreate");
   for (int it=0; it < 3; ++it)
      for (int im=0; im < 49; ++im)
         hmoment[it][im]->Write();
}

void moments_plotter::Terminate()
{
   // The Terminate() function is the last function to be called during
   // a query. It always runs on the client, it can be used to present
   // the results graphically or save the results to file.

}
