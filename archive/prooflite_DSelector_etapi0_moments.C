int NumThreads = 6; // your choice of number of threads to use
char input_file[] = "root://nod25.phys.uconn.edu/Gluex/simulation/moments-6-2023/tree_etapi0__B4_T1_S1_M7_M17_F4.root";
char dselector[] = "DSelector_etapi0_moments.C++";
 
void do_analysis(){
#include "TProof.h"
#include "TProofDebug.h"
  R__LOAD_LIBRARY(libDSelector);
  gROOT->ProcessLine(".x $ROOT_ANALYSIS_HOME/scripts/Load_DSelector.C");
  DPROOFLiteManager *dproof = new DPROOFLiteManager();
  dproof->Process_Tree(input_file, dselector, NumThreads);
}
