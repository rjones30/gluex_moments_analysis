#define HALLD_RECON_HOME "/cvmfs/oasis.opensciencegrid.org/gluex/halld_recon/Linux_CentOS7-x86_64-gcc4.8.5"
#define ROOT_ANALYSIS_HOME "/cvmfs/oasis.opensciencegrid.org/gluex/group/halld/Software/builds/Linux_CentOS7-x86_64-gcc4.8.5-cntr/gluex_root_analysis/gluex_root_analysis-1.12.0^hr4190/Linux_CentOS7-x86_64-gcc4.8.5-cntr"

{
	gSystem->AddIncludePath("-I" HALLD_RECON_HOME "/include/");
	gSystem->AddIncludePath("-I" ROOT_ANALYSIS_HOME "/include/");
	gSystem->Load(ROOT_ANALYSIS_HOME "/lib/libDSelector.so");
	gSystem->Load("${ROOTSYS}/lib/libMathMore.so");
	gROOT->ProcessLine(".L trial_model.C+O");
	gROOT->ProcessLine(".L DSelector_etapi0_moments.C+O");
}
