{
	gSystem->Load("${ROOTSYS}/lib/libMathMore.so");
	gROOT->ProcessLine(".L trial_model.C+O");
}
