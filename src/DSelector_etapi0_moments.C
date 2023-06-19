#include "DSelector_etapi0_moments.h"
#include <Math/SpecFuncMathMore.h>
#include <vector>

void DSelector_etapi0_moments::Init(TTree *locTree)
{
	// USERS: IN THIS FUNCTION, ONLY MODIFY SECTIONS WITH A "USER" OR "EXAMPLE" LABEL. LEAVE THE REST ALONE.

	// The Init() function is called when the selector needs to initialize a new tree or chain.
	// Typically here the branch addresses and branch pointers of the tree will be set.
	// Init() will be called many times when running on PROOF (once per file to be processed).

	//USERS: SET OUTPUT FILE NAME //can be overriden by user in PROOF
	dOutputFileName = ""; //"" for none
	dOutputTreeFileName = ""; //"" for none
	dFlatTreeFileName = "etapi0_moments.root"; //output flat tree (one combo per tree entry), "" for none
	dFlatTreeName = "etapi0_moments"; //if blank, default name will be chosen
	//dSaveDefaultFlatBranches = true; // False: don't save default branches, reduce disk footprint.
	//dSaveTLorentzVectorsAsFundamentaFlatTree = false; // Default (or false): save particles as TLorentzVector objects. True: save as four doubles instead.

	//Because this function gets called for each TTree in the TChain, we must be careful:
		//We need to re-initialize the tree interface & branch wrappers, but don't want to recreate histograms
	bool locInitializedPriorFlag = dInitializedFlag; //save whether have been initialized previously
	DSelector::Init(locTree); //This must be called to initialize wrappers for each new TTree
	//gDirectory now points to the output file with name dOutputFileName (if any)
	if(locInitializedPriorFlag)
		return; //have already created histograms, etc. below: exit

	Get_ComboWrappers();
	dPreviousRunNumber = 0;

	/*********************************** EXAMPLE USER INITIALIZATION: ANALYSIS ACTIONS **********************************/

	// EXAMPLE: Create deque for histogramming particle masses:
	// // For histogramming the phi mass in phi -> K+ K-
	// // Be sure to change this and dAnalyzeCutActions to match reaction
	std::deque<Particle_t> MyPhi;
	MyPhi.push_back(KPlus); MyPhi.push_back(KMinus);

	//ANALYSIS ACTIONS: //Executed in order if added to dAnalysisActions
	//false/true below: use measured/kinfit data

	//PID
	dAnalysisActions.push_back(new DHistogramAction_ParticleID(dComboWrapper, false));
	//below: value: +/- N ns, Unknown: All PIDs, SYS_NULL: all timing systems
	//dAnalysisActions.push_back(new DCutAction_PIDDeltaT(dComboWrapper, false, 0.5, KPlus, SYS_BCAL));

	//PIDFOM (for charged tracks)
	//dAnalysisActions.push_back(new DHistogramAction_PIDFOM(dComboWrapper));
	//dAnalysisActions.push_back(new DCutAction_PIDFOM(dComboWrapper, KPlus, 0.1));
	//dAnalysisActions.push_back(new DCutAction_EachPIDFOM(dComboWrapper, 0.1));

	//MASSES
	//dAnalysisActions.push_back(new DHistogramAction_InvariantMass(dComboWrapper, false, Lambda, 1000, 1.0, 1.2, "Lambda"));
	//dAnalysisActions.push_back(new DHistogramAction_MissingMassSquared(dComboWrapper, false, 1000, -0.1, 0.1));

	//KINFIT RESULTS
	dAnalysisActions.push_back(new DHistogramAction_KinFitResults(dComboWrapper));

	//CUT MISSING MASS
	//dAnalysisActions.push_back(new DCutAction_MissingMassSquared(dComboWrapper, false, -0.03, 0.02));

	//CUT ON SHOWER QUALITY
	//dAnalysisActions.push_back(new DCutAction_ShowerQuality(dComboWrapper, SYS_FCAL, 0.5));

	//BEAM ENERGY
	dAnalysisActions.push_back(new DHistogramAction_BeamEnergy(dComboWrapper, false));
	//dAnalysisActions.push_back(new DCutAction_BeamEnergy(dComboWrapper, false, 8.2, 8.8));  // Coherent peak for runs in the range 30000-59999

	//KINEMATICS
	dAnalysisActions.push_back(new DHistogramAction_ParticleComboKinematics(dComboWrapper, false));

	// ANALYZE CUT ACTIONS
	// // Change MyPhi to match reaction
	dAnalyzeCutActions = new DHistogramAction_AnalyzeCutActions( dAnalysisActions, dComboWrapper, false, 0, MyPhi, 1000, 0.9, 2.4, "CutActionEffect" );

	//INITIALIZE ACTIONS
	//If you create any actions that you want to run manually (i.e. don't add to dAnalysisActions), be sure to initialize them here as well
	Initialize_Actions();
	dAnalyzeCutActions->Initialize(); // manual action, must call Initialize()

	/******************************** EXAMPLE USER INITIALIZATION: STAND-ALONE HISTOGRAMS *******************************/

	//EXAMPLE MANUAL HISTOGRAMS:
	dHist_MissingMassSquared = new TH1I("MissingMassSquared", ";Missing Mass Squared (GeV/c^{2})^{2}", 600, -0.06, 0.06);
	dHist_BeamEnergy = new TH1I("BeamEnergy", ";Beam Energy (GeV)", 600, 0.0, 12.0);

	/************************** EXAMPLE USER INITIALIZATION: CUSTOM OUTPUT BRANCHES - MAIN TREE *************************/

	//EXAMPLE MAIN TREE CUSTOM BRANCHES (OUTPUT ROOT FILE NAME MUST FIRST BE GIVEN!!!! (ABOVE: TOP)):
	//The type for the branch must be included in the brackets
	//1st function argument is the name of the branch
	//2nd function argument is the name of the branch that contains the size of the array (for fundamentals only)
	/*
	dTreeInterface->Create_Branch_Fundamental<Int_t>("my_int"); //fundamental = char, int, float, double, etc.
	dTreeInterface->Create_Branch_FundamentalArray<Int_t>("my_int_array", "my_int");
	dTreeInterface->Create_Branch_FundamentalArray<Float_t>("my_combo_array", "NumCombos");
	dTreeInterface->Create_Branch_NoSplitTObject<TLorentzVector>("my_p4");
	dTreeInterface->Create_Branch_ClonesArray<TLorentzVector>("my_p4_array");
	*/

	/************************** EXAMPLE USER INITIALIZATION: CUSTOM OUTPUT BRANCHES - FLAT TREE *************************/

	// RECOMMENDED: CREATE ACCIDENTAL WEIGHT BRANCH
	// dFlatTreeInterface->Create_Branch_Fundamental<Double_t>("accidweight");

	//EXAMPLE FLAT TREE CUSTOM BRANCHES (OUTPUT ROOT FILE NAME MUST FIRST BE GIVEN!!!! (ABOVE: TOP)):
	//The type for the branch must be included in the brackets
	//1st function argument is the name of the branch
	//2nd function argument is the name of the branch that contains the size of the array (for fundamentals only)
	/*
	dFlatTreeInterface->Create_Branch_Fundamental<Int_t>("flat_my_int"); //fundamental = char, int, float, double, etc.
	dFlatTreeInterface->Create_Branch_FundamentalArray<Int_t>("flat_my_int_array", "flat_my_int");
	dFlatTreeInterface->Create_Branch_NoSplitTObject<TLorentzVector>("flat_my_p4");
	dFlatTreeInterface->Create_Branch_ClonesArray<TLorentzVector>("flat_my_p4_array");
	*/

	dFlatTreeInterface->Create_Branch_Fundamental<Long64_t>("runNo");
	dFlatTreeInterface->Create_Branch_Fundamental<Long64_t>("eventNo");
	dFlatTreeInterface->Create_Branch_Fundamental<Double_t>("weight_");
	dFlatTreeInterface->Create_Branch_Fundamental<Double_t>("sqrts_");
	dFlatTreeInterface->Create_Branch_Fundamental<Double_t>("abst_");
	dFlatTreeInterface->Create_Branch_Fundamental<Double_t>("massEtaPi0_");
	dFlatTreeInterface->Create_Branch_Fundamental<Double_t>("massEta_");
	dFlatTreeInterface->Create_Branch_Fundamental<Double_t>("massPi0_");
	dFlatTreeInterface->Create_Branch_Fundamental<Double_t>("sqrts");
	dFlatTreeInterface->Create_Branch_Fundamental<Double_t>("abst");
	dFlatTreeInterface->Create_Branch_Fundamental<Double_t>("massEtaPi0");
	dFlatTreeInterface->Create_Branch_Fundamental<Double_t>("massEta");
	dFlatTreeInterface->Create_Branch_Fundamental<Double_t>("massPi0");
	dFlatTreeInterface->Create_Branch_Fundamental<Double_t>("phiR_");
	dFlatTreeInterface->Create_Branch_Fundamental<Double_t>("thetaGJ_");
	dFlatTreeInterface->Create_Branch_Fundamental<Double_t>("phiGJ_");
	dFlatTreeInterface->Create_Branch_Fundamental<Double_t>("thetaEta_");
	dFlatTreeInterface->Create_Branch_Fundamental<Double_t>("phiEta_");
	dFlatTreeInterface->Create_Branch_Fundamental<Double_t>("thetaPi0_");
	dFlatTreeInterface->Create_Branch_Fundamental<Double_t>("phiPi0_");
	dFlatTreeInterface->Create_Branch_Fundamental<Double_t>("phiR");
	dFlatTreeInterface->Create_Branch_Fundamental<Double_t>("thetaGJ");
	dFlatTreeInterface->Create_Branch_Fundamental<Double_t>("phiGJ");
	dFlatTreeInterface->Create_Branch_Fundamental<Double_t>("thetaEta");
	dFlatTreeInterface->Create_Branch_Fundamental<Double_t>("phiEta");
	dFlatTreeInterface->Create_Branch_Fundamental<Double_t>("thetaPi0");
	dFlatTreeInterface->Create_Branch_Fundamental<Double_t>("phiPi0");
	dFlatTreeInterface->Create_Branch_Fundamental<Int_t>("momentsGJ");
	dFlatTreeInterface->Create_Branch_Fundamental<Int_t>("momentsEta");
	dFlatTreeInterface->Create_Branch_Fundamental<Int_t>("momentsPi0");
	dFlatTreeInterface->Create_Branch_FundamentalArray<Double_t>("YmomGJ", "momentsGJ");
	dFlatTreeInterface->Create_Branch_FundamentalArray<Double_t>("YmomEta", "momentsEta");
	dFlatTreeInterface->Create_Branch_FundamentalArray<Double_t>("YmomPi0", "momentsPi0");
	dFlatTreeInterface->Create_Branch_FundamentalArray<Double_t>("YmomGJ_", "momentsGJ");
	dFlatTreeInterface->Create_Branch_FundamentalArray<Double_t>("YmomEta_", "momentsEta");
	dFlatTreeInterface->Create_Branch_FundamentalArray<Double_t>("YmomPi0_", "momentsPi0");
	dFlatTreeInterface->Create_Branch_FundamentalArray<Double_t>("model1moment", "momentsGJ");
	dFlatTreeInterface->Create_Branch_FundamentalArray<Double_t>("model1moment_", "momentsGJ");

    model1_generated_kinematics = 0;

	/************************************* ADVANCED EXAMPLE: CHOOSE BRANCHES TO READ ************************************/

	//TO SAVE PROCESSING TIME
		//If you know you don't need all of the branches/data, but just a subset of it, you can speed things up
		//By default, for each event, the data is retrieved for all branches
		//If you know you only need data for some branches, you can skip grabbing data from the branches you don't need
		//Do this by doing something similar to the commented code below

	//dTreeInterface->Clear_GetEntryBranches(); //now get none
	//dTreeInterface->Register_GetEntryBranch("Proton__P4"); //manually set the branches you want

	/************************************** DETERMINE IF ANALYZING SIMULATED DATA *************************************/

	dIsMC = (dTreeInterface->Get_Branch("MCWeight") != NULL);

}

Bool_t DSelector_etapi0_moments::Process(Long64_t locEntry)
{
	// The Process() function is called for each entry in the tree. The entry argument
	// specifies which entry in the currently loaded tree is to be processed.
	//
	// This function should contain the "body" of the analysis. It can contain
	// simple or elaborate selection criteria, run algorithms on the data
	// of the event and typically fill histograms.
	//
	// The processing can be stopped by calling Abort().
	// Use fStatus to set the return value of TTree::Process().
	// The return value is currently not used.

	//CALL THIS FIRST
	DSelector::Process(locEntry); //Gets the data from the tree for the entry
	//cout << "RUN " << Get_RunNumber() << ", EVENT " << Get_EventNumber() << endl;
	//TLorentzVector locProductionX4 = Get_X4_Production();

	/******************************************** GET POLARIZATION ORIENTATION ******************************************/

	//Only if the run number changes
	//RCDB environment must be setup in order for this to work! (Will return false otherwise)
	UInt_t locRunNumber = Get_RunNumber();
	if(locRunNumber != dPreviousRunNumber)
	{
		dIsPolarizedFlag = dAnalysisUtilities.Get_IsPolarizedBeam(locRunNumber, dIsPARAFlag);
		dPreviousRunNumber = locRunNumber;
	}

	/********************************************* SETUP UNIQUENESS TRACKING ********************************************/

	//ANALYSIS ACTIONS: Reset uniqueness tracking for each action
	//For any actions that you are executing manually, be sure to call Reset_NewEvent() on them here
	Reset_Actions_NewEvent();
	dAnalyzeCutActions->Reset_NewEvent(); // manual action, must call Reset_NewEvent()

	//PREVENT-DOUBLE COUNTING WHEN HISTOGRAMMING
		//Sometimes, some content is the exact same between one combo and the next
			//e.g. maybe two combos have different beam particles, but the same data for the final-state
		//When histogramming, you don't want to double-count when this happens: artificially inflates your signal (or background)
		//So, for each quantity you histogram, keep track of what particles you used (for a given combo)
		//Then for each combo, just compare to what you used before, and make sure it's unique

	//EXAMPLE 1: Particle-specific info:
	set<Int_t> locUsedSoFar_BeamEnergy; //Int_t: Unique ID for beam particles. set: easy to use, fast to search

	//EXAMPLE 2: Combo-specific info:
		//In general: Could have multiple particles with the same PID: Use a set of Int_t's
		//In general: Multiple PIDs, so multiple sets: Contain within a map
		//Multiple combos: Contain maps within a set (easier, faster to search)
	set<map<Particle_t, set<Int_t> > > locUsedSoFar_MissingMass;

	//INSERT USER ANALYSIS UNIQUENESS TRACKING HERE

	/**************************************** EXAMPLE: FILL CUSTOM OUTPUT BRANCHES **************************************/

	/*
	Int_t locMyInt = 7;
	dTreeInterface->Fill_Fundamental<Int_t>("my_int", locMyInt);

	TLorentzVector locMyP4(4.0, 3.0, 2.0, 1.0);
	dTreeInterface->Fill_TObject<TLorentzVector>("my_p4", locMyP4);

	for(int loc_i = 0; loc_i < locMyInt; ++loc_i)
		dTreeInterface->Fill_Fundamental<Int_t>("my_int_array", 3*loc_i, loc_i); //2nd argument = value, 3rd = array index
	*/

	/************************************************* LOOP OVER COMBOS *************************************************/

	int best_combo(-1);
    double best_comboCL(0);

	//Loop over combos
	for(UInt_t loc_i = 0; loc_i < Get_NumCombos(); ++loc_i)
	{
		//Set branch array indices for combo and all combo particles
		dComboWrapper->Set_ComboIndex(loc_i);

		// Is used to indicate when combos have been cut
		if(dComboWrapper->Get_IsComboCut()) // Is false when tree originally created
			continue; // Combo has been cut previously

		/********************************************** GET PARTICLE INDICES *********************************************/

		//Used for tracking uniqueness when filling histograms, and for determining unused particles

		//Step 0
		Int_t locBeamID = dComboBeamWrapper->Get_BeamID();
		Int_t locProtonTrackID = dProtonWrapper->Get_TrackID();

		//Step 1
		Int_t locPhoton1NeutralID = dPhoton1Wrapper->Get_NeutralID();
		Int_t locPhoton2NeutralID = dPhoton2Wrapper->Get_NeutralID();

		//Step 2
		Int_t locPhoton3NeutralID = dPhoton3Wrapper->Get_NeutralID();
		Int_t locPhoton4NeutralID = dPhoton4Wrapper->Get_NeutralID();

		/*********************************************** GET FOUR-MOMENTUM **********************************************/

		// Get P4's: //is kinfit if kinfit performed, else is measured
		//dTargetP4 is target p4
		//Step 0
		TLorentzVector locBeamP4 = dComboBeamWrapper->Get_P4();
		TLorentzVector locProtonP4 = dProtonWrapper->Get_P4();
		//Step 1
		TLorentzVector locPhoton1P4 = dPhoton1Wrapper->Get_P4();
		TLorentzVector locPhoton2P4 = dPhoton2Wrapper->Get_P4();
		//Step 2
		TLorentzVector locPhoton3P4 = dPhoton3Wrapper->Get_P4();
		TLorentzVector locPhoton4P4 = dPhoton4Wrapper->Get_P4();

		// Get Measured P4's:
		//Step 0
		TLorentzVector locBeamP4_Measured = dComboBeamWrapper->Get_P4_Measured();
		TLorentzVector locProtonP4_Measured = dProtonWrapper->Get_P4_Measured();
		//Step 1
		TLorentzVector locPhoton1P4_Measured = dPhoton1Wrapper->Get_P4_Measured();
		TLorentzVector locPhoton2P4_Measured = dPhoton2Wrapper->Get_P4_Measured();
		//Step 2
		TLorentzVector locPhoton3P4_Measured = dPhoton3Wrapper->Get_P4_Measured();
		TLorentzVector locPhoton4P4_Measured = dPhoton4Wrapper->Get_P4_Measured();

		/********************************************* GET COMBO RF TIMING INFO *****************************************/

		TLorentzVector locBeamX4_Measured = dComboBeamWrapper->Get_X4_Measured();
		//Double_t locBunchPeriod = dAnalysisUtilities.Get_BeamBunchPeriod(Get_RunNumber());
		// Double_t locDeltaT_RF = dAnalysisUtilities.Get_DeltaT_RF(Get_RunNumber(), locBeamX4_Measured, dComboWrapper);
		// Int_t locRelBeamBucket = dAnalysisUtilities.Get_RelativeBeamBucket(Get_RunNumber(), locBeamX4_Measured, dComboWrapper); // 0 for in-time events, non-zero integer for out-of-time photons
		// Int_t locNumOutOfTimeBunchesInTree = XXX; //YOU need to specify this number
			//Number of out-of-time beam bunches in tree (on a single side, so that total number out-of-time bunches accepted is 2 times this number for left + right bunches) 

		// Bool_t locSkipNearestOutOfTimeBunch = true; // True: skip events from nearest out-of-time bunch on either side (recommended).
		// Int_t locNumOutOfTimeBunchesToUse = locSkipNearestOutOfTimeBunch ? locNumOutOfTimeBunchesInTree-1:locNumOutOfTimeBunchesInTree; 
		// Double_t locAccidentalScalingFactor = dAnalysisUtilities.Get_AccidentalScalingFactor(Get_RunNumber(), locBeamP4.E(), dIsMC); // Ideal value would be 1, but deviations require added factor, which is different for data and MC.
		// Double_t locAccidentalScalingFactorError = dAnalysisUtilities.Get_AccidentalScalingFactorError(Get_RunNumber(), locBeamP4.E()); // Ideal value would be 1, but deviations observed, need added factor.
		// Double_t locHistAccidWeightFactor = locRelBeamBucket==0 ? 1 : -locAccidentalScalingFactor/(2*locNumOutOfTimeBunchesToUse) ; // Weight by 1 for in-time events, ScalingFactor*(1/NBunches) for out-of-time
		// if(locSkipNearestOutOfTimeBunch && abs(locRelBeamBucket)==1) { // Skip nearest out-of-time bunch: tails of in-time distribution also leak in
		// 	dComboWrapper->Set_IsComboCut(true); 
		// 	continue; 
		// } 

		/********************************************* COMBINE FOUR-MOMENTUM ********************************************/

		// DO YOUR STUFF HERE

		// Combine 4-vectors
		TLorentzVector locMissingP4_Measured = locBeamP4_Measured + dTargetP4;
		locMissingP4_Measured -= locProtonP4_Measured + locPhoton1P4_Measured + locPhoton2P4_Measured + locPhoton3P4_Measured + locPhoton4P4_Measured;

		/******************************************** EXECUTE ANALYSIS ACTIONS *******************************************/

		// Loop through the analysis actions, executing them in order for the active particle combo
		dAnalyzeCutActions->Perform_Action(); // Must be executed before Execute_Actions()
		if(!Execute_Actions()) //if the active combo fails a cut, IsComboCutFlag automatically set
			continue;

		//if you manually execute any actions, and it fails a cut, be sure to call:
			//dComboWrapper->Set_IsComboCut(true);

		/**************************************** EXAMPLE: FILL CUSTOM OUTPUT BRANCHES **************************************/

		/*
		TLorentzVector locMyComboP4(8.0, 7.0, 6.0, 5.0);
		//for arrays below: 2nd argument is value, 3rd is array index
		//NOTE: By filling here, AFTER the cuts above, some indices won't be updated (and will be whatever they were from the last event)
			//So, when you draw the branch, be sure to cut on "IsComboCut" to avoid these.
		dTreeInterface->Fill_Fundamental<Float_t>("my_combo_array", -2*loc_i, loc_i);
		dTreeInterface->Fill_TObject<TLorentzVector>("my_p4_array", locMyComboP4, loc_i);
		*/

		/**************************************** EXAMPLE: HISTOGRAM BEAM ENERGY *****************************************/

		//Histogram beam energy (if haven't already)
		if(locUsedSoFar_BeamEnergy.find(locBeamID) == locUsedSoFar_BeamEnergy.end())
		{
			dHist_BeamEnergy->Fill(locBeamP4.E()); // Fills in-time and out-of-time beam photon combos
			//dHist_BeamEnergy->Fill(locBeamP4.E(),locHistAccidWeightFactor); // Alternate version with accidental subtraction

			locUsedSoFar_BeamEnergy.insert(locBeamID);
		}

		/************************************ EXAMPLE: HISTOGRAM MISSING MASS SQUARED ************************************/

		//Missing Mass Squared
		double locMissingMassSquared = locMissingP4_Measured.M2();

		//Uniqueness tracking: Build the map of particles used for the missing mass
			//For beam: Don't want to group with final-state photons. Instead use "Unknown" PID (not ideal, but it's easy).
		map<Particle_t, set<Int_t> > locUsedThisCombo_MissingMass;
		locUsedThisCombo_MissingMass[Unknown].insert(locBeamID); //beam
		locUsedThisCombo_MissingMass[Proton].insert(locProtonTrackID);
		locUsedThisCombo_MissingMass[Gamma].insert(locPhoton1NeutralID);
		locUsedThisCombo_MissingMass[Gamma].insert(locPhoton2NeutralID);
		locUsedThisCombo_MissingMass[Gamma].insert(locPhoton3NeutralID);
		locUsedThisCombo_MissingMass[Gamma].insert(locPhoton4NeutralID);

		//compare to what's been used so far
		if(locUsedSoFar_MissingMass.find(locUsedThisCombo_MissingMass) == locUsedSoFar_MissingMass.end())
		{
			//unique missing mass combo: histogram it, and register this combo of particles
			dHist_MissingMassSquared->Fill(locMissingMassSquared); // Fills in-time and out-of-time beam photon combos
			//dHist_MissingMassSquared->Fill(locMissingMassSquared,locHistAccidWeightFactor); // Alternate version with accidental subtraction

			locUsedSoFar_MissingMass.insert(locUsedThisCombo_MissingMass);
		}

		//E.g. Cut
		//if((locMissingMassSquared < -0.04) || (locMissingMassSquared > 0.04))
		//{
		//	dComboWrapper->Set_IsComboCut(true);
		//	continue;
		//}

		/****************************************** FILL FLAT TREE (IF DESIRED) ******************************************/

		// RECOMMENDED: FILL ACCIDENTAL WEIGHT
		// dFlatTreeInterface->Fill_Fundamental<Double_t>("accidweight",locHistAccidWeightFactor);

		/*
		//FILL ANY CUSTOM BRANCHES FIRST!!
		Int_t locMyInt_Flat = 7;
		dFlatTreeInterface->Fill_Fundamental<Int_t>("flat_my_int", locMyInt_Flat);

		TLorentzVector locMyP4_Flat(4.0, 3.0, 2.0, 1.0);
		dFlatTreeInterface->Fill_TObject<TLorentzVector>("flat_my_p4", locMyP4_Flat);

		for(int loc_j = 0; loc_j < locMyInt_Flat; ++loc_j)
		{
			dFlatTreeInterface->Fill_Fundamental<Int_t>("flat_my_int_array", 3*loc_j, loc_j); //2nd argument = value, 3rd = array index
			TLorentzVector locMyComboP4_Flat(8.0, 7.0, 6.0, 5.0);
			dFlatTreeInterface->Fill_TObject<TLorentzVector>("flat_my_p4_array", locMyComboP4_Flat, loc_j);
		}
		*/

		//FILL FLAT TREE
		//Fill_FlatTree(); //for the active combo

        if (dComboWrapper->Get_ConfidenceLevel_KinFit() < 0.01 ||
            dComboWrapper->Get_NumUnusedShowers() > 0)
		{
           dComboWrapper->Set_IsComboCut(true);
           continue;
		}
		else if (dComboWrapper->Get_ConfidenceLevel_KinFit() > best_comboCL) {
			best_comboCL = dComboWrapper->Get_ConfidenceLevel_KinFit();
			best_combo = loc_i;

            std::vector<TLorentzVector> pp{locBeamP4, dTargetP4, locPhoton1P4, locPhoton2P4,
                                           locPhoton3P4, locPhoton4P4, locProtonP4};
			sqrts = (pp[0] + pp[1]).Mag();
			abst = -(pp[1] - pp[6]).Mag2();
            massEtaPi0 = (pp[2] + pp[3] + pp[4] + pp[5]).Mag();
            massEta = (pp[2] + pp[3]).Mag();
            massPi0 = (pp[4] + pp[5]).Mag();
			phiR = atan2(pp[6](1), pp[6](0));
			for (int i=0; i < 7; ++i)
				pp[i].RotateZ(M_PI - phiR);
			TLorentzVector betaGJ{pp[2] + pp[3] + pp[4] + pp[5]};
            betaGJ *= 1 / betaGJ(3);
			for (int i=0; i < 7; ++i)
				pp[i].Boost(-betaGJ(0), -betaGJ(1), -betaGJ(2));
			double thetaR = atan2(pp[0](0), pp[0](2));
			for (int i=0; i < 7; ++i)
				pp[i].RotateY(thetaR);
			TLorentzVector pEtaGJ{pp[2] + pp[3]};
			thetaGJ = pEtaGJ.Theta();
			phiGJ = pEtaGJ.Phi();
			for (int i=2; i < 4; ++i) {
				pp[i].RotateZ(-phiGJ);
				pp[i].RotateY(thetaGJ);
			}
			TLorentzVector beta12{pp[2] + pp[3]};
			beta12 *= 1 / beta12(3);
			for (int i=2; i < 4; ++i)
				pp[i].Boost(-beta12(0), -beta12(1), -beta12(2));
			for (int i=4; i < 6; ++i) {
				pp[i].RotateZ(-phiGJ);
				pp[i].RotateY(M_PI - thetaGJ);
			}
			TLorentzVector beta34{pp[4] + pp[5]};
			beta34 *= 1 / beta34(3);
			for (int i=4; i < 6; ++i)
				pp[i].Boost(-beta34(0), -beta34(1), -beta34(2));
			thetaEta = pp[2].Theta();
			phiEta = pp[2].Phi();
			thetaPi0 = pp[4].Theta();
			phiPi0 = pp[4].Phi();
		}
	} // end of combo loop

	//FILL HISTOGRAMS: Num combos / events surviving actions
	Fill_NumCombosSurvivedHists();

	/******************************************* LOOP OVER THROWN DATA (OPTIONAL) ***************************************/
/*
	//Thrown beam: just use directly
	if(dThrownBeam != NULL)
		double locEnergy = dThrownBeam->Get_P4().E();

	//Loop over throwns
	for(UInt_t loc_i = 0; loc_i < Get_NumThrown(); ++loc_i)
	{
		//Set branch array indices corresponding to this particle
		dThrownWrapper->Set_ArrayIndex(loc_i);

		//Do stuff with the wrapper here ...
	}
*/
	if (best_comboCL > 0.01) {
		if (Get_NumThrown() == 5) {
            std::vector<TLorentzVector> pp{dThrownBeam->Get_P4(), dTargetP4};
            for (int i=0; i < 5; ++i) {
               dThrownWrapper->Set_ArrayIndex(i);
               pp.push_back(dThrownWrapper->Get_P4());
			}
			sqrts_ = (pp[0] + pp[1]).Mag();
			abst_ = -(pp[1] - pp[6]).Mag2();
            massEtaPi0_ = (pp[2] + pp[3] + pp[4] + pp[5]).Mag();
            massEta_ = (pp[2] + pp[3]).Mag();
            massPi0_ = (pp[4] + pp[5]).Mag();
			phiR_ = atan2(pp[6](1), pp[6](0));
			for (int i=0; i < 7; ++i)
				pp[i].RotateZ(M_PI - phiR_);
			TLorentzVector betaGJ{pp[2] + pp[3] + pp[4] + pp[5]};
            betaGJ *= 1 / betaGJ(3);
			for (int i=0; i < 7; ++i)
				pp[i].Boost(-betaGJ(0), -betaGJ(1), -betaGJ(2));
			double thetaR = atan2(pp[0](0), pp[0](2));
			for (int i=0; i < 7; ++i)
				pp[i].RotateY(thetaR);
			TLorentzVector pEtaGJ{pp[2] + pp[3]};
			thetaGJ_ = pEtaGJ.Theta();
			phiGJ_ = pEtaGJ.Phi();
			for (int i=2; i < 4; ++i) {
				pp[i].RotateZ(-phiGJ_);
				pp[i].RotateY(thetaGJ_);
			}
			TLorentzVector beta12{pp[2] + pp[3]};
			beta12 *= 1 / beta12(3);
			for (int i=2; i < 4; ++i)
				pp[i].Boost(-beta12(0), -beta12(1), -beta12(2));
			for (int i=4; i < 6; ++i) {
				pp[i].RotateZ(-phiGJ_);
				pp[i].RotateY(M_PI - thetaGJ_);
			}
			TLorentzVector beta34{pp[4] + pp[5]};
			beta34 *= 1 / beta34(3);
			for (int i=4; i < 6; ++i)
				pp[i].Boost(-beta34(0), -beta34(1), -beta34(2));
			thetaEta_ = pp[2].Theta();
			phiEta_ = pp[2].Phi();
			thetaPi0_ = pp[4].Theta();
			phiPi0_ = pp[4].Phi();
		}

		runNo = dComboWrapper->Get_RunNumber();
		eventNo = dComboWrapper->Get_EventNumber();
		weight = dComboWrapper->Get_MCWeight();
		dFlatTreeInterface->Fill_Fundamental<ULong64_t>("runNo", runNo);
		dFlatTreeInterface->Fill_Fundamental<ULong64_t>("eventNo", eventNo);
		dFlatTreeInterface->Fill_Fundamental<Double_t>("weight", weight);
		dFlatTreeInterface->Fill_Fundamental<Double_t>("sqrts_", sqrts_);
		dFlatTreeInterface->Fill_Fundamental<Double_t>("abst_", abst_);
		dFlatTreeInterface->Fill_Fundamental<Double_t>("massEtaPi0_", massEtaPi0_);
		dFlatTreeInterface->Fill_Fundamental<Double_t>("massEta_", massEta_);
		dFlatTreeInterface->Fill_Fundamental<Double_t>("massPi0_", massPi0_);
		dFlatTreeInterface->Fill_Fundamental<Double_t>("sqrts", sqrts);
		dFlatTreeInterface->Fill_Fundamental<Double_t>("abst", abst);
		dFlatTreeInterface->Fill_Fundamental<Double_t>("massEtaPi0", massEtaPi0);
		dFlatTreeInterface->Fill_Fundamental<Double_t>("massEta", massEta);
		dFlatTreeInterface->Fill_Fundamental<Double_t>("massPi0", massPi0);
		dFlatTreeInterface->Fill_Fundamental<Double_t>("phiR_", phiR_);
		dFlatTreeInterface->Fill_Fundamental<Double_t>("thetaGJ_", thetaGJ_);
		dFlatTreeInterface->Fill_Fundamental<Double_t>("phiGJ_", phiGJ_);
		dFlatTreeInterface->Fill_Fundamental<Double_t>("thetaEta_", thetaEta_);
		dFlatTreeInterface->Fill_Fundamental<Double_t>("phiEta_", phiEta_);
		dFlatTreeInterface->Fill_Fundamental<Double_t>("thetaPi0_", thetaPi0_);
		dFlatTreeInterface->Fill_Fundamental<Double_t>("phiPi0_", phiPi0_);
		dFlatTreeInterface->Fill_Fundamental<Double_t>("thetaGJ", thetaGJ);
		dFlatTreeInterface->Fill_Fundamental<Double_t>("phiR", phiR);
		dFlatTreeInterface->Fill_Fundamental<Double_t>("phiGJ", phiGJ);
		dFlatTreeInterface->Fill_Fundamental<Double_t>("thetaEta", thetaEta);
		dFlatTreeInterface->Fill_Fundamental<Double_t>("phiEta", phiEta);
		dFlatTreeInterface->Fill_Fundamental<Double_t>("thetaPi0", thetaPi0);
		dFlatTreeInterface->Fill_Fundamental<Double_t>("phiPi0", phiPi0);
		momentsGJ = 0;
		for (int L=0; L < 13; L += 1) {
			for (int M=-L; M <= L; M += 1) {
				YmomGJ[momentsGJ] = angular_moment(L, M, thetaGJ, phiGJ);
				YmomGJ_[momentsGJ] = angular_moment(L, M, thetaGJ_, phiGJ_);

				model1_use_generated(0);
				model1moment[momentsGJ] = model1_moment(L, M);
				model1_use_generated(1);
				model1moment_[momentsGJ] = model1_moment(L, M);

				++momentsGJ;
			}
		}
		dFlatTreeInterface->Fill_Fundamental<UInt_t>("momentsGJ", momentsGJ);
        for (unsigned int m=0; m < momentsGJ; ++m) {
			dFlatTreeInterface->Fill_Fundamental<Double_t>("YmomGJ", YmomGJ[m], m);
			dFlatTreeInterface->Fill_Fundamental<Double_t>("YmomGJ_", YmomGJ_[m], m);

			dFlatTreeInterface->Fill_Fundamental<Double_t>("model1moment", model1moment[m], m);
			dFlatTreeInterface->Fill_Fundamental<Double_t>("model1moment_", model1moment_[m], m);
		}
        momentsEta = 0;
		for (int L=0; L < 13; L += 2) {
			for (int M=-L; M <= L; M += 1) {
				YmomEta[momentsEta] = angular_moment(L, M, thetaEta, phiEta);
				YmomEta_[momentsEta] = angular_moment(L, M, thetaEta_, phiEta_);
				++momentsEta;
			}
		}
		dFlatTreeInterface->Fill_Fundamental<UInt_t>("momentsEta", momentsEta);
        for (unsigned int m=0; m < momentsGJ; ++m) {
			dFlatTreeInterface->Fill_Fundamental<Double_t>("YmomEta", YmomEta[m], m);
			dFlatTreeInterface->Fill_Fundamental<Double_t>("YmomEta_", YmomEta_[m], m);
		}
        momentsPi0 = 0;
		for (int L=0; L < 13; L += 2) {
			for (int M=-L; M <= L; M += 1) {
				YmomPi0[momentsPi0] = angular_moment(L, M, thetaPi0, phiPi0);
				YmomPi0_[momentsPi0] = angular_moment(L, M, thetaPi0_, phiPi0_);
				++momentsPi0;
			}
		}
		dFlatTreeInterface->Fill_Fundamental<UInt_t>("momentsPi0", momentsPi0);
        for (unsigned int m=0; m < momentsGJ; ++m) {
			dFlatTreeInterface->Fill_Fundamental<Double_t>("YmomPi0", YmomPi0[m], m);
			dFlatTreeInterface->Fill_Fundamental<Double_t>("YmomPi0_", YmomPi0_[m], m);
		}
		Fill_FlatTree();
	}

	/****************************************** LOOP OVER OTHER ARRAYS (OPTIONAL) ***************************************/
/*
	//Loop over beam particles (note, only those appearing in combos are present)
	for(UInt_t loc_i = 0; loc_i < Get_NumBeam(); ++loc_i)
	{
		//Set branch array indices corresponding to this particle
		dBeamWrapper->Set_ArrayIndex(loc_i);

		//Do stuff with the wrapper here ...
	}

	//Loop over charged track hypotheses
	for(UInt_t loc_i = 0; loc_i < Get_NumChargedHypos(); ++loc_i)
	{
		//Set branch array indices corresponding to this particle
		dChargedHypoWrapper->Set_ArrayIndex(loc_i);

		//Do stuff with the wrapper here ...
	}

	//Loop over neutral particle hypotheses
	for(UInt_t loc_i = 0; loc_i < Get_NumNeutralHypos(); ++loc_i)
	{
		//Set branch array indices corresponding to this particle
		dNeutralHypoWrapper->Set_ArrayIndex(loc_i);

		//Do stuff with the wrapper here ...
	}
*/

	/************************************ EXAMPLE: FILL CLONE OF TTREE HERE WITH CUTS APPLIED ************************************/
/*
	Bool_t locIsEventCut = true;
	for(UInt_t loc_i = 0; loc_i < Get_NumCombos(); ++loc_i) {
		//Set branch array indices for combo and all combo particles
		dComboWrapper->Set_ComboIndex(loc_i);
		// Is used to indicate when combos have been cut
		if(dComboWrapper->Get_IsComboCut())
			continue;
		locIsEventCut = false; // At least one combo succeeded
		break;
	}
	if(!locIsEventCut && dOutputTreeFileName != "")
		Fill_OutputTree();
*/

	return kTRUE;
}

void DSelector_etapi0_moments::Finalize(void)
{
	//Save anything to output here that you do not want to be in the default DSelector output ROOT file.

	//Otherwise, don't do anything else (especially if you are using PROOF).
		//If you are using PROOF, this function is called on each thread,
		//so anything you do will not have the combined information from the various threads.
		//Besides, it is best-practice to do post-processing (e.g. fitting) separately, in case there is a problem.

	//DO YOUR STUFF HERE

	//CALL THIS LAST
	DSelector::Finalize(); //Saves results to the output file
}

double DSelector_etapi0_moments::angular_moment(int L, int M, double theta, double phi) const
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

double DSelector_etapi0_moments::model1_moment(int L, int M) const
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

double DSelector_etapi0_moments::model1_Lmax() const
{
   return 3;
}

std::vector<double> DSelector_etapi0_moments::model1_amplitude(int L, int M) const
{
   // Returns the complex amplitude specified in model 1 for partial wave L,M
   // in the Gottfried Jackson frame, represented as a two-component vector with
   // real part first, imaginary part second.

#define SQR(x) ((x)*(x))

   double mX(massEtaPi0);
   double t(abst);
   if (model1_generated_kinematics) {
      mX = massEtaPi0_;
      t = abst_;
   }

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

double DSelector_etapi0_moments::model1_density(int source) const
{
   // Returns the model 1 density at Gottfried-Jackson angles thetaGJ,phiGJ
   // computed either from the complex amplitude sum (source==0) or from
   // a sum over computed moments (source==1).

   double theta(thetaGJ);
   double phi(phiGJ);
   if (model1_generated_kinematics) {
      theta = thetaGJ_;
      phi = phiGJ_;
   }

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

void DSelector_etapi0_moments::model1_use_generated(int yesno)
{
   model1_generated_kinematics = yesno;
}
