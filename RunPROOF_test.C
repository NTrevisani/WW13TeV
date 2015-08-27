///////////////////////////////////////////////////////////////////////
//
//    FILE: RunProof.C
// AUTHORS: I. Gonzalez Caballero, A.Y. Rodriguez Marrero
//    DATE: 2010
//
// CONTENT: Main macro to run over MiniTrees or TESCO Trees using PROOF
//          in PROOF-Lite, PROOF-Cluster or Sequential mode
///////////////////////////////////////////////////////////////////////

#include<iostream>
#include<string>

TProof* proof = 0;

/////////////////////////////////////////////////////////////
/////////////           TOP ANALYSIS            /////////////
/////////////////////////////////////////////////////////////

void RunPROOF_test(double luminosity,
		   const char* data, 
		   TString  theSample,  
		   Int_t JetChannel,
		   TString FlavorChannel,
		   TString proofMode_,
		   TString SameSign,
		   TString MuonID) 
{
 
  // This loads all the PROOF Analysis Framework utilities
  gROOT->LoadMacro("$PAFPATH/PAF.C");
  
  //  double luminosity = 19468;
  TString dataInfo;
  TString mcInfo;
  
  Bool_t runAtOviedo = false;
  Bool_t runAtIfca   = !runAtOviedo;
  
  /////////////////////////////////////////////////////////////////////////
  //
  // PROOF SETTINGS
  // ==============
  //
  // (You may inspect scripts/PAFOptions.h to see all the posible settings)
  //
  // Edit the lines below to select tree type, input data sample, output
  // file name ...
  //
  
  ///////////////////////////////
  // PROOF MODE
  // Defines the mode in which you want to run PROOF. Read the documentation
  // for details:
  // * kSequential: No PROOF. Plain sequential code
  // * kLite: PROOF Lite mode
  // * kCluster: PROOF Cluter mode
  // * kPoD: PROOF on Demand mode
  if (proofMode_ == "Cluster" || proofMode_ == "kCluster" || proofMode_ == "cluster" || proofMode_ == "kcluster")
    gPAFOptions->proofMode = kCluster;
  else if (proofMode_ == "Sequential" || proofMode_ == "kSequential" || proofMode_ == "sequential" || proofMode_ == "ksequential")
    gPAFOptions->proofMode = kSequential;
  else if (proofMode_ == "Lite" || proofMode_ == "kLite" || proofMode_ == "lite" || proofMode_ == "klite")
    gPAFOptions->proofMode = kLite;
  else{
    cout<<"Please select a valid PROOF operating mode: Cluster, Sequential, or Lite"<<endl;
    return;
  }
  //
  // Optional parameters for PROOF Cluster (kCluster):
  //   + The number of slots you would like to use (default is 10)
  //gPAFOptions->NSlots = 10; 
  //   + Proof server and port (default are proof.ifca.es and 1093) 
   gPAFOptions->proofServer = "proof.ifca.es";
   gPAFOptions->proofServerPort = 1093;
  //   + Maximum number of slaves per node (default is 9999, i.e. all)
  // gPAFOptions->maxSlavesPerNode = 9999;
  //
  // Start PROOF
  //
  cout << ">> Starting PROOF..." << endl;
  proof = InitProof(); 
  if (!proof && gPAFOptions->proofMode != kSequential) {
    cerr << "ERROR: I could not initialise a PROOF session!" << endl;
    return;
  }
 
  if (SameSign != "SS" && SameSign != "OS"){
    cout<<"Please choose if you want to use same sign (SS) or opposite sign (OS) selections"<<endl;
    return;
  }

  ///////////////////////////////
  // TREE TYPE
  // Defines the data formats that may be used: MiniTrees (default) or TESCO
  gPAFOptions->SetTreeType(kMiniTrees);
  //gPAFOptions->treeType = kTESCO;


  //  gPAFOptions->SetTreeDir("demo");
  gPAFOptions->SetTreeName("latino");   

  ///////////////////////////////
  // INPUT DATA SAMPLE
  //
  //TString dataPath="/gpfs/csic_projects/cms/arodrig/"; //IFCA   (gridui)
  // 1) Set the base path to files
  // TString dataPath="/gpfs/csic_projects/cms/calderon/TreesCSA14"; //IFCA   (gridui)
  //  TString dataPath="/gpfs/csic_projects/tier3data/LatinosSkims/ReducedTrees/DiferentialXSection";
  // TString dataPath="/gpfs/csic_projects/tier3data"; //IFCA   (gridui)
  //        TString dataPath="/hadoop";                      //UniOvi (fanaeui)
  //   TString dataPath="/pool/data1/MiniTrees";        //CERN   (cmsovd02)

  TString filesPath = "/gpfs/csic_projects/tier3data/LatinosSkims/ReducedTrees/Systematics2013_nominals_fromMaiko/";
  
  // 2) Asign the files to gPAFOptions->dataFiles (a vector<TString>)
  // Ex. MiniTree...
  //  gPAFOptions->dataFiles.push_back(dataPath + "/Data/Data7TeVRun2010A/Tree_Mu_Run2010A_Sep27ReReco_Skim2LPt1010.root");
  //
  // Note: You may consider using DatasetManager for this. See the
  //       documentation in the wiki

  TString Signal = data;
  
  bool isdata = true;
  int nEventsInTheSample = 99555; 
  double xSection =  1.0 ;
  int whichRun = 2; 
  

  // *********** DATA
  
  if (Signal=="test") {
    
    //50ns
    //**********************************************************************************************************************

    if (theSample == "Data2015"){
      gPAFOptions->dataFiles.push_back("/gpfs/csic_projects/cms/trevisanin/latinoNoElement/Data/latino_DoubleEG.root");
      gPAFOptions->dataFiles.push_back("/gpfs/csic_projects/cms/trevisanin/latinoNoElement/Data/latino_MuonEG.root");
      gPAFOptions->dataFiles.push_back("/gpfs/csic_projects/cms/trevisanin/latinoNoElement/Data/latino_SingleElectron.root");
      gPAFOptions->dataFiles.push_back("/gpfs/csic_projects/cms/trevisanin/latinoNoElement/Data/latino_SingleMuon.root");
      //gPAFOptions->dataFiles.push_back("/gpfs/csic_projects/cms/trevisanin/latinoNoElement/Data/latino_SingleMu.root");
    }

    else if (theSample == "WW50"){
      gPAFOptions->dataFiles.push_back("/gpfs/csic_projects/cms/trevisanin/latinoNoElement/50ns/latino_WWTo2L2Nu.root");
    }

    else if (theSample == "WJets50"){
      gPAFOptions->dataFiles.push_back("/gpfs/csic_projects/cms/trevisanin/latinoNoElement/50ns/latino_WJetsToLNu.root");
    }

    else if (theSample == "VV50"){
      gPAFOptions->dataFiles.push_back("/gpfs/csic_projects/cms/trevisanin/latinoNoElement/50ns/latino_WZ.root");
      gPAFOptions->dataFiles.push_back("/gpfs/csic_projects/cms/trevisanin/latinoNoElement/50ns/latino_ZZ.root");
    }

    else if (theSample == "WZ50"){
      gPAFOptions->dataFiles.push_back("/gpfs/csic_projects/cms/trevisanin/latinoNoElement/50ns/latino_WZ.root");
    }

    else if (theSample == "ZZ50"){
      gPAFOptions->dataFiles.push_back("/gpfs/csic_projects/cms/trevisanin/latinoNoElement/50ns/latino_ZZ.root");
    }

    else if (theSample == "TTJets50"){
      gPAFOptions->dataFiles.push_back("/gpfs/csic_projects/cms/trevisanin/latinoNoElement/50ns/latino_TTJets.root");
    }
 
    else if (theSample == "TTbar50"){
      gPAFOptions->dataFiles.push_back("/gpfs/csic_projects/cms/trevisanin/newLatino/50ns/latino_TTTo2L2Nu.root");
    }

    else if (theSample == "DY50"){
      gPAFOptions->dataFiles.push_back("/gpfs/csic_projects/cms/trevisanin/latinoNoElement/50ns/latino_DYJetsToLL.root");
    }

    //25ns
    //**********************************************************************************************************************

    else if (theSample == "WW"){
      gPAFOptions->dataFiles.push_back("/gpfs/csic_projects/cms/trevisanin/latinoNoElement/25ns/latino_WWTo2L2Nu.root");
    }

    else if (theSample == "WJets"){
      gPAFOptions->dataFiles.push_back("/gpfs/csic_projects/cms/trevisanin/latinoNoElement/25ns/latino_WJetsToLNu.root");
    }

    else if (theSample == "VV"){
      gPAFOptions->dataFiles.push_back("/gpfs/csic_projects/cms/trevisanin/latinoNoElement/25ns/latino_WZTo3LNu.root");
      gPAFOptions->dataFiles.push_back("/gpfs/csic_projects/cms/trevisanin/latinoNoElement/25ns/latino_ZZ4l.root");
    }

    else if (theSample == "TTJets"){
      gPAFOptions->dataFiles.push_back("/gpfs/csic_projects/cms/trevisanin/latinoNoElement/25ns/latino_TTJets.root");
    }
 
    else if (theSample == "SingleTop"){
      gPAFOptions->dataFiles.push_back("/gpfs/csic_projects/cms/trevisanin/latinoNoElement/25ns/latino_ST_t-channel_antitop.root");
      gPAFOptions->dataFiles.push_back("/gpfs/csic_projects/cms/trevisanin/latinoNoElement/25ns/latino_ST_t-channel_top.root");
      gPAFOptions->dataFiles.push_back("/gpfs/csic_projects/cms/trevisanin/latinoNoElement/25ns/latino_ST_tW_antitop.root");
      gPAFOptions->dataFiles.push_back("/gpfs/csic_projects/cms/trevisanin/latinoNoElement/25ns/latino_ST_tW_top.root");
    }

    else if (theSample == "Top"){
      gPAFOptions->dataFiles.push_back("/gpfs/csic_projects/cms/trevisanin/latinoNoElement/25ns/latino_ST_t-channel_antitop.root");
      gPAFOptions->dataFiles.push_back("/gpfs/csic_projects/cms/trevisanin/latinoNoElement/25ns/latino_ST_t-channel_top.root");
      gPAFOptions->dataFiles.push_back("/gpfs/csic_projects/cms/trevisanin/latinoNoElement/25ns/latino_ST_tW_antitop.root");
      gPAFOptions->dataFiles.push_back("/gpfs/csic_projects/cms/trevisanin/latinoNoElement/25ns/latino_ST_tW_top.root");
      gPAFOptions->dataFiles.push_back("/gpfs/csic_projects/cms/trevisanin/latinoNoElement/25ns/latino_TTJets.root");
    }

    else if (theSample == "DY"){
      gPAFOptions->dataFiles.push_back("/gpfs/csic_projects/cms/trevisanin/latinoNoElement/25ns/latino_DYJetsToLL_M-10to50.root");
    }
    else if (theSample == "HWW"){
      gPAFOptions->dataFiles.push_back("/gpfs/csic_projects/cms/trevisanin/latinoNoElement/25ns/latino_ggHWW120.root");
    }

    //**********************************************************************************************************************

    else if (theSample == "VBF"){
      gPAFOptions->dataFiles.push_back("/gpfs/csic_projects/cms/trevisanin/newLatino/latino_VBF125.root");
    }

    else if (theSample == "QCD"){
      gPAFOptions->dataFiles.push_back("/gpfs/csic_projects/cms/piedra/latino/RunII/MC_Spring15/25ns/latino_QCD15to20EM.root");
    }

    else if (theSample == "WJets"){
      gPAFOptions->dataFiles.push_back("/gpfs/csic_projects/cms/trevisanin/newLatino/25ns/latino_WJetsToLNu.root");
    }

    else if (theSample == "ZZ"){
      gPAFOptions->dataFiles.push_back("/gpfs/csic_projects/cms/trevisanin/newLatino/25ns/latino_ZZ.root");
    }

    else if (theSample == "TTJets"){
      gPAFOptions->dataFiles.push_back("/gpfs/csic_projects/cms/trevisanin/newLatino/latino_TTJets.root");
    }

    else if (theSample == "WJets8TeV"){
      gPAFOptions->dataFiles.push_back("/gpfs/csic_projects/tier3data/LatinosSkims/ReducedTrees/R53X_S1_V08_S2_V09_S3_V13/MC_LooseLoose/4L/latino_080_WJetsToLNuMad.root");
      }    
  }
    else
      return;
    

///////////////////////////////
// OUTPUT FILE NAME
// Specify the name of the file where you want your histograms to be saved


std::ostringstream out;

TString outTest = out.str();

//TString output = TString("csa14_GluGluToHToWWTo2LAndTau2Nu_M-125_13TeV_PUS14_Vertex.root"); 
//TString output = TString("csa14_W1234JetsToLNu_Tune4C_13TeV_PUS14.root");
//TString output = TString("csa14_WToMuNu_Tune4C_13TeV_PUS14.root");
//TString output = TString("csa14_GluGluToHToWWTo2LAndTau2Nu_M-125_13TeV_PU20bx25.root");
//TString output = TString("csa14_WJetsToLNu_13TeV-madgraph_PU20bx25.root"); 
//TString output = TString("csa14_WJets_Madgraph_8Tev.root");
//TString output = TString("csa14_HWW125_8Tev.root");

TString jets = Form("%d",JetChannel);

//  TString path = Form("rootfiles/%djet/%s/", JetChannel, FlavorChannel.Data());

TString output = TString(theSample+".root");



//TString output = TString(Signal+".root"); 

gPAFOptions->outputFile=output;



///////////////////////////////
// PARAMETERS FOR THE ANALYSIS
// This parameters are passed to the analysis class and can be use there.
// They are stored in a InputParameters object. They are saved to the 
// output file.
// See packages/InputParameters/InputParameters.h for information on how
// to use this class.

//  std::string sample = theSample;

gPAFOptions->inputParameters = new InputParameters();

gPAFOptions->inputParameters->SetNamedBool("IsDATA", isdata);
gPAFOptions->inputParameters->SetNamedString("Signal", data);
gPAFOptions->inputParameters->SetNamedDouble("XSection", xSection);
gPAFOptions->inputParameters->SetNamedDouble("Luminosity", luminosity);
gPAFOptions->inputParameters->SetNamedInt("NEvents", nEventsInTheSample); // all
gPAFOptions->inputParameters->SetNamedFloat("luminosityPU", 19468.3);  
gPAFOptions->inputParameters->SetNamedInt("WhichRun", whichRun);
gPAFOptions->inputParameters->SetNamedString("theSample", theSample.Data());
gPAFOptions->inputParameters->SetNamedString("FlavorChannel", FlavorChannel.Data());
gPAFOptions->inputParameters->SetNamedString("SameSign", SameSign.Data());
gPAFOptions->inputParameters->SetNamedString("MuonID", MuonID.Data());
gPAFOptions->inputParameters->SetNamedInt("jetChannel", JetChannel);

////// I.G.
//Find the total number of entries in the dataset and send it to the input parameters
/*  TChain* chain = new TChain("Tree", "Tree");
    for (unsigned int i = 0; i < dataFiles.size(); i++) 
    chain->Add(dataFiles[i]);
    gPAFOptions->inputParameters->SetNamedInt("NEventsTotal", chain->GetEntries()); //after skimming
    TString eventsfile(gSystem->pwd());
    eventsfile+="/";
    eventsfile+=Signal;
    eventsfile+="_events.log";
    gPAFOptions->inputParameters->SetNamedString("fFileList", (const char*) eventsfile);
*/

///////////////////////////////
// DYNAMIC HISTOGRAMS
// Specify the name of the histograms you would like to monitor as they are
// filled by PROOF
//
//  gPAFOptions->dynamicHistograms.push_back("myHistogram");
//...

/////////////////////////////
// NUMBER OF EVENTS 
// Specify the number (Long64_t) of events to process.
// Set it to -1 to use the full sample.

gPAFOptions->SetNEvents(-1);

//
/////////////////////////////////////////////////////////////////////////







/////////////////////////////////////////////////////////////////////////
//
// EXTRA proof settings:
// ====================
//
// It is unlikely that you need to edit anything below. At least at the
// beginning of your PAF experience. However we provide a couple of hooks
// for extensions.
//

///////////////////////////////
// NAME OF ANALYSIS CLASS. 
// If 0 the default name schema will be used, i.e. depending on the value
// of gPAFOptions->treeType: MyAnalysisTESCO or MyAnalsyisMiniTrees
//

gPAFOptions->SetAnalysis("test");


///////////////////////////////
// ADDITIONAL PACKAGES TO BE UPLOADED TO PROOF.
// The mandatory ones are added automatically in PAFOptions
//


gPAFOptions->AddPackage("PUWeight");
//  gPAFOptions->AddPackage("MuonIsoMVA");


///////////////////////////////
// CONTROL OUTPUT AND CHECKS
// + If true (default) PAF checks for new version in CVS every time
// gPAFOptions->checkVersion = true;
// + If true (default) the output file is reopened so the objects in the
//   file can be interactively accessed. The object in the output are also
//   listed
// gPAFOptions->reopenOutputFile = false;

//
/////////////////////////////////////////////////////////////////////////




/////////////////////////////////////////////////////////////////////////
//
// RUN THE ANALYSIS
// ================
//
// If something needs to be edited below (or inside), contact the 
// developers.
//
// Run the analysis
//
//gPAFOptions->reopenOutputFile = false;
gPAFOptions->reopenOutputFileRemoved= false;

if (!RunAnalysis())
  cerr << "ERROR: There was a problem running the analysis!" << endl;
//
/////////////////////////////////////////////////////////////////////////


}
