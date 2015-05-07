#include <TFile.h>
#include "MuonIsoMVA.h"
#include "TMVA/Tools.h"
#include "TMVA/Reader.h"
#include <math.h>

//#define DEBUG
using namespace std;

//--------------------------------------------------------------------------------------------------
MuonIsoMVA::MuonIsoMVA() :
  fMethodname("BDT method"),
  fIsInitialized(kFALSE){
  
  cout << "Creating Constructor" << endl;
  // Constructor.
  fTMVAReader = vector<TMVA::Reader*>(0);
  cout << "Constructed" << endl;
  //  for(UInt_t i=0; i<4; ++i) {    fTMVAReader[i] = 0;  }
}


//--------------------------------------------------------------------------------------------------
MuonIsoMVA::~MuonIsoMVA(){
  for(UInt_t i=0; i<fTMVAReader.size(); ++i) {
    if (fTMVAReader[i]) delete fTMVAReader[i];
  }
}

//--------------------------------------------------------------------------------------------------
void MuonIsoMVA::Initialize( string methodName,
			     MVAType type, 
			     string Pt10To20BarrelWeights,
			     string Pt10To20EndcapWeights,
			     string Pt20ToInfBarrelWeights,
			     string Pt20ToInfEndcapWeights			     ){
  
#ifdef DEBUG
  cout << "Initialise( " 
       << methodName             << ", " << type                    << ", " 
       << Pt10To20BarrelWeights  << ", " << Pt10To20BarrelWeights   << ", "
       << Pt20ToInfBarrelWeights << ", " <<  Pt20ToInfEndcapWeights << " )" <<endl;
#endif
  
  for (unsigned int i=0; i<fTMVAReader.size(); ++i){
    if (fTMVAReader[i])  delete fTMVAReader[i];
  }
  fTMVAReader.clear();
  
  fIsInitialized = kTRUE;
  fMethodname = methodName;
  
  unsigned int ExpectedNBins = 0;
  if (type == kIsoRings) ExpectedNBins = 6;
  else                   ExpectedNBins = 4;
  
  for(UInt_t i=0; i<ExpectedNBins; ++i) {
    TMVA::Reader *tmpTMVAReader = new TMVA::Reader( "!Color:!Silent:Error" );  
    tmpTMVAReader->SetVerbose(kTRUE);
    
    if (type == kIsoDeltaRV1){
      tmpTMVAReader->AddVariable("PFChargedDR04",     &fMVAVar_MuRelIsoPFChargedDR04  );
      tmpTMVAReader->AddVariable("PFNePhCorrDR04",    &fMVAVar_MuRelIsoPFNePhCorrDR04 );
      tmpTMVAReader->AddVariable("SumDeltaRCharged",  &fMVAVar_MuDeltaRSumCharged     );
      tmpTMVAReader->AddVariable("DeltaRMeanCharged", &fMVAVar_MuDeltaRMeanCharged    );
    }
    if (type == kIsoDeltaR){
      tmpTMVAReader->AddVariable("PFCharged",  &fMVAVar_MuRelIsoPFCharged );
      tmpTMVAReader->AddVariable("PFNeutral",  &fMVAVar_MuRelIsoPFNeutral );
      tmpTMVAReader->AddVariable("PFPhotons",  &fMVAVar_MuRelIsoPFPhotons );
      tmpTMVAReader->AddVariable("SumDeltaR",  &fMVAVar_MuDeltaRSum       );
      tmpTMVAReader->AddVariable("DeltaRMean", &fMVAVar_MuDeltaRMean      );
      tmpTMVAReader->AddVariable("Density",    &fMVAVar_MuDensity         );
    }
    if (type == kIsoRings){
      tmpTMVAReader->AddVariable("PFChargedDR0p0to0p1", &fMVAVar_ChargedIso_DR0p0To0p1); 
      tmpTMVAReader->AddVariable("PFChargedDR0p1to0p2", &fMVAVar_ChargedIso_DR0p1To0p2);
      tmpTMVAReader->AddVariable("PFChargedDR0p2to0p3", &fMVAVar_ChargedIso_DR0p2To0p3);
      tmpTMVAReader->AddVariable("PFChargedDR0p3to0p4", &fMVAVar_ChargedIso_DR0p3To0p4);
      tmpTMVAReader->AddVariable("PFChargedDR0p4to0p5", &fMVAVar_ChargedIso_DR0p4To0p5);
      
      tmpTMVAReader->AddVariable("PFNeutralDR0p0to0p1", &fMVAVar_NeutralHadronIso_DR0p0To0p1);
      tmpTMVAReader->AddVariable("PFNeutralDR0p1to0p2", &fMVAVar_NeutralHadronIso_DR0p1To0p2);
      tmpTMVAReader->AddVariable("PFNeutralDR0p2to0p3", &fMVAVar_NeutralHadronIso_DR0p2To0p3);
      tmpTMVAReader->AddVariable("PFNeutralDR0p3to0p4", &fMVAVar_NeutralHadronIso_DR0p3To0p4);
      tmpTMVAReader->AddVariable("PFNeutralDR0p4to0p5", &fMVAVar_NeutralHadronIso_DR0p4To0p5);
     
      tmpTMVAReader->AddVariable("PFPhotonsDR0p0to0p1", &fMVAVar_GammaIso_DR0p0To0p1);
      tmpTMVAReader->AddVariable("PFPhotonsDR0p1to0p2", &fMVAVar_GammaIso_DR0p1To0p2);
      tmpTMVAReader->AddVariable("PFPhotonsDR0p2to0p3", &fMVAVar_GammaIso_DR0p2To0p3);
      tmpTMVAReader->AddVariable("PFPhotonsDR0p3to0p4", &fMVAVar_GammaIso_DR0p3To0p4);
      tmpTMVAReader->AddVariable("PFPhotonsDR0p4to0p5", &fMVAVar_GammaIso_DR0p4To0p5);
     
      tmpTMVAReader->AddSpectator("Pt",        &fMVAVar_Pt);
      tmpTMVAReader->AddSpectator("Eta",       &fMVAVar_Eta);
    }
    if (type == kIsoCombined){
      tmpTMVAReader->AddVariable("PFCharged",  &fMVAVar_MuRelIsoPFCharged );
      tmpTMVAReader->AddVariable("PFNeutral",  &fMVAVar_MuRelIsoPFNeutral );
      tmpTMVAReader->AddVariable("PFPhotons",  &fMVAVar_MuRelIsoPFPhotons );
      
      tmpTMVAReader->AddVariable("PFChargedDR0p0to0p1", &fMVAVar_ChargedIso_DR0p0To0p1); 
      tmpTMVAReader->AddVariable("PFChargedDR0p1to0p2", &fMVAVar_ChargedIso_DR0p1To0p2);
      tmpTMVAReader->AddVariable("PFChargedDR0p2to0p3", &fMVAVar_ChargedIso_DR0p2To0p3);
      tmpTMVAReader->AddVariable("PFChargedDR0p3to0p4", &fMVAVar_ChargedIso_DR0p3To0p4);
      tmpTMVAReader->AddVariable("PFChargedDR0p4to0p5", &fMVAVar_ChargedIso_DR0p4To0p5);
      
      tmpTMVAReader->AddVariable("PFNeutralDR0p0to0p1", &fMVAVar_NeutralHadronIso_DR0p0To0p1);
      tmpTMVAReader->AddVariable("PFNeutralDR0p1to0p2", &fMVAVar_NeutralHadronIso_DR0p1To0p2);
      tmpTMVAReader->AddVariable("PFNeutralDR0p2to0p3", &fMVAVar_NeutralHadronIso_DR0p2To0p3);
      tmpTMVAReader->AddVariable("PFNeutralDR0p3to0p4", &fMVAVar_NeutralHadronIso_DR0p3To0p4);
      tmpTMVAReader->AddVariable("PFNeutralDR0p4to0p5", &fMVAVar_NeutralHadronIso_DR0p4To0p5);
     
      tmpTMVAReader->AddVariable("PFPhotonsDR0p0to0p1", &fMVAVar_GammaIso_DR0p0To0p1);
      tmpTMVAReader->AddVariable("PFPhotonsDR0p1to0p2", &fMVAVar_GammaIso_DR0p1To0p2);
      tmpTMVAReader->AddVariable("PFPhotonsDR0p2to0p3", &fMVAVar_GammaIso_DR0p2To0p3);
      tmpTMVAReader->AddVariable("PFPhotonsDR0p3to0p4", &fMVAVar_GammaIso_DR0p3To0p4);
      tmpTMVAReader->AddVariable("PFPhotonsDR0p4to0p5", &fMVAVar_GammaIso_DR0p4To0p5);
     
      tmpTMVAReader->AddVariable("SumDeltaR",  &fMVAVar_MuDeltaRSum       );
      tmpTMVAReader->AddVariable("DeltaRMean", &fMVAVar_MuDeltaRMean      );
      tmpTMVAReader->AddVariable("Density",    &fMVAVar_MuDensity         );
      
      tmpTMVAReader->AddSpectator("Pt",        &fMVAVar_Pt);
      tmpTMVAReader->AddSpectator("Eta",       &fMVAVar_Eta);
    }
    if (i==0) tmpTMVAReader->BookMVA(fMethodname , Pt10To20BarrelWeights  );
    if (i==1) tmpTMVAReader->BookMVA(fMethodname , Pt10To20EndcapWeights  );
    if (i==2) tmpTMVAReader->BookMVA(fMethodname , Pt20ToInfBarrelWeights );
    if (i==3) tmpTMVAReader->BookMVA(fMethodname , Pt20ToInfEndcapWeights );
    
    fTMVAReader.push_back(tmpTMVAReader);
  }

#ifdef DEBUG
  cout << "Muon ID MVA Initialization\n";
  cout << "MethodName : " << fMethodname << endl;
  cout << "Load weights file : " << Pt10To20BarrelWeights  << endl;
  cout << "Load weights file : " << Pt10To20EndcapWeights  << endl;
  cout << "Load weights file : " << Pt20ToInfBarrelWeights << endl;
  cout << "Load weights file : " << Pt20ToInfEndcapWeights << endl;
#endif
}

void MuonIsoMVA::Initialize( string methodName,
			     MVAType type, 
			     std::vector<std::string> weightfiles){
  
  for (unsigned int i=0; i<fTMVAReader.size(); ++i){
    if (fTMVAReader[i])  delete fTMVAReader[i];
  }
  fTMVAReader.clear();
  fIsInitialized = kTRUE;
  fMethodname = methodName;
  
  unsigned int ExpectedNBins = 0;
  if (type == kIsoRings) ExpectedNBins = 6;
  else                   ExpectedNBins = 4;
  
  for(UInt_t i=0; i<ExpectedNBins; ++i) {
    //    if (fTMVAReader[i]) delete fTMVAReader[i];
    
    TMVA::Reader *tmpTMVAReader = new TMVA::Reader( "!Color:!Silent:Error" );  
    tmpTMVAReader->SetVerbose(kTRUE);
    
    if (type == kIsoDeltaR){
      tmpTMVAReader->AddVariable("PFCharged",  &fMVAVar_MuRelIsoPFCharged );
      tmpTMVAReader->AddVariable("PFNeutral",  &fMVAVar_MuRelIsoPFNeutral );
      tmpTMVAReader->AddVariable("PFPhotons",  &fMVAVar_MuRelIsoPFPhotons );
      tmpTMVAReader->AddVariable("SumDeltaR",  &fMVAVar_MuDeltaRSum       );
      tmpTMVAReader->AddVariable("DeltaRMean", &fMVAVar_MuDeltaRMean      );
      tmpTMVAReader->AddVariable("Density",    &fMVAVar_MuDensity         );
    }
    if (type == kIsoDeltaRV1){
      tmpTMVAReader->AddVariable("PFChargedDR04",     &fMVAVar_MuRelIsoPFChargedDR04  );
      tmpTMVAReader->AddVariable("PFNePhCorrDR04",    &fMVAVar_MuRelIsoPFNePhCorrDR04 );
      tmpTMVAReader->AddVariable("SumDeltaRCharged",  &fMVAVar_MuDeltaRSumCharged     );
      tmpTMVAReader->AddVariable("DeltaRMeanCharged", &fMVAVar_MuDeltaRMeanCharged    );
    }
    if (type == kIsoRings){
      tmpTMVAReader->AddVariable("ChargedIso_DR0p0To0p1", &fMVAVar_ChargedIso_DR0p0To0p1); 
      tmpTMVAReader->AddVariable("ChargedIso_DR0p1To0p2", &fMVAVar_ChargedIso_DR0p1To0p2);
      tmpTMVAReader->AddVariable("ChargedIso_DR0p2To0p3", &fMVAVar_ChargedIso_DR0p2To0p3);
      tmpTMVAReader->AddVariable("ChargedIso_DR0p3To0p4", &fMVAVar_ChargedIso_DR0p3To0p4);
      tmpTMVAReader->AddVariable("ChargedIso_DR0p4To0p5", &fMVAVar_ChargedIso_DR0p4To0p5);
          
      tmpTMVAReader->AddVariable("GammaIso_DR0p0To0p1", &fMVAVar_GammaIso_DR0p0To0p1);
      tmpTMVAReader->AddVariable("GammaIso_DR0p1To0p2", &fMVAVar_GammaIso_DR0p1To0p2);
      tmpTMVAReader->AddVariable("GammaIso_DR0p2To0p3", &fMVAVar_GammaIso_DR0p2To0p3);
      tmpTMVAReader->AddVariable("GammaIso_DR0p3To0p4", &fMVAVar_GammaIso_DR0p3To0p4);
      tmpTMVAReader->AddVariable("GammaIso_DR0p4To0p5", &fMVAVar_GammaIso_DR0p4To0p5);
     
      tmpTMVAReader->AddVariable("NeutralHadronIso_DR0p0To0p1", &fMVAVar_NeutralHadronIso_DR0p0To0p1);
      tmpTMVAReader->AddVariable("NeutralHadronIso_DR0p1To0p2", &fMVAVar_NeutralHadronIso_DR0p1To0p2);
      tmpTMVAReader->AddVariable("NeutralHadronIso_DR0p2To0p3", &fMVAVar_NeutralHadronIso_DR0p2To0p3);
      tmpTMVAReader->AddVariable("NeutralHadronIso_DR0p3To0p4", &fMVAVar_NeutralHadronIso_DR0p3To0p4);
      tmpTMVAReader->AddVariable("NeutralHadronIso_DR0p4To0p5", &fMVAVar_NeutralHadronIso_DR0p4To0p5);

      //      tmpTMVAReader->AddSpectator("Pt",        &fMVAVar_Pt);
      //      tmpTMVAReader->AddSpectator("Eta",       &fMVAVar_Eta);
    }
    if (type == kIsoCombined){
      tmpTMVAReader->AddVariable("PFCharged",  &fMVAVar_MuRelIsoPFCharged );
      tmpTMVAReader->AddVariable("PFNeutral",  &fMVAVar_MuRelIsoPFNeutral );
      tmpTMVAReader->AddVariable("PFPhotons",  &fMVAVar_MuRelIsoPFPhotons );
      
      tmpTMVAReader->AddVariable("PFChargedDR0p0to0p1", &fMVAVar_ChargedIso_DR0p0To0p1); 
      tmpTMVAReader->AddVariable("PFChargedDR0p1to0p2", &fMVAVar_ChargedIso_DR0p1To0p2);
      tmpTMVAReader->AddVariable("PFChargedDR0p2to0p3", &fMVAVar_ChargedIso_DR0p2To0p3);
      tmpTMVAReader->AddVariable("PFChargedDR0p3to0p4", &fMVAVar_ChargedIso_DR0p3To0p4);
      tmpTMVAReader->AddVariable("PFChargedDR0p4to0p5", &fMVAVar_ChargedIso_DR0p4To0p5);
      
      tmpTMVAReader->AddVariable("PFNeutralDR0p0to0p1", &fMVAVar_NeutralHadronIso_DR0p0To0p1);
      tmpTMVAReader->AddVariable("PFNeutralDR0p1to0p2", &fMVAVar_NeutralHadronIso_DR0p1To0p2);
      tmpTMVAReader->AddVariable("PFNeutralDR0p2to0p3", &fMVAVar_NeutralHadronIso_DR0p2To0p3);
      tmpTMVAReader->AddVariable("PFNeutralDR0p3to0p4", &fMVAVar_NeutralHadronIso_DR0p3To0p4);
      tmpTMVAReader->AddVariable("PFNeutralDR0p4to0p5", &fMVAVar_NeutralHadronIso_DR0p4To0p5);
     
      tmpTMVAReader->AddVariable("PFPhotonsDR0p0to0p1", &fMVAVar_GammaIso_DR0p0To0p1);
      tmpTMVAReader->AddVariable("PFPhotonsDR0p1to0p2", &fMVAVar_GammaIso_DR0p1To0p2);
      tmpTMVAReader->AddVariable("PFPhotonsDR0p2to0p3", &fMVAVar_GammaIso_DR0p2To0p3);
      tmpTMVAReader->AddVariable("PFPhotonsDR0p3to0p4", &fMVAVar_GammaIso_DR0p3To0p4);
      tmpTMVAReader->AddVariable("PFPhotonsDR0p4to0p5", &fMVAVar_GammaIso_DR0p4To0p5);
     
      tmpTMVAReader->AddVariable("SumDeltaR",  &fMVAVar_MuDeltaRSum       );
      tmpTMVAReader->AddVariable("DeltaRMean", &fMVAVar_MuDeltaRMean      );
      tmpTMVAReader->AddVariable("Density",    &fMVAVar_MuDensity         );
      
      tmpTMVAReader->AddSpectator("Pt",        &fMVAVar_Pt);
      tmpTMVAReader->AddSpectator("Eta",       &fMVAVar_Eta);
    }
    
    tmpTMVAReader->BookMVA(fMethodname , weightfiles[i]);
   
    fTMVAReader.push_back(tmpTMVAReader);
  }
}


//--------------------------------------------------------------------------------------------------
Double_t MuonIsoMVA::MVAValue(Double_t MuPt, Double_t MuEta,
			      Double_t MuDZ,
			      Double_t MuIP2d, 
			      Double_t MuRelIsoPFCharged,
			      Double_t MuRelIsoPFNeutral,
			      Double_t MuRelIsoPFPhotons,
			      Double_t MuDeltaRSum,
			      Double_t MuDeltaRMean,
			      Double_t MuDensity){
  if (!fIsInitialized) { 
    cout << "Error: MuonIsoMVA not properly initialized.\n"; 
    return -9999;
  }
  
  //set all input variables
  fMVAVar_MuRelIsoPFCharged = MuRelIsoPFCharged;
  fMVAVar_MuRelIsoPFNeutral = MuRelIsoPFNeutral;
  fMVAVar_MuRelIsoPFPhotons = MuRelIsoPFPhotons;
  fMVAVar_MuDeltaRMean      = MuDeltaRMean;
  fMVAVar_MuDeltaRSum       = MuDeltaRSum; 
  fMVAVar_MuDensity         = MuDensity;
  fMVAVar_MuDZ              = MuDZ;
  fMVAVar_MuIP2d            = MuIP2d;

  Double_t mva = -9999;  
  TMVA::Reader *reader = 0;
  Int_t MVABin = -1;
  if (MuPt < 20  && fabs(MuEta) <  1.479) MVABin = 0;
  if (MuPt < 20  && fabs(MuEta) >= 1.479) MVABin = 1;
  if (MuPt >= 20 && fabs(MuEta) <  1.479) MVABin = 2;
  if (MuPt >= 20 && fabs(MuEta) >= 1.479) MVABin = 3;
  
  if (MVABin < 0) return mva;
  if (MVABin > 3) return mva;

  reader = fTMVAReader[MVABin];
  mva = reader->EvaluateMVA( fMethodname );
  
  float muoniso = (mva + 0.7)*0.4/1.4;
  //  float muoniso = (mva + 1.0)*0.4/2.0;
  
#ifdef DEBUG
  cout << "Muon with pT = " << MuPt << " and eta = " << MuEta << endl;
  cout << "MVA value    = " << mva  << endl;
#endif
  return muoniso;
}


Double_t MuonIsoMVA::MVAValue(Double_t MuPt, Double_t MuEta,
			      Double_t MuRelIsoPFCharged,
			      Double_t MuRelIsoPFNeutral,
			      Double_t MuRelIsoPFPhotons,
			      Double_t MuDeltaRSum,
			      Double_t MuDeltaRMean,
			      Double_t MuDensity){
  if (!fIsInitialized) { 
    cout << "Error: MuonIsoMVA not properly initialized.\n"; 
    return -9999;
  }
  
  //set all input variables
  fMVAVar_MuRelIsoPFCharged = MuRelIsoPFCharged;
  fMVAVar_MuRelIsoPFNeutral = MuRelIsoPFNeutral;
  fMVAVar_MuRelIsoPFPhotons = MuRelIsoPFPhotons;
  fMVAVar_MuDeltaRMean      = MuDeltaRMean;
  fMVAVar_MuDeltaRSum       = MuDeltaRSum; 
  fMVAVar_MuDensity         = MuDensity;

  Double_t mva = -9999;  
  TMVA::Reader *reader = 0;
  Int_t MVABin = -1;
  if (MuPt < 20  && fabs(MuEta) <  1.479) MVABin = 0;
  if (MuPt < 20  && fabs(MuEta) >= 1.479) MVABin = 1;
  if (MuPt >= 20 && fabs(MuEta) <  1.479) MVABin = 2;
  if (MuPt >= 20 && fabs(MuEta) >= 1.479) MVABin = 3;
  
  reader = fTMVAReader[MVABin];
  mva = reader->EvaluateMVA( fMethodname );
  
  float muoniso = (mva + 0.6)*0.4/1.2;
  
#ifdef DEBUG
  cout << "Muon with pT = " << MuPt << " and eta = " << MuEta << endl;
  cout << "MVA value    = " << mva  << endl;
#endif
  return muoniso;
}
Double_t MuonIsoMVA::MVAValue(Double_t MuPt, Double_t MuEta,
			      Double_t MuRelIsoPFCharged,
			      Double_t MuRelIsoPFNePh,
			      Double_t MuDeltaRSumCharged,
			      Double_t MuDeltaRMeanCharged){

  if (!fIsInitialized) { 
    cout << "Error: MuonIsoMVA not properly initialized.\n"; 
    return -9999;
  }
  
  //set all input variables
  fMVAVar_MuRelIsoPFChargedDR04  = MuRelIsoPFCharged;
  fMVAVar_MuRelIsoPFNePhCorrDR04 = MuRelIsoPFNePh;
  fMVAVar_MuDeltaRMeanCharged    = MuDeltaRMeanCharged;
  fMVAVar_MuDeltaRSumCharged     = MuDeltaRSumCharged; 


  Double_t mva = -9999;  
  TMVA::Reader *reader = 0;
  Int_t MVABin = -1;
  if (MuPt < 20  && fabs(MuEta) <  1.479) MVABin = 0;
  if (MuPt < 20  && fabs(MuEta) >= 1.479) MVABin = 1;
  if (MuPt >= 20 && fabs(MuEta) <  1.479) MVABin = 2;
  if (MuPt >= 20 && fabs(MuEta) >= 1.479) MVABin = 3;
  
  reader = fTMVAReader[MVABin];
  mva = reader->EvaluateMVA( fMethodname );
   
  float muoniso = (mva + 1.0)*0.4/2.0;
  
#ifdef DEBUG
  cout << "Muon with pT = " << MuPt << " and eta = " << MuEta << endl;
  cout << "MVA value    = " << mva  << endl;
#endif
  return muoniso;
}

double MuonIsoMVA::MVAValue(double MuPt, double MuEta,
			    Bool_t   IsGlobal, Bool_t IsTracker,
			    double ch00to01, double ch01to02, double ch02to03, double ch03to04, double ch04to05, 
			    double ne00to01, double ne01to02, double ne02to03, double ne03to04, double ne04to05, 
			    double ph00to01, double ph01to02, double ph02to03, double ph03to04, double ph04to05){

  if (!fIsInitialized) { 
    cout << "Error: MuonIsoMVA not properly initialized.\n"; 
    return -9999;
  }

  fMVAVar_ChargedIso_DR0p0To0p1       = ch00to01;
  fMVAVar_ChargedIso_DR0p1To0p2       = ch01to02;
  fMVAVar_ChargedIso_DR0p2To0p3       = ch02to03;
  fMVAVar_ChargedIso_DR0p3To0p4       = ch03to04;
  fMVAVar_ChargedIso_DR0p4To0p5       = ch04to05;
  fMVAVar_GammaIso_DR0p0To0p1         = ph00to01;
  fMVAVar_GammaIso_DR0p1To0p2         = ph01to02;
  fMVAVar_GammaIso_DR0p2To0p3         = ph02to03; 
  fMVAVar_GammaIso_DR0p3To0p4         = ph03to04;
  fMVAVar_GammaIso_DR0p4To0p5         = ph04to05;
  fMVAVar_NeutralHadronIso_DR0p0To0p1 = ne00to01;
  fMVAVar_NeutralHadronIso_DR0p1To0p2 = ne01to02;
  fMVAVar_NeutralHadronIso_DR0p2To0p3 = ne02to03;
  fMVAVar_NeutralHadronIso_DR0p3To0p4 = ne03to04;
  fMVAVar_NeutralHadronIso_DR0p4To0p5 = ne04to05;

  Double_t mva = -9999;  
  TMVA::Reader *reader = 0;
  Int_t MVABin = -1;
  if (IsGlobal && IsTracker){
    if (MuPt < 10  && fabs(MuEta) <  1.479) MVABin = 0;
    if (MuPt < 10  && fabs(MuEta) >= 1.479) MVABin = 1;
    if (MuPt >= 10 && fabs(MuEta) <  1.479) MVABin = 2;
    if (MuPt >= 10 && fabs(MuEta) >= 1.479) MVABin = 3;
  }
  else if (IsGlobal) {
    MVABin = 4;
  }
  else if (IsTracker){
    MVABin = 5;
  }
  
  reader = fTMVAReader[MVABin];
  mva = reader->EvaluateMVA( fMethodname );
  
  //  float muoniso = (mva + 0.6)*0.4/1.2;
  float muoniso = (mva + 1.0)*0.4/2.0;
  return mva;
}

double MuonIsoMVA::MVAValue(double MuPt, double MuEta,
			    double ch00to01, double ch01to02, double ch02to03, double ch03to04, double ch04to05, 
			    double ne00to01, double ne01to02, double ne02to03, double ne03to04, double ne04to05, 
			    double ph00to01, double ph01to02, double ph02to03, double ph03to04, double ph04to05,
			    double drsum,    double drmean,   double density){

  if (!fIsInitialized) { 
    cout << "Error: MuonIsoMVA not properly initialized.\n"; 
    return -9999;
  }

  fMVAVar_MuRelIsoPFCharged = TMath::Min(2.5, (double)(ch00to01 + ch01to02 + ch02to03 + ch03to04 + ch04to05));
  fMVAVar_MuRelIsoPFNeutral = TMath::Min(2.5, (double)(ne00to01 + ne01to02 + ne02to03 + ne03to04 + ne04to05));
  fMVAVar_MuRelIsoPFPhotons = TMath::Min(2.5, (double)(ph00to01 + ph01to02 + ph02to03 + ph03to04 + ph04to05));
  fMVAVar_MuDeltaRMean      = drsum;
  fMVAVar_MuDeltaRSum       = drmean;
  fMVAVar_MuDensity         = density;

  fMVAVar_ChargedIso_DR0p0To0p1       = ch00to01;
  fMVAVar_ChargedIso_DR0p1To0p2       = ch01to02;
  fMVAVar_ChargedIso_DR0p2To0p3       = ch02to03;
  fMVAVar_ChargedIso_DR0p3To0p4       = ch03to04;
  fMVAVar_ChargedIso_DR0p4To0p5       = ch04to05;
  fMVAVar_GammaIso_DR0p0To0p1         = ph00to01;
  fMVAVar_GammaIso_DR0p1To0p2         = ph01to02;
  fMVAVar_GammaIso_DR0p2To0p3         = ph02to03; 
  fMVAVar_GammaIso_DR0p3To0p4         = ph03to04;
  fMVAVar_GammaIso_DR0p4To0p5         = ph04to05;
  fMVAVar_NeutralHadronIso_DR0p0To0p1 = ne00to01;
  fMVAVar_NeutralHadronIso_DR0p1To0p2 = ne01to02;
  fMVAVar_NeutralHadronIso_DR0p2To0p3 = ne02to03;
  fMVAVar_NeutralHadronIso_DR0p3To0p4 = ne03to04;
  fMVAVar_NeutralHadronIso_DR0p4To0p5 = ne04to05;


  Double_t mva = -9999;  
  TMVA::Reader *reader = 0;
  Int_t MVABin = -1;
  if (MuPt < 20  && fabs(MuEta) <  1.479) MVABin = 0;
  if (MuPt < 20  && fabs(MuEta) >= 1.479) MVABin = 1;
  if (MuPt >= 20 && fabs(MuEta) <  1.479) MVABin = 2;
  if (MuPt >= 20 && fabs(MuEta) >= 1.479) MVABin = 3;

  if (MVABin < 0) return mva;
  if (MVABin > 3) return mva;

  reader = fTMVAReader[MVABin];
  mva = reader->EvaluateMVA( fMethodname );
  
  float muoniso = (mva + 0.6)*0.4/1.2;
  //  float muoniso = (mva + 1.0)*0.4/2.0;
  return muoniso;
}
