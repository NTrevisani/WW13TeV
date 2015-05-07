//--------------------------------------------------------------------------------------------------
// $Id $
//
// MuonIsoMVA
//
// Helper Class for applying MVA muon Iso selection
//
// Authors: S. Folgueras, E. Di Marco
//--------------------------------------------------------------------------------------------------

#ifndef MUONISOMVA_H
#define MUONISOMVA_H

#include "TMVA/Tools.h"
#include "TMVA/Reader.h"

class MuonIsoMVA {
 public:
  enum MVAType {
    kIsoDeltaR,
    kIsoDeltaRV1,
    kIsoRings,
    kIsoCombined
  };
  enum MVAbin{
    kLowPtBarrel,
    kLowPtEndcap,
    kHighPtBarrel,
    kHighPtEndCap
  };
  
  MuonIsoMVA();
  ~MuonIsoMVA(); 
  
  void   Initialize(std::string methodName,
		    MVAType type,		    
		    std::string Pt10To20BarrelWeights , 
		    std::string Pt10To20EndcapWeights , 
		    std::string Pt20ToInfBarrelWeights , 
		    std::string Pt20ToInfEndcapWeights);

  void   Initialize(std::string methodName,
		    MVAType type,		    
		    std::vector<std::string> weightfiles);
  
  Bool_t IsInitialized() const { return fIsInitialized; }
  
  double MVAValue(double MuPt, double MuEta,
		  double MuDZ,
		  double MuIP2d,
		  double MuRelIsoPFCharged,
		  double MuRelIsoPFNeutral,
		  double MuRelIsoPFPhotons,
		  double MuDeltaRSum,
		  double MuDeltaRMean,
		  double MuDensity);
  
  // DeltaR only
  double MVAValue(double MuPt, double MuEta,
		  double MuRelIsoPFCharged,
		  double MuRelIsoPFNeutral,
		  double MuRelIsoPFPhotons,
		  double MuDeltaRSum,
		  double MuDeltaRMean,
		  double MuDensity);
  
  // DeltaRV1
  double MVAValue(double MuPt, double MuEta,
		  double MuRelIsoPFCharged,
		  double MuRelIsoPFNePhCorr,
		  double MuDeltaRSumCharged,
		  double MuDeltaRMeanCharged);
  
  
  // IsoRings only
  double MVAValue(double MuPt, double MuEta,
		  bool,   bool,
		  double, double, double, double, double, 
		  double, double, double, double, double, 
		  double, double, double, double, double);

  
  // IsoRings + DeltaR
  double MVAValue(double MuPt, double MuEta,
		  double, double, double, double, double, 
		  double, double, double, double, double,
		  double, double, double, double, double,
		  double, double, double); 

 protected:
  std::vector<TMVA::Reader*> fTMVAReader;
  std::string                fMethodname;
  
  Bool_t                    fIsInitialized;

  //isolation
  Float_t                   fMVAVar_Pt;
  Float_t                   fMVAVar_Eta;

  Float_t                   fMVAVar_ChargedIso_DR0p0To0p1;
  Float_t                   fMVAVar_ChargedIso_DR0p1To0p2;
  Float_t                   fMVAVar_ChargedIso_DR0p2To0p3;
  Float_t                   fMVAVar_ChargedIso_DR0p3To0p4;
  Float_t                   fMVAVar_ChargedIso_DR0p4To0p5;
  Float_t                   fMVAVar_GammaIso_DR0p0To0p1;
  Float_t                   fMVAVar_GammaIso_DR0p1To0p2;
  Float_t                   fMVAVar_GammaIso_DR0p2To0p3;
  Float_t                   fMVAVar_GammaIso_DR0p3To0p4;
  Float_t                   fMVAVar_GammaIso_DR0p4To0p5;
  Float_t                   fMVAVar_NeutralHadronIso_DR0p0To0p1;
  Float_t                   fMVAVar_NeutralHadronIso_DR0p1To0p2;
  Float_t                   fMVAVar_NeutralHadronIso_DR0p2To0p3;
  Float_t                   fMVAVar_NeutralHadronIso_DR0p3To0p4;
  Float_t                   fMVAVar_NeutralHadronIso_DR0p4To0p5;
    
  Float_t                   fMVAVar_MuRelIsoPFCharged;
  Float_t                   fMVAVar_MuRelIsoPFNeutral;
  Float_t                   fMVAVar_MuRelIsoPFPhotons;
  Float_t                   fMVAVar_MuDeltaRMean;
  Float_t                   fMVAVar_MuDeltaRSum;
  Float_t                   fMVAVar_MuDensity;
  Float_t                   fMVAVar_MuDZ;
  Float_t                   fMVAVar_MuIP2d;
  
  Float_t                   fMVAVar_MuRelIsoPFChargedDR04;
  Float_t                   fMVAVar_MuRelIsoPFNePhCorrDR04;
  Float_t                   fMVAVar_MuDeltaRSumCharged;
  Float_t                   fMVAVar_MuDeltaRMeanCharged;
};

#endif
