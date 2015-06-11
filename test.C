////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////

/////           DUMMY                                      /////

////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////


#include "test.h"
#include "TH1.h"
#include "TH2.h"
#include "TKey.h"
#include "TFile.h"
#include "TChain.h"
#include "TTree.h"
#include "TRandom.h"
#include "TString.h"

#include "TLorentzVector.h"
#include <vector>
#include "TROOT.h"
#include <iostream>

#include "TDatabasePDG.h"

//lepton ID
bool IsTightLepton(int k,
		   std::vector<float>* vector_lepton_id,  
		   std::vector<float>* vector_lepton_isTightMuon,
		   std::vector<float>* vector_electron_deltaEtaIn,
		   std::vector<float>* vector_electron_sigmaIetaIeta,
		   std::vector<float>* vector_electron_HoE,
		   std::vector<float>* vector_electron_d0,
		   std::vector<float>* vector_electron_dz,
		   std::vector<float>* vector_electron_ooEooP,
		   std::vector<float>* vector_electron_passConversion,
		   std::vector<float>* vector_electron_scEta,
		   std::vector<float>* vector_electron_deltaPhiIn
		   )
{
  bool is_tight_lepton = false;

  // Muon tight ID
  if (fabs(vector_lepton_id->at(k)) == 13)
    {
      is_tight_lepton = vector_lepton_isTightMuon->at(k);
    }
  // Electron cut based medium ID
  else if (fabs(vector_lepton_id->at(k)) == 11)
    {
      float aeta = fabs(vector_electron_scEta->at(k));

      if (aeta <= 1.479)
	{                                                   //V1         //V2
	  if (fabs(vector_electron_deltaEtaIn->at(k)) < 0.007641 &&  //0.008925 &&
	      fabs(vector_electron_deltaPhiIn->at(k)) < 0.032643 &&  //0.035973 &&
	      vector_electron_sigmaIetaIeta->at(k)    < 0.010399 &&  //0.009996 &&
	      vector_electron_HoE->at(k)              < 0.060662 &&  //0.050537 &&
	      fabs(vector_electron_d0->at(k))         < 0.011811 &&  //0.012235 &&
	      fabs(vector_electron_dz->at(k))         < 0.070775 &&  //0.042020 &&
	      fabs(vector_electron_ooEooP->at(k))     < 0.153897 &&  //0.091942 &&
	      !vector_electron_passConversion->at(k))  // Includes expectedMissingInnerHits
	    {
	      is_tight_lepton = true;
	    }
	}
      else if (aeta > 1.479 && aeta < 2.5)
	{                                                   //V1         //V2
	  if (fabs(vector_electron_deltaEtaIn->at(k)) < 0.009285 &&  //0.007429 &&
	      fabs(vector_electron_deltaPhiIn->at(k)) < 0.042447 &&  //0.067879 &&
	      vector_electron_sigmaIetaIeta->at(k)    < 0.029524 &&  //0.030135 &&
	      vector_electron_HoE->at(k)              < 0.104263 &&  //0.086782 &&
	      fabs(vector_electron_d0->at(k))         < 0.051682 &&  //0.036719 &&
	      fabs(vector_electron_dz->at(k))         < 0.180720 &&  //0.138142 &&
	      fabs(vector_electron_ooEooP->at(k))     < 0.137468 &&  //0.100683 &&
	      !vector_electron_passConversion->at(k))  // Includes expectedMissingInnerHits
	    {
	      is_tight_lepton = true;
	    }
	}
    }

  return is_tight_lepton;
}

//lepton isolation
bool IsIsolatedLepton(int k,
		      std::vector<float>* vector_lepton_chargedHadronIso,
		      std::vector<float>* vector_lepton_neutralHadronIso,
		      std::vector<float>* vector_lepton_sumPUPt,
		      std::vector<float>* vector_lepton_pt,
		      std::vector<float>* vector_lepton_id,
		      std::vector<float>* vector_lepton_eta,
		      std::vector<float>* vector_lepton_photonIso,
		      float jetRho_
)
{
  float pt = vector_lepton_pt->at(k);
  float id = vector_lepton_id->at(k);

  float isolation = 999;

  bool is_isolated_lepton = false;

  if (fabs(id) == 13)
    {
      isolation =
	vector_lepton_chargedHadronIso->at(k) +
	max(float(0.0),
	    float(vector_lepton_photonIso->at(k) +
		  vector_lepton_neutralHadronIso->at(k) -
		  0.5*vector_lepton_sumPUPt->at(k)));

      is_isolated_lepton = (isolation/pt < 0.12);
    }
  else if (fabs(id) == 11)
    {
      float aeta = fabs(vector_lepton_eta->at(k));
      
      float effective_area = -999;
      
      if      (aeta >  2.2)               effective_area = 0.2680;
      else if (aeta >= 2.0 && aeta < 2.2) effective_area = 0.1565;
      else if (aeta >= 1.3 && aeta < 2.0) effective_area = 0.1077;
      else if (aeta >= 0.8 && aeta < 1.3) effective_area = 0.1734;
      else if (aeta <  0.8)               effective_area = 0.1830;

      isolation =
	vector_lepton_chargedHadronIso->at(k) +
	max(float(0.0),
	    float(vector_lepton_photonIso->at(k) +
		  vector_lepton_neutralHadronIso->at(k) -
		  jetRho_*effective_area));

      is_isolated_lepton = (isolation/pt < 0.15);
    }
  
  return is_isolated_lepton;
}


test::test(TTree* tree):
  PAFAnalysis(tree) {
}


void test::Initialise() {

  GetInputParameters()->TheNamedBool("IsDATA", IsDATA);
  GetInputParameters()->TheNamedInt("NEvents", NEvents);
  GetInputParameters()->TheNamedDouble("Luminosity", Luminosity);
  GetInputParameters()->TheNamedDouble("XSection", XSection);
  GetInputParameters()->TheNamedInt("WhichRun", WhichRun);
  GetInputParameters()->TheNamedInt("jetChannel", jetChannel);
  
  TString TheSample = GetInputParameters()->TheNamedString("theSample");
  TString flavorChannel = GetInputParameters()->TheNamedString("FlavorChannel");
  
  cout<<"la stringa che voglio passare si chiama: "<<TheSample<<endl;
  
  
  // To do only once
  float LuminosityPU = 0;
  //fInputParameters->TheNamedFloat("LuminosityPU",LuminosityPU);
  //fPUWeight = new PUWeight(LuminosityPU, Spring11);//Summer11InTime);
  
  // Set the channel                                                                                                                                                               
  //----------------------------------------------------------------------------                                                                                                   
  SelectedChannel = -999;
  
  if      (flavorChannel == "MuMu") SelectedChannel =  0;
  else if (flavorChannel == "EE"  ) SelectedChannel =  1;
  else if (flavorChannel == "EMu" ) SelectedChannel =  2;
  else if (flavorChannel == "MuE" ) SelectedChannel =  3;
  else if (flavorChannel == "OF" )  SelectedChannel =  4;
  else if (flavorChannel == "SF" )  SelectedChannel =  5;
  
  else if (flavorChannel == "All" ) SelectedChannel = -1;
  
  G_Debug_DefineAnalysisVariables = false;
  
  //------------------------------------------------------------------------------
  // Create histos
  //------------------------------------------------------------------------------
  
  h_n_PV = CreateH1F("h_n_PV","h_n_PV",50,0.,10.);

  // Counting histograms    
  //---------------------------------------------------------------------------- 
  
  hWTrigger     = CreateH1F("hWTrigger",     "", 10, 0, 10);
  hWMetCut      = CreateH1F("hWMetCut",      "", 10, 0, 10);
  hWLowMinv     = CreateH1F("hWLowMinv",     "", 10, 0, 10);
  hWZVeto       = CreateH1F("hWZVeto",       "", 10, 0, 10);
  hWpMetCut     = CreateH1F("hWpMetCut",     "", 10, 0, 10);
  hWJetVeto     = CreateH1F("hWJetVeto",     "", 10, 0, 10);
  hWnJets       = CreateH1F("hWnJets",         "", 10, 0, 10);
  hWeffnJets    = CreateH1F("hWeffnJets",      "", 10, 0, 10);
  hWnBtaggedJets     = CreateH1F("hWnBtaggedJets",          "", 10, 0, 10);
  hWeffnBtaggedJets  = CreateH1F("hWeffnBtaggedJets",       "", 10, 0, 10);
  hWnJetsBveto    = CreateH1F("hWnJetsBveto",         "", 10, 0, 10);
  hWeffnJetsBveto = CreateH1F("hWeffnJetsBveto",      "", 10, 0, 10);
  hWnJetsBvetoAfterHt    = CreateH1F("hWnJetsBvetoAfterHt",         "", 10, 0, 10);
  hWeffnJetsBvetoAfterHt = CreateH1F("hWeffnJetsBvetoAfterHt",      "", 10, 0, 10);
  
  hWDeltaPhiJet = CreateH1F("hWDeltaPhiJet", "", 10, 0, 10);
  hWSoftMuVeto  = CreateH1F("hWSoftMuVeto",  "", 10, 0, 10);
  hWExtraLepton = CreateH1F("hWExtraLepton", "", 10, 0, 10);
  hWPtll        = CreateH1F("hWPtll",        "", 10, 0, 10);
  hWTopTagging  = CreateH1F("hWTopTagging",  "", 10, 0, 10);

  hWeffTrigger     = CreateH1F("hWeffTrigger",     "", 10, 0, 10);
  hWeffMetCut      = CreateH1F("hWeffMetCut",      "", 10, 0, 10);
  hWeffLowMinv     = CreateH1F("hWeffLowMinv",     "", 10, 0, 10);
  hWeffZVeto       = CreateH1F("hWeffZVeto",       "", 10, 0, 10);
  hWeffpMetCut     = CreateH1F("hWeffpMetCut",     "", 10, 0, 10);
  hWeffJetVeto     = CreateH1F("hWeffJetVeto",     "", 10, 0, 10);
  hWeffDeltaPhiJet = CreateH1F("hWeffDeltaPhiJet", "", 10, 0, 10);
  hWeffSoftMuVeto  = CreateH1F("hWeffSoftMuVeto",  "", 10, 0, 10);
  hWeffExtraLepton = CreateH1F("hWeffExtraLepton", "", 10, 0, 10);
  hWeffPtll        = CreateH1F("hWeffPtll",        "", 10, 0, 10);
  hWeffTopTagging  = CreateH1F("hWeffTopTagging",  "", 10, 0, 10);
  
  // WW level histograms     
  //----------------------------------------------------------------------------  
  
  for (Int_t nC=0; nC<4; nC++) {

    hPtLepton1WWLevel[nC]       = CreateH1F(Form("hPtLepton1WWLevel%.1i", nC),       "", 200, 0, 200);
    hPtLepton2WWLevel[nC]       = CreateH1F(Form("hPtLepton2WWLevel%.1i", nC),       "", 200, 0, 200);
    hPtDiLeptonWWLevel[nC]      = CreateH1F(Form("hPtDiLeptonWWLevel%.1i", nC),      "", 200, 0, 200);
    hMinvWWLevel[nC]            = CreateH1F(Form("hMinvWWLevel%.1i", nC),            "", 200, 0, 200);
    hMtWWLevel[nC]              = CreateH1F(Form("hMtWWLevel%.1i", nC),              "", 250, 0, 250);
    hNJets30WWLevel[nC]         = CreateH1F(Form("hNJetsPF30WWLevel%.1i", nC),       "",  10, 0,  10);
    hpfMetWWLevel[nC]           = CreateH1F(Form("hpfMetWWLevel%.1i", nC),           "", 150, 0, 150);
    hppfMetWWLevel[nC]          = CreateH1F(Form("hppfMetWWLevel%.1i", nC),          "", 150, 0, 150);
    hchMetWWLevel[nC]           = CreateH1F(Form("hchMetWWLevel%.1i", nC),           "", 150, 0, 150);
    hpchMetWWLevel[nC]          = CreateH1F(Form("hpchMetWWLevel%.1i", nC),          "", 150, 0, 150);
    hpminMetWWLevel[nC]         = CreateH1F(Form("hpminMetWWLevel%.1i", nC),         "", 150, 0, 150);
    hDeltaRLeptonsWWLevel[nC]   = CreateH1F(Form("hDeltaRLeptonsWWLevel%.1i", nC),   "",  50, 0,   5);
    hDeltaPhiLeptonsWWLevel[nC] = CreateH1F(Form("hDeltaPhiLeptonsWWLevel%.1i", nC), "",  32, 0, 3.2);
    hDPhiPtllJetWWLevel[nC]     = CreateH1F(Form("hDPhiPtllJetWWLevel%.1i", nC),     "",  32, 0, 3.2);
    hDataEvents[nC]             = CreateH1F(Form("hDataEvents%.1i", nC),             "",  32, 0, 3.2);
    hBackgroundEvents[nC]       = CreateH1F(Form("hBackgroundEvents%.1i", nC),       "",  32, 0, 3.2);
    hEff[nC]                    = CreateH1F(Form("hEff%.1i", nC),                     "",  100, 0, 10);
    hSigEl[nC]                  = CreateH1F(Form("hSigEl%.1i", nC),                   "", 1000,0.,500);
    hSigMu[nC]                  = CreateH1F(Form("hSigMu%.1i", nC),                   "", 1000,0.,500);

    hPtLepton1WWLevelNoHt[nC]       = CreateH1F(Form("hPtLepton1WWLevelNoHt%.1i", nC),       "", 200, 0, 200);
    hPtLepton2WWLevelNoHt[nC]       = CreateH1F(Form("hPtLepton2WWLevelNoHt%.1i", nC),       "", 200, 0, 200);
    hPtDiLeptonWWLevelNoHt[nC]      = CreateH1F(Form("hPtDiLeptonWWLevelNoHt%.1i", nC),      "", 200, 0, 200);
    hMinvWWLevelNoHt[nC]            = CreateH1F(Form("hMinvWWLevelNoHt%.1i", nC),            "", 200, 0, 200);
    hMtWWLevelNoHt[nC]              = CreateH1F(Form("hMtWWLevelNoHt%.1i", nC),              "", 250, 0, 250);
    hNJets30WWLevelNoHt[nC]         = CreateH1F(Form("hNJetsPF30WWLevelNoHt%.1i", nC),       "",  10, 0,  10);
    hpfMetWWLevelNoHt[nC]           = CreateH1F(Form("hpfMetWWLevelNoHt%.1i", nC),           "", 150, 0, 150);
    hppfMetWWLevelNoHt[nC]          = CreateH1F(Form("hppfMetWWLevelNoHt%.1i", nC),          "", 150, 0, 150);
    hchMetWWLevelNoHt[nC]           = CreateH1F(Form("hchMetWWLevelNoHt%.1i", nC),           "", 150, 0, 150);
    hpchMetWWLevelNoHt[nC]          = CreateH1F(Form("hpchMetWWLevelNoHt%.1i", nC),          "", 150, 0, 150);
    hpminMetWWLevelNoHt[nC]         = CreateH1F(Form("hpminMetWWLevelNoHt%.1i", nC),         "", 150, 0, 150);
    hDeltaRLeptonsWWLevelNoHt[nC]   = CreateH1F(Form("hDeltaRLeptonsWWLevelNoHt%.1i", nC),   "",  50, 0,   5);
    hDeltaPhiLeptonsWWLevelNoHt[nC] = CreateH1F(Form("hDeltaPhiLeptonsWWLevelNoHt%.1i", nC), "",  32, 0, 3.2);
    hDPhiPtllJetWWLevelNoHt[nC]     = CreateH1F(Form("hDPhiPtllJetWWLevelNoHt%.1i", nC),     "",  32, 0, 3.2);
    hDataEventsNoHt[nC]             = CreateH1F(Form("hDataEventsNoHt%.1i", nC),             "",  32, 0, 3.2);
    hBackgroundEventsNoHt[nC]       = CreateH1F(Form("hBackgroundEventsNoHt%.1i", nC),       "",  32, 0, 3.2);
    hEffNoHt[nC]                    = CreateH1F(Form("hEffNoHt%.1i", nC),                     "",  100, 0, 10);
    hSigElNoHt[nC]                  = CreateH1F(Form("hSigElNoHt%.1i", nC),                   "", 1000,0.,500);
    hSigMuNoHt[nC]                  = CreateH1F(Form("hSigMuNoHt%.1i", nC),                   "", 1000,0.,500);

    hPtLepton1WWLevelHtPlus[nC]       = CreateH1F(Form("hPtLepton1WWLevelHtPlus%.1i", nC),       "", 200, 0, 200);
    hPtLepton2WWLevelHtPlus[nC]       = CreateH1F(Form("hPtLepton2WWLevelHtPlus%.1i", nC),       "", 200, 0, 200);
    hPtDiLeptonWWLevelHtPlus[nC]      = CreateH1F(Form("hPtDiLeptonWWLevelHtPlus%.1i", nC),      "", 200, 0, 200);
    hMinvWWLevelHtPlus[nC]            = CreateH1F(Form("hMinvWWLevelHtPlus%.1i", nC),            "", 200, 0, 200);
    hMtWWLevelHtPlus[nC]              = CreateH1F(Form("hMtWWLevelHtPlus%.1i", nC),              "", 250, 0, 250);
    hNJets30WWLevelHtPlus[nC]         = CreateH1F(Form("hNJetsPF30WWLevelHtPlus%.1i", nC),       "",  10, 0,  10);
    hpfMetWWLevelHtPlus[nC]           = CreateH1F(Form("hpfMetWWLevelHtPlus%.1i", nC),           "", 150, 0, 150);
    hppfMetWWLevelHtPlus[nC]          = CreateH1F(Form("hppfMetWWLevelHtPlus%.1i", nC),          "", 150, 0, 150);
    hchMetWWLevelHtPlus[nC]           = CreateH1F(Form("hchMetWWLevelHtPlus%.1i", nC),           "", 150, 0, 150);
    hpchMetWWLevelHtPlus[nC]          = CreateH1F(Form("hpchMetWWLevelHtPlus%.1i", nC),          "", 150, 0, 150);
    hpminMetWWLevelHtPlus[nC]         = CreateH1F(Form("hpminMetWWLevelHtPlus%.1i", nC),         "", 150, 0, 150);
    hDeltaRLeptonsWWLevelHtPlus[nC]   = CreateH1F(Form("hDeltaRLeptonsWWLevelHtPlus%.1i", nC),   "",  50, 0,   5);
    hDeltaPhiLeptonsWWLevelHtPlus[nC] = CreateH1F(Form("hDeltaPhiLeptonsWWLevelHtPlus%.1i", nC), "",  32, 0, 3.2);
    hDPhiPtllJetWWLevelHtPlus[nC]     = CreateH1F(Form("hDPhiPtllJetWWLevelHtPlus%.1i", nC),     "",  32, 0, 3.2);
    hDataEventsHtPlus[nC]             = CreateH1F(Form("hDataEventsHtPlus%.1i", nC),             "",  32, 0, 3.2);
    hBackgroundEventsHtPlus[nC]       = CreateH1F(Form("hBackgroundEventsHtPlus%.1i", nC),       "",  32, 0, 3.2);
    hEffHtPlus[nC]                    = CreateH1F(Form("hEffHtPlus%.1i", nC),                     "",  100, 0, 10);
    hSigElHtPlus[nC]                  = CreateH1F(Form("hSigElHtPlus%.1i", nC),                   "", 1000,0.,500);
    hSigMuHtPlus[nC]                  = CreateH1F(Form("hSigMuHtPlus%.1i", nC),                   "", 1000,0.,500);

    hHt[nC]                     = CreateH1F(Form("hHt%.1i",               nC),       "", 3000, 0, 3000);
    hHtAfter[nC]                = CreateH1F(Form("hHtAfter%.1i",          nC),       "", 3000, 0, 3000);

  }

  h_WWLevel_TightFailEvents  = CreateH1F("h_WWLevel_TightFailEvents",  "", 10, 0 , 10);
  h_WWLevel_TightTightEvents = CreateH1F("h_WWLevel_TightTightEvents", "", 10, 0 , 10);
  h_WWLevel_TightLooseEvents = CreateH1F("h_WWLevel_TightLooseEvents", "", 10, 0 , 10);
  
  
  // TwoLeptons level histograms   
  //----------------------------------------------------------------------------
  
  hPtLepton1TwoLeptonsLevel       = CreateH1F("hPtLepton1TwoLeptonsLevel",       "", 200, 0, 200);
  hPtLepton2TwoLeptonsLevel       = CreateH1F("hPtLepton2TwoLeptonsLevel",       "", 200, 0, 200);
  hPtDiLeptonTwoLeptonsLevel      = CreateH1F("hPtDiLeptonTwoLeptonsLevel",      "", 200, 0, 200);
  hMinvTwoLeptonsLevel            = CreateH1F("hMinvTwoLeptonsLevel",            "", 200, 0, 200);
  hMtTwoLeptonsLevel              = CreateH1F("hMtTwoLeptonsLevel",              "", 250, 0, 250);
  hNJets30TwoLeptonsLevel         = CreateH1F("hNJets30TwoLeptonsLevel",         "",  10, 0,  10);
  hpfMetTwoLeptonsLevel           = CreateH1F("hpfMetTwoLeptonsLevel",           "", 150, 0, 150);
  hppfMetTwoLeptonsLevel          = CreateH1F("hppfMetTwoLeptonsLevel",          "", 150, 0, 150);
  hchMetTwoLeptonsLevel           = CreateH1F("hchMetTwoLeptonsLevel",           "", 150, 0, 150);
  hpchMetTwoLeptonsLevel          = CreateH1F("hpchMetTwoLeptonsLevel",          "", 150, 0, 150);
  hpminMetTwoLeptonsLevel         = CreateH1F("hpminMetTwoLeptonsLevel",         "", 150, 0, 150);
  hDeltaRLeptonsTwoLeptonsLevel   = CreateH1F("hDeltaRLeptonsTwoLeptonsLevel",   "",  50, 0,   5);
  hDeltaPhiLeptonsTwoLeptonsLevel = CreateH1F("hDeltaPhiLeptonsTwoLeptonsLevel", "",  32, 0, 3.2);
  hDPhiPtllJetTwoLeptonsLevel     = CreateH1F("hDPhiPtllJetTwoLeptonsLevel",     "",  32, 0, 3.2);
  
  h_TwoLeptons_TightFailEvents  = CreateH1F("h_TwoLeptons_TightFailEvents",  "", 10, 0 , 10);
  h_TwoLeptons_TightTightEvents = CreateH1F("h_TwoLeptons_TightTightEvents", "", 10, 0 , 10);
  h_TwoLeptons_TightLooseEvents = CreateH1F("h_TwoLeptons_TightLooseEvents", "", 10, 0 , 10);

  //Isolation Plots
  //-----------------------------------------------------------------------------

  hIsoMu = CreateH1F("hIsoMu","", 100000, 0., 100);
  hIsoEl = CreateH1F("hIsoEl","", 100000, 0., 100);

  hchHadronMu = CreateH1F("hchHadronMu","", 10000, 0., 100);
  hchHadronEl = CreateH1F("hchHadronEl","", 10000, 0., 100);

  hneuHadronMu = CreateH1F("hneuHadronMu","", 10000, 0., 100);
  hneuHadronEl = CreateH1F("hnueHadronEl","", 10000, 0., 100);

  hphotonMu = CreateH1F("hphotonMu","", 10000, 0., 100);
  hphotonEl = CreateH1F("hphotonEl","", 10000, 0., 100);

  hPUMu = CreateH1F("hPUMu","", 10000, 0., 100);
  hPUEl = CreateH1F("hPUEl","", 10000, 0., 100);
}

// The InsideLoop() function is called for each entry in the tree to be processed  
void test::InsideLoop() {

  TString TheSample = GetInputParameters()->TheNamedString("theSample");
  
  Double_t efficiencyW = effW * triggW;
  Double_t totalW      = -999;
  
  efficiencyW = puW * effW * triggW ;
  totalW      = (1 + 0.5 * (dataset >= 82 && dataset <= 84)) * baseW * efficiencyW * Luminosity;
  
  h_n_PV -> Fill(1,efficiencyW);  
  
  Int_t dphiv = (njet <= 1 || (njet > 1 && dphilljetjet < 165.*TMath::DegToRad()));
  
  Float_t jetbin = njet;
  
  Float_t dphimin = (min(dphilmet1,dphilmet2));
  Float_t fullpmet = 0;
  Float_t trkpmet = 0;
  
  if (dphimin < TMath::Pi() / 2)
    fullpmet = pfType1Met * sin(dphimin);
  else 
    fullpmet = pfType1Met;

  if (dphimin < TMath::Pi() / 2)
    trkpmet = trkMet * sin(dphimin);
  else
    trkpmet = trkMet;
  
   Float_t mpmet = min(trkpmet,fullpmet);

   Float_t metvar = (njet <= 1) ? mpmet : pfType1Met;
  
   Float_t Ht = std_vector_lepton_pt->at(0) + std_vector_lepton_pt->at(1) + pfType1Met;
   
   if(njet > 10) njet = 10;
   for (int i = 0; i < njet; ++i)
     if(std_vector_jet_pt->at(i) > 0)
       Ht += std_vector_jet_pt->at(i);

   //Defining Isolation
   //--------------------------------------------------------------------------
   /*
   Float_t isoOne = 1.;
   Float_t isoZero = 1.;
   
   if(std_vector_lepton_pt->size() > 1){
     if(std_vector_lepton_neutralHadronIso->at(0) + std_vector_lepton_photonIso->at(0) - 0.5 * std_vector_lepton_sumPUPt->at(0) > 0)
       isoZero = 
	 std_vector_lepton_chargedHadronIso->at(0) + 
	 std_vector_lepton_neutralHadronIso->at(0) + 
	 std_vector_lepton_photonIso->at(0) - 
	 0.5 * std_vector_lepton_sumPUPt->at(0);
     else
       isoZero = std_vector_lepton_chargedHadronIso->at(0);
     isoZero = isoZero / std_vector_lepton_pt->at(0);
     cout<<"isoZero: "<<isoZero<<endl;
     
     if(std_vector_lepton_neutralHadronIso->at(1) + std_vector_lepton_photonIso->at(1) - 0.5 * std_vector_lepton_sumPUPt->at(1) > 0)
       isoOne = 
	 std_vector_lepton_chargedHadronIso->at(1) + 
	 std_vector_lepton_neutralHadronIso->at(1) +
	 std_vector_lepton_photonIso->at(1) -
	 0.5 * std_vector_lepton_sumPUPt->at(1);
     else
       isoOne = std_vector_lepton_chargedHadronIso->at(1);
     isoOne = isoOne / std_vector_lepton_pt->at(1);
     cout<<"isoOne :"<<isoOne<<endl;
     cout<<" "<<endl;
   }
   
   if(fabs(std_vector_lepton_id->at(0)) == 13){
     //hIsoMu->Fill(isoZero);
     hchHadronMu->Fill(std_vector_lepton_chargedHadronIso->at(0));
     hneuHadronMu->Fill(std_vector_lepton_neutralHadronIso->at(0));
     hphotonMu->Fill(std_vector_lepton_photonIso->at(0));
     hPUMu->Fill(std_vector_lepton_sumPUPt->at(0));
   }

   if(fabs(std_vector_lepton_id->at(1)) == 13){
     hIsoMu->Fill(isoOne);
     hchHadronMu->Fill(std_vector_lepton_chargedHadronIso->at(1));
     hneuHadronMu->Fill(std_vector_lepton_neutralHadronIso->at(1));
     hphotonMu->Fill(std_vector_lepton_photonIso->at(1));
     hPUMu->Fill(std_vector_lepton_sumPUPt->at(1));
   }
   
   if(fabs(std_vector_lepton_id->at(0)) == 11){
     hIsoEl->Fill(isoZero);
     hchHadronEl->Fill(std_vector_lepton_chargedHadronIso->at(0));
     hneuHadronEl->Fill(std_vector_lepton_neutralHadronIso->at(0));
     hphotonEl->Fill(std_vector_lepton_photonIso->at(0));
     hPUEl->Fill(std_vector_lepton_sumPUPt->at(0));
   }

   if(fabs(std_vector_lepton_id->at(1)) == 11){
     hIsoEl->Fill(isoOne);
     hchHadronEl->Fill(std_vector_lepton_chargedHadronIso->at(1));
     hneuHadronEl->Fill(std_vector_lepton_neutralHadronIso->at(1));
     hphotonEl->Fill(std_vector_lepton_photonIso->at(1));
     hPUEl->Fill(std_vector_lepton_sumPUPt->at(1));
   }
   */
   /*
   //------------------------------------------------------------------------------
   //reject background prompt-prompt events
   Float_t prompt = 1;
   if(TheSample == "TTJets" || TheSample == "Top" || TheSample == "QCD" || TheSample == "WJets")
   if(fabs(std_vector_leptonGen_mpid -> at(0)) == 24 && 
      fabs(std_vector_leptonGen_mstatus -> at(0)) > 20 && 
      fabs(std_vector_leptonGen_mstatus -> at(0)) < 30)
	   if(fabs(std_vector_leptonGen_mpid -> at(1)) == 24 && 
	      fabs(std_vector_leptonGen_mstatus -> at(1)) > 20 && 
	      fabs(std_vector_leptonGen_mstatus -> at(1)) < 30)
		   prompt = 0;
   */
		   
   // The selection begins here
   //--------------------------------------------------------------------------
   if (std_vector_lepton_pt->at(0) > 20)
     if (std_vector_lepton_pt->at(1) > 20) 
       if (ch1*ch2 > 0)
	 if (IsTightLepton(0,
			   std_vector_lepton_id,
			   std_vector_lepton_isTightMuon,
			   std_vector_electron_deltaEtaIn,
			   std_vector_electron_sigmaIetaIeta,
			   std_vector_electron_HoE,
			   std_vector_electron_d0,
			   std_vector_electron_dz,
			   std_vector_electron_ooEooP,
			   std_vector_electron_passConversion,
			   std_vector_electron_scEta,
			   std_vector_electron_deltaPhiIn
			   ))  
	   if(IsIsolatedLepton(0,
			       std_vector_lepton_chargedHadronIso,
			       std_vector_lepton_neutralHadronIso,
			       std_vector_lepton_sumPUPt,
			       std_vector_lepton_pt,
			       std_vector_lepton_id,
			       std_vector_lepton_eta,
			       std_vector_lepton_photonIso,
			       jetRho
			       ))

	     if (IsTightLepton(1,
			       std_vector_lepton_id,
			       std_vector_lepton_isTightMuon,
			       std_vector_electron_deltaEtaIn,
			       std_vector_electron_sigmaIetaIeta,
			       std_vector_electron_HoE,
			       std_vector_electron_d0,
			       std_vector_electron_dz,
			       std_vector_electron_ooEooP,
			       std_vector_electron_passConversion,
			       std_vector_electron_scEta,
			       std_vector_electron_deltaPhiIn
			       ))  
	       if(IsIsolatedLepton(1,
				   std_vector_lepton_chargedHadronIso,
				   std_vector_lepton_neutralHadronIso,
				   std_vector_lepton_sumPUPt,
				   std_vector_lepton_pt,
				   std_vector_lepton_id,
				   std_vector_lepton_eta,
				   std_vector_lepton_photonIso,
				   jetRho
				   ))
	 //if (std_vector_lepton_BestTrackdxy -> at(0) < 0.01 && std_vector_lepton_BestTrackdz -> at(0) < 0.1)
	 //if (std_vector_lepton_BestTrackdxy -> at(1) < 0.01 && std_vector_lepton_BestTrackdz -> at(1) < 0.1)
	 //if (prompt == 1)
	 //if (std_vector_lepton_isTightMuon->at(0) == 1)// && std_vector_lepton_isTightMuon->at(1) == 0 || 
	 //if (std_vector_lepton_isTightMuon->at(1) == 1)// && std_vector_lepton_isTightMuon->at(1) == 1)
	 //if (isoZero < 0.12)
	 //if (isoOne < 0.12)
		       if ( (SelectedChannel == -1)                                   || 
			    (channel == SelectedChannel)                              || 
			    (SelectedChannel == 4 && (channel == 2 || channel == 3) ) || 
			    (SelectedChannel == 5 && (channel == 0 || channel == 1) ) 
			    )
			 {
			   hNJets30TwoLeptonsLevel        ->Fill(njet,      totalW);
			     
		  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		  //
		  // Main analisis
		  //
		  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		  
		  hWTrigger   ->Fill(1, totalW); 
		  hWeffTrigger->Fill(1, efficiencyW);
		  
		  hPtLepton1TwoLeptonsLevel      ->Fill(pt1,       totalW);
		  hPtLepton2TwoLeptonsLevel      ->Fill(pt2,       totalW);
		  hPtDiLeptonTwoLeptonsLevel     ->Fill(ptll,      totalW);
		  hMinvTwoLeptonsLevel           ->Fill(mll,       totalW);
		  hMtTwoLeptonsLevel             ->Fill(mth,       totalW);

		  hpfMetTwoLeptonsLevel          ->Fill(pfType1Met,totalW);
		  //hppfMetTwoLeptonsLevel         ->Fill(ppfmet,    totalW);
		  hchMetTwoLeptonsLevel          ->Fill(chmet,     totalW);
		  hpchMetTwoLeptonsLevel         ->Fill(pchmet,    totalW);
		  hpminMetTwoLeptonsLevel        ->Fill(mpmet,     totalW);
		  hDeltaRLeptonsTwoLeptonsLevel  ->Fill(drll,      totalW);
		  hDeltaPhiLeptonsTwoLeptonsLevel->Fill(dphill,    totalW);
		  hDPhiPtllJetTwoLeptonsLevel    ->Fill(dphilljet, totalW);
		  
		  
		  h_TwoLeptons_TightFailEvents ->Fill(1, totalW); 
		  
		  if (nextra == 0) {
		    
		    hWExtraLepton->Fill(1, totalW);
		    hWeffExtraLepton->Fill(1, efficiencyW);
		    
		    if (pfType1Met > 20 ) { // removed for differential xsec
		      
		      hWMetCut->Fill(1, totalW);
		      hWeffMetCut->Fill(1, efficiencyW);
		      
		      if (mll > 12) {
			
			hWLowMinv->Fill(1, totalW);
			hWeffLowMinv->Fill(1, efficiencyW);

			if (mpmet > 20){

			  if (dphiv || channel == 2 || channel == 3) {
			      
			  hWDeltaPhiJet->Fill(1, totalW);
			  hWeffDeltaPhiJet->Fill(1, efficiencyW);
			 
			  if ( ptll>30 /*&& (!sameflav || ptll>45)*/ ) {
			    
			    hWPtll->Fill(1, totalW);			    
			    hWeffPtll->Fill(1, efficiencyW);			    

			    hWnJets->Fill(njet, totalW);
			    hWeffnJets->Fill(njet, efficiencyW);
			    
			    hWnBtaggedJets->Fill(nbjet, totalW);
			    hWeffnBtaggedJets->Fill(nbjet, efficiencyW);
			    
			    /*if (bveto_mu)*/ {

			    for (Int_t jetNumber = 0; jetNumber < 3 ; ++jetNumber){
			      if (jetbin >= 3) jetbin = 2;
			      if(jetNumber == jetbin){
				hHt[jetNumber]->Fill(Ht,totalW);				    
			      }
			    }
			  
			  hHt[3]->Fill(Ht,totalW);				    
			  /*
			  hWTopTagging->Fill(1, totalW);
			  hWeffTopTagging->Fill(1, efficiencyW);
			  
			  hWSoftMuVeto->Fill(1, totalW);
			  hWeffSoftMuVeto->Fill(1, efficiencyW);
			  */
			  hHtAfter[3]->Fill(Ht,totalW);				    
			  
			  hPtLepton1WWLevelNoHt[3]      ->Fill(pt1,       totalW);
			  hPtLepton2WWLevelNoHt[3]      ->Fill(pt2,       totalW);
			  hPtDiLeptonWWLevelNoHt[3]     ->Fill(ptll,      totalW);
			  hMinvWWLevelNoHt[3]           ->Fill(mll,       totalW);
			  hMtWWLevelNoHt[3]             ->Fill(mth,       totalW);
			  hNJets30WWLevelNoHt[3]        ->Fill(jetbin,    totalW);
			  hpfMetWWLevelNoHt[3]          ->Fill(pfType1Met,totalW);
			  //hppfMetWWLevelNoHt[3]         ->Fill(ppfmet,    totalW);
			  hchMetWWLevelNoHt[3]          ->Fill(chmet,     totalW);
			  hpchMetWWLevelNoHt[3]         ->Fill(pchmet,    totalW);
			  hpminMetWWLevelNoHt[3]        ->Fill(mpmet,     totalW);
			  hDeltaRLeptonsWWLevelNoHt[3]  ->Fill(drll,      totalW);
			  hDeltaPhiLeptonsWWLevelNoHt[3]->Fill(dphill,    totalW);
			  hDPhiPtllJetWWLevelNoHt[3]    ->Fill(dphilljet, totalW);
			  hWnJetsBveto                  ->Fill(njet,      totalW);
			  hWeffnJetsBveto               ->Fill(njet, efficiencyW);
			  hEffNoHt[3]                   ->Fill(1,    efficiencyW);
			  hSigMuNoHt[3]                 ->Fill(std_vector_lepton_muSIP3D->at(0),totalW);
			  hSigElNoHt[3]                 ->Fill(std_vector_lepton_elSIP3D->at(0),totalW);
			  
			  //bveto Ht 
			  if(Ht < 250){
			    
			    hPtLepton1WWLevel[3]      ->Fill(pt1,       totalW);
			    hPtLepton2WWLevel[3]      ->Fill(pt2,       totalW);
			    hPtDiLeptonWWLevel[3]     ->Fill(ptll,      totalW);
			    hMinvWWLevel[3]           ->Fill(mll,       totalW);
			    hMtWWLevel[3]             ->Fill(mth,       totalW);
			    hNJets30WWLevel[3]        ->Fill(jetbin,    totalW);
			    hpfMetWWLevel[3]          ->Fill(pfType1Met,totalW);
			    //hppfMetWWLevel[3]         ->Fill(ppfmet,    totalW);
			    hchMetWWLevel[3]          ->Fill(chmet,     totalW);
			    hpchMetWWLevel[3]         ->Fill(pchmet,    totalW);
			    hpminMetWWLevel[3]        ->Fill(mpmet,     totalW);
			    hDeltaRLeptonsWWLevel[3]  ->Fill(drll,      totalW);
			    hDeltaPhiLeptonsWWLevel[3]->Fill(dphill,    totalW);
			    hDPhiPtllJetWWLevel[3]    ->Fill(dphilljet, totalW);
			    hWnJetsBvetoAfterHt       ->Fill(njet,      totalW);
			    hWeffnJetsBvetoAfterHt    ->Fill(njet, efficiencyW);					
			    hEff[3]                   ->Fill(1,    efficiencyW);					
			    hSigMu[3]                 ->Fill(std_vector_lepton_muSIP3D->at(0),totalW);
			    hSigEl[3]                 ->Fill(std_vector_lepton_elSIP3D->at(0),totalW);
			  }
			  
			  //bveto Ht 
			  if(Ht > 250){
			    
			    hPtLepton1WWLevelHtPlus[3]      ->Fill(pt1,       totalW);
			    hPtLepton2WWLevelHtPlus[3]      ->Fill(pt2,       totalW);
			    hPtDiLeptonWWLevelHtPlus[3]     ->Fill(ptll,      totalW);
			    hMinvWWLevelHtPlus[3]           ->Fill(mll,       totalW);
			    hMtWWLevelHtPlus[3]             ->Fill(mth,       totalW);
			    hNJets30WWLevelHtPlus[3]        ->Fill(jetbin,    totalW);
			    hpfMetWWLevelHtPlus[3]          ->Fill(pfType1Met,totalW);
			    //hppfMetWWLevelHtPlus[3]         ->Fill(ppfmet,    totalW);
			    hchMetWWLevelHtPlus[3]          ->Fill(chmet,     totalW);
			    hpchMetWWLevelHtPlus[3]         ->Fill(pchmet,    totalW);
			    hpminMetWWLevelHtPlus[3]        ->Fill(mpmet,     totalW);
			    hDeltaRLeptonsWWLevelHtPlus[3]  ->Fill(drll,      totalW);
			    hDeltaPhiLeptonsWWLevelHtPlus[3]->Fill(dphill,    totalW);
			    hDPhiPtllJetWWLevelHtPlus[3]    ->Fill(dphilljet, totalW);
			    hEffHtPlus[3]                   ->Fill(1,    efficiencyW);					
			    hSigMuHtPlus[3]                 ->Fill(std_vector_lepton_muSIP3D->at(0),totalW);
			    hSigElHtPlus[3]                 ->Fill(std_vector_lepton_elSIP3D->at(0),totalW);
			  }
			  
			  for (Int_t jetNumber = 0; jetNumber < 3 ; ++jetNumber){
			    if (jetbin >= 3) jetbin = 2;
			    if(jetNumber == jetbin){
			      
			      hHtAfter[jetNumber]                   ->Fill(Ht,        totalW);				    
			      hPtLepton1WWLevelNoHt[jetNumber]      ->Fill(pt1,       totalW);
			      hPtLepton2WWLevelNoHt[jetNumber]      ->Fill(pt2,       totalW);
			      hPtDiLeptonWWLevelNoHt[jetNumber]     ->Fill(ptll,      totalW);
			      hMinvWWLevelNoHt[jetNumber]           ->Fill(mll,       totalW);
			      hMtWWLevelNoHt[jetNumber]             ->Fill(mth,       totalW);
			      hNJets30WWLevelNoHt[jetNumber]        ->Fill(jetNumber, totalW);
			      hpfMetWWLevelNoHt[jetNumber]          ->Fill(pfType1Met,totalW);
			      //hppfMetWWLevelNoHt[jetNumber]         ->Fill(ppfmet,    totalW);
			      hchMetWWLevelNoHt[jetNumber]          ->Fill(chmet,     totalW);
			      hpchMetWWLevelNoHt[jetNumber]         ->Fill(pchmet,    totalW);
			      hpminMetWWLevelNoHt[jetNumber]        ->Fill(mpmet,     totalW);
			      hDeltaRLeptonsWWLevelNoHt[jetNumber]  ->Fill(drll,      totalW);
			      hDeltaPhiLeptonsWWLevelNoHt[jetNumber]->Fill(dphill,    totalW);
			      hDPhiPtllJetWWLevelNoHt[jetNumber]    ->Fill(dphilljet, totalW);
			      hEffNoHt[jetNumber]                   ->Fill(1,    efficiencyW);
			      hSigMuNoHt[jetNumber]                 ->Fill(std_vector_lepton_muSIP3D->at(0),totalW);
			      hSigElNoHt[jetNumber]                 ->Fill(std_vector_lepton_elSIP3D->at(0),totalW);
			      
			      //bveto Ht  
			      if(Ht < 250){
				
				hPtLepton1WWLevel[jetNumber]      ->Fill(pt1,       totalW);
				hPtLepton2WWLevel[jetNumber]      ->Fill(pt2,       totalW);
				hPtDiLeptonWWLevel[jetNumber]     ->Fill(ptll,      totalW);
				hMinvWWLevel[jetNumber]           ->Fill(mll,       totalW);
				hMtWWLevel[jetNumber]             ->Fill(mth,       totalW);
				hNJets30WWLevel[jetNumber]        ->Fill(jetNumber, totalW);
				hpfMetWWLevel[jetNumber]          ->Fill(pfType1Met,totalW);
				//hppfMetWWLevel[jetNumber]         ->Fill(ppfmet,    totalW);
				hchMetWWLevel[jetNumber]          ->Fill(chmet,     totalW);
				hpchMetWWLevel[jetNumber]         ->Fill(pchmet,    totalW);
				hpminMetWWLevel[jetNumber]        ->Fill(mpmet,     totalW);
				hDeltaRLeptonsWWLevel[jetNumber]  ->Fill(drll,      totalW);
				hDeltaPhiLeptonsWWLevel[jetNumber]->Fill(dphill,    totalW);
				hDPhiPtllJetWWLevel[jetNumber]    ->Fill(dphilljet, totalW);
				hEff[jetNumber]                   ->Fill(1,    efficiencyW);
				hSigMu[jetNumber]                 ->Fill(std_vector_lepton_muSIP3D->at(0),totalW);
				hSigEl[jetNumber]                 ->Fill(std_vector_lepton_elSIP3D->at(0),totalW);
			      }
			      
			      //bveto Ht  
			      if(Ht > 250){
				
				hPtLepton1WWLevelHtPlus[jetNumber]      ->Fill(pt1,       totalW);
				hPtLepton2WWLevelHtPlus[jetNumber]      ->Fill(pt2,       totalW);
				hPtDiLeptonWWLevelHtPlus[jetNumber]     ->Fill(ptll,      totalW);
				hMinvWWLevelHtPlus[jetNumber]           ->Fill(mll,       totalW);
				hMtWWLevelHtPlus[jetNumber]             ->Fill(mth,       totalW);
				hNJets30WWLevelHtPlus[jetNumber]        ->Fill(jetNumber, totalW);
				hpfMetWWLevelHtPlus[jetNumber]          ->Fill(pfType1Met,totalW);
				//hppfMetWWLevelHtPlus[jetNumber]         ->Fill(ppfmet,    totalW);
				hchMetWWLevelHtPlus[jetNumber]          ->Fill(chmet,     totalW);
				hpchMetWWLevelHtPlus[jetNumber]         ->Fill(pchmet,    totalW);
				hpminMetWWLevelHtPlus[jetNumber]        ->Fill(mpmet,     totalW);
				hDeltaRLeptonsWWLevelHtPlus[jetNumber]  ->Fill(drll,      totalW);
				hDeltaPhiLeptonsWWLevelHtPlus[jetNumber]->Fill(dphill,    totalW);
				hDPhiPtllJetWWLevelHtPlus[jetNumber]    ->Fill(dphilljet, totalW);
				hEffHtPlus[jetNumber]                   ->Fill(1,    efficiencyW);
				hSigMuHtPlus[jetNumber]                 ->Fill(std_vector_lepton_muSIP3D->at(0),totalW);
				hSigElHtPlus[jetNumber]                 ->Fill(std_vector_lepton_elSIP3D->at(0),totalW);
			      }
			    }  					
			  }
			  h_WWLevel_TightFailEvents ->Fill(1, totalW); 
			    }
			  }
			}
			}
		      }
		    }
		  }
			 }		       
		  // Define Normalization Factor for MC samples 
		  
		  //------------------------------------------------------------------------------
		  // Define weights
		  //------------------------------------------------------------------------------
		  
		  float pileupweight = 1;
		  
		  //  if (!IsDATA)
		  // pileupweight = fPUWeight->GetWeight(T_Event_nPU);
		  
		  
		  
		  
		  double factN = 1;
		  
		  if (XSection > 0) factN = XSection * Luminosity / NEvents;
		  
		  
		  //factN = factN*pileupweight;
		  
		  
		  //------------------------------------------------------------------------------
		  // Init variables
		  //------------------------------------------------------------------------------
		  
		  
			   
		       
}// end inside Loop





void test::SetDataMembersAtTermination() {

  h_n_PV = ((TH1F*) FindOutput("h_n_PV")); 
  
  // Counting histograms                             
  //----------------------------------------------------------------------------       
  
  hWTrigger     = ((TH1F*) FindOutput("hWTrigger"));
  hWMetCut      = ((TH1F*) FindOutput("hWMetCut"));
  hWLowMinv     = ((TH1F*) FindOutput("hWLowMinv"));
  hWZVeto       = ((TH1F*) FindOutput("hWZVeto"));
  hWpMetCut     = ((TH1F*) FindOutput("hWpMetCut"));
  hWJetVeto     = ((TH1F*) FindOutput("hWJetVeto"));
  hWnJets       = ((TH1F*) FindOutput("hWnJets"));
  hWeffnJets    = ((TH1F*) FindOutput("hWeffnJets"));
  hWnBtaggedJets = ((TH1F*) FindOutput("hWnBtaggedJets"));
  hWeffnBtaggedJets = ((TH1F*) FindOutput("hWeffnBtaggedJets"));
  hWnJetsBveto      = ((TH1F*) FindOutput("hWnJetsBveto"));
  hWeffnJetsBveto   = ((TH1F*) FindOutput("hWeffnJetsBveto"));
  hWnJetsBvetoAfterHt       = ((TH1F*) FindOutput("hWnJetsBvetoAfterHt"));
  hWeffnJetsBvetoAfterHt       = ((TH1F*) FindOutput("hWeffnJetsBvetoAfterHt"));

  hWDeltaPhiJet = ((TH1F*) FindOutput("hWDeltaPhiJet"));
  hWSoftMuVeto  = ((TH1F*) FindOutput("hWSoftMuVeto"));
  hWExtraLepton = ((TH1F*) FindOutput("hWExtraLepton"));
  hWPtll        = ((TH1F*) FindOutput("hWPtll"));
  hWTopTagging  = ((TH1F*) FindOutput("hWTopTagging"));

  hWeffTrigger     = ((TH1F*) FindOutput("hWeffTrigger"));
  hWeffMetCut      = ((TH1F*) FindOutput("hWeffMetCut"));
  hWeffLowMinv     = ((TH1F*) FindOutput("hWeffLowMinv"));
  hWeffZVeto       = ((TH1F*) FindOutput("hWeffZVeto"));
  hWeffpMetCut     = ((TH1F*) FindOutput("hWeffpMetCut"));
  hWeffJetVeto     = ((TH1F*) FindOutput("hWeffJetVeto"));
  hWeffDeltaPhiJet = ((TH1F*) FindOutput("hWeffDeltaPhiJet"));
  hWeffSoftMuVeto  = ((TH1F*) FindOutput("hWeffSoftMuVeto"));
  hWeffExtraLepton = ((TH1F*) FindOutput("hWeffExtraLepton"));
  hWeffPtll        = ((TH1F*) FindOutput("hWeffPtll"));
  hWeffTopTagging  = ((TH1F*) FindOutput("hWeffTopTagging"));

  

  // WW level histograms                                                                                                                                                           
  //----------------------------------------------------------------------------                                                                                                   
  for (Int_t qq = 0; qq < 4; ++qq){
  hPtLepton1WWLevel[qq]       = ((TH1F*) FindOutput("hPtLepton1WWLevel%.1i",qq));
  hPtLepton2WWLevel[qq]       = ((TH1F*) FindOutput("hPtLepton2WWLevel%.1i",qq));
  hPtDiLeptonWWLevel[qq]      = ((TH1F*) FindOutput("hPtDiLeptonWWLevel%.1i",qq));
  hMinvWWLevel[qq]            = ((TH1F*) FindOutput("hMinvWWLevel%.1i",qq));     // %.1i",qq 
  hMtWWLevel[qq]              = ((TH1F*) FindOutput("hMtWWLevel%.1i",qq));
  hNJets30WWLevel[qq]         = ((TH1F*) FindOutput("hNJetsPF30WWLevel%.1i",qq));
  hpfMetWWLevel[qq]           = ((TH1F*) FindOutput("hpfMetWWLevel%.1i",qq));
  hppfMetWWLevel[qq]          = ((TH1F*) FindOutput("hppfMetWWLevel%.1i",qq));
  hchMetWWLevel[qq]           = ((TH1F*) FindOutput("hchMetWWLevel%.1i",qq));
  hpchMetWWLevel[qq]          = ((TH1F*) FindOutput("hpchMetWWLevel%.1i",qq));
  hpminMetWWLevel[qq]         = ((TH1F*) FindOutput("hpminMetWWLevel%.1i",qq));
  hDeltaRLeptonsWWLevel[qq]   = ((TH1F*) FindOutput("hDeltaRLeptonsWWLevel%.1i",qq));
  hDeltaPhiLeptonsWWLevel[qq] = ((TH1F*) FindOutput("hDeltaPhiLeptonsWWLevel%.1i",qq));
  hDPhiPtllJetWWLevel[qq]     = ((TH1F*) FindOutput("hDPhiPtllJetWWLevel%.1i",qq));
  hDataEvents[qq]             = ((TH1F*) FindOutput("hDataEvents%.1i",qq));
  hBackgroundEvents[qq]       = ((TH1F*) FindOutput("hBackgroundEvents%.1i",qq));
  hEff[qq]                    = ((TH1F*) FindOutput("hEff%.1i",qq));
  hSigEl[qq]                  = ((TH1F*) FindOutput("hSigEl%.1i",qq));
  hSigMu[qq]                  = ((TH1F*) FindOutput("hSigMu%.1i",qq));

  hPtLepton1WWLevelNoHt[qq]       = ((TH1F*) FindOutput("hPtLepton1WWLevelNoHt%.1i",qq));
  hPtLepton2WWLevelNoHt[qq]       = ((TH1F*) FindOutput("hPtLepton2WWLevelNoHt%.1i",qq));
  hPtDiLeptonWWLevelNoHt[qq]      = ((TH1F*) FindOutput("hPtDiLeptonWWLevelNoHt%.1i",qq));
  hMinvWWLevelNoHt[qq]            = ((TH1F*) FindOutput("hMinvWWLevelNoHt%.1i",qq));     // %.1i",qq 
  hMtWWLevelNoHt[qq]              = ((TH1F*) FindOutput("hMtWWLevelNoHt%.1i",qq));
  hNJets30WWLevelNoHt[qq]         = ((TH1F*) FindOutput("hNJetsPF30WWLevelNoHt%.1i",qq));
  hpfMetWWLevelNoHt[qq]           = ((TH1F*) FindOutput("hpfMetWWLevelNoHt%.1i",qq));
  hppfMetWWLevelNoHt[qq]          = ((TH1F*) FindOutput("hppfMetWWLevelNoHt%.1i",qq));
  hchMetWWLevelNoHt[qq]           = ((TH1F*) FindOutput("hchMetWWLevelNoHt%.1i",qq));
  hpchMetWWLevelNoHt[qq]          = ((TH1F*) FindOutput("hpchMetWWLevelNoHt%.1i",qq));
  hpminMetWWLevelNoHt[qq]         = ((TH1F*) FindOutput("hpminMetWWLevelNoHt%.1i",qq));
  hDeltaRLeptonsWWLevelNoHt[qq]   = ((TH1F*) FindOutput("hDeltaRLeptonsWWLevelNoHt%.1i",qq));
  hDeltaPhiLeptonsWWLevelNoHt[qq] = ((TH1F*) FindOutput("hDeltaPhiLeptonsWWLevelNoHt%.1i",qq));
  hDPhiPtllJetWWLevelNoHt[qq]     = ((TH1F*) FindOutput("hDPhiPtllJetWWLevelNoHt%.1i",qq));
  hDataEventsNoHt[qq]             = ((TH1F*) FindOutput("hDataEventsNoHt%.1i",qq));
  hBackgroundEventsNoHt[qq]       = ((TH1F*) FindOutput("hBackgroundEventsNoHt%.1i",qq));
  hEffNoHt[qq]                    = ((TH1F*) FindOutput("hEffNoHt%.1i",qq));
  hSigElNoHt[qq]                  = ((TH1F*) FindOutput("hSigElNoHt%.1i",qq));
  hSigMuNoHt[qq]                  = ((TH1F*) FindOutput("hSigMuNoHt%.1i",qq));

  hPtLepton1WWLevelHtPlus[qq]       = ((TH1F*) FindOutput("hPtLepton1WWLevelHtPlus%.1i",qq));
  hPtLepton2WWLevelHtPlus[qq]       = ((TH1F*) FindOutput("hPtLepton2WWLevelHtPlus%.1i",qq));
  hPtDiLeptonWWLevelHtPlus[qq]      = ((TH1F*) FindOutput("hPtDiLeptonWWLevelHtPlus%.1i",qq));
  hMinvWWLevelHtPlus[qq]            = ((TH1F*) FindOutput("hMinvWWLevelHtPlus%.1i",qq));     // HtPlus%.1i",qq 
  hMtWWLevelHtPlus[qq]              = ((TH1F*) FindOutput("hMtWWLevelHtPlus%.1i",qq));
  hNJets30WWLevelHtPlus[qq]         = ((TH1F*) FindOutput("hNJetsPF30WWLevelHtPlus%.1i",qq));
  hpfMetWWLevelHtPlus[qq]           = ((TH1F*) FindOutput("hpfMetWWLevelHtPlus%.1i",qq));
  hppfMetWWLevelHtPlus[qq]          = ((TH1F*) FindOutput("hppfMetWWLevelHtPlus%.1i",qq));
  hchMetWWLevelHtPlus[qq]           = ((TH1F*) FindOutput("hchMetWWLevelHtPlus%.1i",qq));
  hpchMetWWLevelHtPlus[qq]          = ((TH1F*) FindOutput("hpchMetWWLevelHtPlus%.1i",qq));
  hpminMetWWLevelHtPlus[qq]         = ((TH1F*) FindOutput("hpminMetWWLevelHtPlus%.1i",qq));
  hDeltaRLeptonsWWLevelHtPlus[qq]   = ((TH1F*) FindOutput("hDeltaRLeptonsWWLevelHtPlus%.1i",qq));
  hDeltaPhiLeptonsWWLevelHtPlus[qq] = ((TH1F*) FindOutput("hDeltaPhiLeptonsWWLevelHtPlus%.1i",qq));
  hDPhiPtllJetWWLevelHtPlus[qq]     = ((TH1F*) FindOutput("hDPhiPtllJetWWLevelHtPlus%.1i",qq));
  hDataEventsHtPlus[qq]             = ((TH1F*) FindOutput("hDataEventsHtPlus%.1i",qq));
  hBackgroundEventsHtPlus[qq]       = ((TH1F*) FindOutput("hBackgroundEventsHtPlus%.1i",qq));
  hEffHtPlus[qq]                    = ((TH1F*) FindOutput("hEffHtPlus%.1i",qq));
  hSigElHtPlus[qq]                  = ((TH1F*) FindOutput("hSigElHtPlus%.1i",qq));
  hSigMuHtPlus[qq]                  = ((TH1F*) FindOutput("hSigMuHtPlus%.1i",qq));

  hHt[qq]                     = ((TH1F*) FindOutput("hHt%.1i",qq));
  hHtAfter[qq]                = ((TH1F*) FindOutput("hHtAfter%.1i",qq));
 }

  h_WWLevel_TightFailEvents = ((TH1F*) FindOutput("h_WWLevel_TightFailEvents"));
  h_WWLevel_TightTightEvents = ((TH1F*) FindOutput("h_WWLevel_TightTightEvents"));
  h_WWLevel_TightLooseEvents = ((TH1F*) FindOutput("h_WWLevel_TightLooseEvents"));


  // TwoLeptons level histograms                                                                                                                                                   
  //----------------------------------------------------------------------------                                                                                                   
  hPtLepton1TwoLeptonsLevel       = ((TH1F*) FindOutput("hPtLepton1TwoLeptonsLevel"));
  hPtLepton2TwoLeptonsLevel       = ((TH1F*) FindOutput("hPtLepton2TwoLeptonsLevel"));
  hPtDiLeptonTwoLeptonsLevel      = ((TH1F*) FindOutput("hPtDiLeptonTwoLeptonsLevel"));
  hMinvTwoLeptonsLevel            = ((TH1F*) FindOutput("hMinvTwoLeptonsLevel"));
  hMtTwoLeptonsLevel              = ((TH1F*) FindOutput("hMtTwoLeptonsLevel"));
  hNJets30TwoLeptonsLevel         = ((TH1F*) FindOutput("hNJets30TwoLeptonsLevel"));
  hpfMetTwoLeptonsLevel           = ((TH1F*) FindOutput("hpfMetTwoLeptonsLevel"));
  hppfMetTwoLeptonsLevel          = ((TH1F*) FindOutput("hppfMetTwoLeptonsLevel"));
  hchMetTwoLeptonsLevel           = ((TH1F*) FindOutput("hchMetTwoLeptonsLevel"));
  hpchMetTwoLeptonsLevel          = ((TH1F*) FindOutput("hpchMetTwoLeptonsLevel"));
  hpminMetTwoLeptonsLevel         = ((TH1F*) FindOutput("hpminMetTwoLeptonsLevel"));
  hDeltaRLeptonsTwoLeptonsLevel   = ((TH1F*) FindOutput("hDeltaRLeptonsTwoLeptonsLevel"));
  hDeltaPhiLeptonsTwoLeptonsLevel = ((TH1F*) FindOutput("hDeltaPhiLeptonsTwoLeptonsLevel"));
  hDPhiPtllJetTwoLeptonsLevel     = ((TH1F*) FindOutput("hDPhiPtllJetTwoLeptonsLevel"));

  h_TwoLeptons_TightFailEvents  = ((TH1F*) FindOutput("h_TwoLeptons_TightFailEvents"));
  h_TwoLeptons_TightTightEvents = ((TH1F*) FindOutput("h_TwoLeptons_TightTightEvents"));
  h_TwoLeptons_TightLooseEvents = ((TH1F*) FindOutput("h_TwoLeptons_TightLooseEvents"));


  //Isolation Plots
  //-----------------------------------------------------------------------------

  hIsoMu = ((TH1F*) FindOutput("hIsoMu"));
  hIsoEl = ((TH1F*) FindOutput("hIsoEl"));

  hchHadronMu = ((TH1F*) FindOutput("hchHadronMu"));
  hchHadronEl = ((TH1F*) FindOutput("hchHadronEl"));

  hneuHadronMu = ((TH1F*) FindOutput("hneuHadronMu"));
  hneuHadronEl = ((TH1F*) FindOutput("hnueHadronEl"));

  hphotonMu = ((TH1F*) FindOutput("hphotonMu"));
  hphotonEl = ((TH1F*) FindOutput("hphotonEl"));

  hPUMu = ((TH1F*) FindOutput("hPUMu"));
  hPUEl = ((TH1F*) FindOutput("hPUEl"));

}


void test::Summary() {

  cout << " ---------------------------------------------------" << endl;
  cout << " " << endl;
  
  InitialiseParameters();

  cout << " Number of Events::  " << NEvents  << endl;

  double factN = 1.; 
  if (XSection > 0) factN = XSection * Luminosity / NEvents; //fractionoftotalevents;
  
  cout << " Normalization factor: " << factN << endl;
  cout << endl;
  cout << " ---------------------------------------------------" << endl;

}
