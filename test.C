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
  TString sameSign = GetInputParameters()->TheNamedString("SameSign");

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
  hNjetsTwoLeptonsLevel    = CreateH1F("hNjetsTwoLeptonsLevel",         "", 10, 0, 10);
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
    hpfMetWWLevel[nC]           = CreateH1F(Form("hpfMetWWLevel%.1i", nC),           "", 150, 0, 150);
    hpminMetWWLevel[nC]         = CreateH1F(Form("hpminMetWWLevel%.1i", nC),         "", 150, 0, 150);
    hDeltaRLeptonsWWLevel[nC]   = CreateH1F(Form("hDeltaRLeptonsWWLevel%.1i", nC),   "",  50, 0,   5);
    hDeltaPhiLeptonsWWLevel[nC] = CreateH1F(Form("hDeltaPhiLeptonsWWLevel%.1i", nC), "",  32, 0, 3.2);
    hDPhiPtllJetWWLevel[nC]     = CreateH1F(Form("hDPhiPtllJetWWLevel%.1i", nC),     "",  32, 0, 3.2);
    hSigEl[nC]                  = CreateH1F(Form("hSigEl%.1i", nC),                   "", 1000,0.,500);
    hSigMu[nC]                  = CreateH1F(Form("hSigMu%.1i", nC),                   "", 1000,0.,500);

    hPtLepton1WWLevelNoHt[nC]       = CreateH1F(Form("hPtLepton1WWLevelNoHt%.1i", nC),       "", 200, 0, 200);
    hPtLepton2WWLevelNoHt[nC]       = CreateH1F(Form("hPtLepton2WWLevelNoHt%.1i", nC),       "", 200, 0, 200);
    hPtDiLeptonWWLevelNoHt[nC]      = CreateH1F(Form("hPtDiLeptonWWLevelNoHt%.1i", nC),      "", 200, 0, 200);
    hMinvWWLevelNoHt[nC]            = CreateH1F(Form("hMinvWWLevelNoHt%.1i", nC),            "", 200, 0, 200);
    hMtWWLevelNoHt[nC]              = CreateH1F(Form("hMtWWLevelNoHt%.1i", nC),              "", 250, 0, 250);
    hpfMetWWLevelNoHt[nC]           = CreateH1F(Form("hpfMetWWLevelNoHt%.1i", nC),           "", 150, 0, 150);
    hpminMetWWLevelNoHt[nC]         = CreateH1F(Form("hpminMetWWLevelNoHt%.1i", nC),         "", 150, 0, 150);
    hDeltaRLeptonsWWLevelNoHt[nC]   = CreateH1F(Form("hDeltaRLeptonsWWLevelNoHt%.1i", nC),   "",  50, 0,   5);
    hDeltaPhiLeptonsWWLevelNoHt[nC] = CreateH1F(Form("hDeltaPhiLeptonsWWLevelNoHt%.1i", nC), "",  32, 0, 3.2);
    hDPhiPtllJetWWLevelNoHt[nC]     = CreateH1F(Form("hDPhiPtllJetWWLevelNoHt%.1i", nC),     "",  32, 0, 3.2);
    hSigElNoHt[nC]                  = CreateH1F(Form("hSigElNoHt%.1i", nC),                   "", 1000,0.,500);
    hSigMuNoHt[nC]                  = CreateH1F(Form("hSigMuNoHt%.1i", nC),                   "", 1000,0.,500);

    hPtLepton1WWLevelHtPlus[nC]       = CreateH1F(Form("hPtLepton1WWLevelHtPlus%.1i", nC),       "", 200, 0, 200);
    hPtLepton2WWLevelHtPlus[nC]       = CreateH1F(Form("hPtLepton2WWLevelHtPlus%.1i", nC),       "", 200, 0, 200);
    hPtDiLeptonWWLevelHtPlus[nC]      = CreateH1F(Form("hPtDiLeptonWWLevelHtPlus%.1i", nC),      "", 200, 0, 200);
    hMinvWWLevelHtPlus[nC]            = CreateH1F(Form("hMinvWWLevelHtPlus%.1i", nC),            "", 200, 0, 200);
    hMtWWLevelHtPlus[nC]              = CreateH1F(Form("hMtWWLevelHtPlus%.1i", nC),              "", 250, 0, 250);
    hpfMetWWLevelHtPlus[nC]           = CreateH1F(Form("hpfMetWWLevelHtPlus%.1i", nC),           "", 150, 0, 150);
    hpminMetWWLevelHtPlus[nC]         = CreateH1F(Form("hpminMetWWLevelHtPlus%.1i", nC),         "", 150, 0, 150);
    hDeltaRLeptonsWWLevelHtPlus[nC]   = CreateH1F(Form("hDeltaRLeptonsWWLevelHtPlus%.1i", nC),   "",  50, 0,   5);
    hDeltaPhiLeptonsWWLevelHtPlus[nC] = CreateH1F(Form("hDeltaPhiLeptonsWWLevelHtPlus%.1i", nC), "",  32, 0, 3.2);
    hDPhiPtllJetWWLevelHtPlus[nC]     = CreateH1F(Form("hDPhiPtllJetWWLevelHtPlus%.1i", nC),     "",  32, 0, 3.2);
    hSigElHtPlus[nC]                  = CreateH1F(Form("hSigElHtPlus%.1i", nC),                   "", 1000,0.,500);
    hSigMuHtPlus[nC]                  = CreateH1F(Form("hSigMuHtPlus%.1i", nC),                   "", 1000,0.,500);

    hHt[nC]                     = CreateH1F(Form("hHt%.1i",               nC),       "", 3000, 0, 3000);
    hHtAfter[nC]                = CreateH1F(Form("hHtAfter%.1i",          nC),       "", 3000, 0, 3000);

  }

  
  // TwoLeptons level histograms   
  //----------------------------------------------------------------------------
  
  hPtLepton1TwoLeptonsLevel       = CreateH1F("hPtLepton1TwoLeptonsLevel",       "", 200, 0, 200);
  hPtLepton2TwoLeptonsLevel       = CreateH1F("hPtLepton2TwoLeptonsLevel",       "", 200, 0, 200);
  hPtDiLeptonTwoLeptonsLevel      = CreateH1F("hPtDiLeptonTwoLeptonsLevel",      "", 200, 0, 200);
  hMinvTwoLeptonsLevel            = CreateH1F("hMinvTwoLeptonsLevel",            "", 200, 0, 200);
  hMtTwoLeptonsLevel              = CreateH1F("hMtTwoLeptonsLevel",              "", 250, 0, 250);
  hpfMetTwoLeptonsLevel           = CreateH1F("hpfMetTwoLeptonsLevel",           "", 150, 0, 150);
  hpminMetTwoLeptonsLevel         = CreateH1F("hpminMetTwoLeptonsLevel",         "", 150, 0, 150);
  hDeltaRLeptonsTwoLeptonsLevel   = CreateH1F("hDeltaRLeptonsTwoLeptonsLevel",   "",  50, 0,   5);
  hDeltaPhiLeptonsTwoLeptonsLevel = CreateH1F("hDeltaPhiLeptonsTwoLeptonsLevel", "",  32, 0, 3.2);
  hDPhiPtllJetTwoLeptonsLevel     = CreateH1F("hDPhiPtllJetTwoLeptonsLevel",     "",  32, 0, 3.2);
  
}

// The InsideLoop() function is called for each entry in the tree to be processed  
void test::InsideLoop() {

  int lep_size = std_vector_lepton_pt->size();

  //for (int i=0; i<lep_size; i++) printf("testing %f %f\n", Testing(i), std_vector_lepton_pt->at(i)+105.);

  TString TheSample = GetInputParameters()->TheNamedString("theSample");
  TString sameSign  = GetInputParameters()->TheNamedString("SameSign");
  
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

   // The selection begins here
   //--------------------------------------------------------------------------
   if (std_vector_lepton_pt->at(0) > 20)
     if (std_vector_lepton_pt->at(1) > 20) 
       if ((sameSign == "SS" && ch1*ch2 > 0) || (sameSign == "OS" && ch1*ch2 < 0))
	 if (IsTightLepton(0))
	   if (IsIsolatedLepton(0))
	     if (IsTightLepton(1))
	       if (IsIsolatedLepton(1))
		 if ( (SelectedChannel == -1)                                   || 
		      (channel == SelectedChannel)                              || 
		      (SelectedChannel == 4 && (channel == 2 || channel == 3) ) || 
		      (SelectedChannel == 5 && (channel == 0 || channel == 1) ) 
		      )
		   {

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
		     hpminMetTwoLeptonsLevel        ->Fill(mpmet,     totalW);
		     hDeltaRLeptonsTwoLeptonsLevel  ->Fill(drll,      totalW);
		     hDeltaPhiLeptonsTwoLeptonsLevel->Fill(dphill,    totalW);
		     hDPhiPtllJetTwoLeptonsLevel    ->Fill(dphilljet, totalW);
		     hNjetsTwoLeptonsLevel          ->Fill(njet,      totalW);
		     
		     if (nextra == 0) {
		       
		       hWExtraLepton->Fill(1, totalW);
		       hWeffExtraLepton->Fill(1, efficiencyW);
		       
		       if (pfType1Met > 20 ) { // removed for differential xsec
			 
			 hWMetCut->Fill(1, totalW);
			 hWeffMetCut->Fill(1, efficiencyW);
			 
			 if (mll > 12) {
			   
			   hWLowMinv->Fill(1, totalW);
			   hWeffLowMinv->Fill(1, efficiencyW);
			   
			   //zveto (in case of same flavour)
			   if ( fabs(ZMASS - mll) > 15 || 
				channel == 0           ||
				channel == 1           ){
			     
			     hWZVeto->Fill(1, totalW);
			     hWeffZVeto->Fill(1, efficiencyW);

			     if (mpmet > 20){

			       hWpMetCut->Fill(1, totalW);
			       hWeffpMetCut->Fill(1, efficiencyW);
			       
			       if (dphiv || channel == 2 || channel == 3) {
				 
				 hWDeltaPhiJet->Fill(1, totalW);
				 hWeffDeltaPhiJet->Fill(1, efficiencyW);
				 
				 if ( ptll>30 && (channel == 2 || channel == 3 || ptll>45) ) {
				   
				   hWPtll->Fill(1, totalW);			    
				   hWeffPtll->Fill(1, efficiencyW);			    
				   
				   hWnJets->Fill(njet, totalW);
				   hWeffnJets->Fill(njet, efficiencyW);
				   
				   hWnBtaggedJets->Fill(nbjet, totalW);
				   hWeffnBtaggedJets->Fill(nbjet, efficiencyW);
				   
				   hHt[3]->Fill(Ht,totalW);				    
				   
				   for (Int_t jetNumber = 0; jetNumber < 3 ; ++jetNumber){
				     if (jetbin >= 3) jetbin = 2;
				     if(jetNumber == jetbin){
				       hHt[jetNumber]->Fill(Ht,totalW);				    
				     }
				   }
				   
				   //b-veto
				   if (bveto_ip == 1 && nbjettche == 0) {
				     
				     hWTopTagging->Fill(1, totalW);
				     hWeffTopTagging->Fill(1, efficiencyW);
				     
				     //b-veto
				     if (bveto_mu == 1) {
				       
				       hWSoftMuVeto->Fill(1, totalW);
				       hWeffSoftMuVeto->Fill(1, efficiencyW);
				       
				       hHtAfter[3]->Fill(Ht,totalW);				    
				       
				       hPtLepton1WWLevelNoHt[3]      ->Fill(pt1,       totalW);
				       hPtLepton2WWLevelNoHt[3]      ->Fill(pt2,       totalW);
				       hPtDiLeptonWWLevelNoHt[3]     ->Fill(ptll,      totalW);
				       hMinvWWLevelNoHt[3]           ->Fill(mll,       totalW);
				       hMtWWLevelNoHt[3]             ->Fill(mth,       totalW);
				       hpfMetWWLevelNoHt[3]          ->Fill(pfType1Met,totalW);
				       hpminMetWWLevelNoHt[3]        ->Fill(mpmet,     totalW);
				       hDeltaRLeptonsWWLevelNoHt[3]  ->Fill(drll,      totalW);
				       hDeltaPhiLeptonsWWLevelNoHt[3]->Fill(dphill,    totalW);
				       hDPhiPtllJetWWLevelNoHt[3]    ->Fill(dphilljet, totalW);
				       hWnJetsBveto                  ->Fill(njet,      totalW);
				       hWeffnJetsBveto               ->Fill(njet, efficiencyW);
				       hSigMuNoHt[3]                 ->Fill(std_vector_lepton_muSIP3D->at(0),totalW);
				       hSigElNoHt[3]                 ->Fill(std_vector_lepton_elSIP3D->at(0),totalW);
				       
				       //bveto Ht 
				       if(Ht < 250){
					 
					 hPtLepton1WWLevel[3]      ->Fill(pt1,       totalW);
					 hPtLepton2WWLevel[3]      ->Fill(pt2,       totalW);
					 hPtDiLeptonWWLevel[3]     ->Fill(ptll,      totalW);
					 hMinvWWLevel[3]           ->Fill(mll,       totalW);
					 hMtWWLevel[3]             ->Fill(mth,       totalW);
					 hpfMetWWLevel[3]          ->Fill(pfType1Met,totalW);
					 hpminMetWWLevel[3]        ->Fill(mpmet,     totalW);
					 hDeltaRLeptonsWWLevel[3]  ->Fill(drll,      totalW);
					 hDeltaPhiLeptonsWWLevel[3]->Fill(dphill,    totalW);
					 hDPhiPtllJetWWLevel[3]    ->Fill(dphilljet, totalW);
					 hWeffnJetsBvetoAfterHt    ->Fill(njet, efficiencyW);					
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
					 hpfMetWWLevelHtPlus[3]          ->Fill(pfType1Met,totalW);
					 hpminMetWWLevelHtPlus[3]        ->Fill(mpmet,     totalW);
					 hDeltaRLeptonsWWLevelHtPlus[3]  ->Fill(drll,      totalW);
					 hDeltaPhiLeptonsWWLevelHtPlus[3]->Fill(dphill,    totalW);
					 hDPhiPtllJetWWLevelHtPlus[3]    ->Fill(dphilljet, totalW);
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
					   hpfMetWWLevelNoHt[jetNumber]          ->Fill(pfType1Met,totalW);
					   hpminMetWWLevelNoHt[jetNumber]        ->Fill(mpmet,     totalW);
					   hDeltaRLeptonsWWLevelNoHt[jetNumber]  ->Fill(drll,      totalW);
					   hDeltaPhiLeptonsWWLevelNoHt[jetNumber]->Fill(dphill,    totalW);
					   hDPhiPtllJetWWLevelNoHt[jetNumber]    ->Fill(dphilljet, totalW);
					   hSigMuNoHt[jetNumber]                 ->Fill(std_vector_lepton_muSIP3D->at(0),totalW);
					   hSigElNoHt[jetNumber]                 ->Fill(std_vector_lepton_elSIP3D->at(0),totalW);
					   
					   //bveto Ht  
					   if(Ht < 250){
					     
					     hPtLepton1WWLevel[jetNumber]      ->Fill(pt1,       totalW);
					     hPtLepton2WWLevel[jetNumber]      ->Fill(pt2,       totalW);
					     hPtDiLeptonWWLevel[jetNumber]     ->Fill(ptll,      totalW);
					     hMinvWWLevel[jetNumber]           ->Fill(mll,       totalW);
					     hMtWWLevel[jetNumber]             ->Fill(mth,       totalW);
					     hpfMetWWLevel[jetNumber]          ->Fill(pfType1Met,totalW);
					     hpminMetWWLevel[jetNumber]        ->Fill(mpmet,     totalW);
					     hDeltaRLeptonsWWLevel[jetNumber]  ->Fill(drll,      totalW);
					     hDeltaPhiLeptonsWWLevel[jetNumber]->Fill(dphill,    totalW);
					     hDPhiPtllJetWWLevel[jetNumber]    ->Fill(dphilljet, totalW);
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
					     hpfMetWWLevelHtPlus[jetNumber]          ->Fill(pfType1Met,totalW);
					     hpminMetWWLevelHtPlus[jetNumber]        ->Fill(mpmet,     totalW);
					     hDeltaRLeptonsWWLevelHtPlus[jetNumber]  ->Fill(drll,      totalW);
					     hDeltaPhiLeptonsWWLevelHtPlus[jetNumber]->Fill(dphill,    totalW);
					     hDPhiPtllJetWWLevelHtPlus[jetNumber]    ->Fill(dphilljet, totalW);
					     hSigMuHtPlus[jetNumber]                 ->Fill(std_vector_lepton_muSIP3D->at(0),totalW);
					     hSigElHtPlus[jetNumber]                 ->Fill(std_vector_lepton_elSIP3D->at(0),totalW);
					   }
					 }  					
				       }
				     }
				   }
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
  hNjetsTwoLeptonsLevel       = ((TH1F*) FindOutput("hNjetsTwoLeptonsLevel"));
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
  hpfMetWWLevel[qq]           = ((TH1F*) FindOutput("hpfMetWWLevel%.1i",qq));
  hpminMetWWLevel[qq]         = ((TH1F*) FindOutput("hpminMetWWLevel%.1i",qq));
  hDeltaRLeptonsWWLevel[qq]   = ((TH1F*) FindOutput("hDeltaRLeptonsWWLevel%.1i",qq));
  hDeltaPhiLeptonsWWLevel[qq] = ((TH1F*) FindOutput("hDeltaPhiLeptonsWWLevel%.1i",qq));
  hDPhiPtllJetWWLevel[qq]     = ((TH1F*) FindOutput("hDPhiPtllJetWWLevel%.1i",qq));
  hSigEl[qq]                  = ((TH1F*) FindOutput("hSigEl%.1i",qq));
  hSigMu[qq]                  = ((TH1F*) FindOutput("hSigMu%.1i",qq));

  hPtLepton1WWLevelNoHt[qq]       = ((TH1F*) FindOutput("hPtLepton1WWLevelNoHt%.1i",qq));
  hPtLepton2WWLevelNoHt[qq]       = ((TH1F*) FindOutput("hPtLepton2WWLevelNoHt%.1i",qq));
  hPtDiLeptonWWLevelNoHt[qq]      = ((TH1F*) FindOutput("hPtDiLeptonWWLevelNoHt%.1i",qq));
  hMinvWWLevelNoHt[qq]            = ((TH1F*) FindOutput("hMinvWWLevelNoHt%.1i",qq));     // %.1i",qq 
  hMtWWLevelNoHt[qq]              = ((TH1F*) FindOutput("hMtWWLevelNoHt%.1i",qq));
  hpfMetWWLevelNoHt[qq]           = ((TH1F*) FindOutput("hpfMetWWLevelNoHt%.1i",qq));
  hpminMetWWLevelNoHt[qq]         = ((TH1F*) FindOutput("hpminMetWWLevelNoHt%.1i",qq));
  hDeltaRLeptonsWWLevelNoHt[qq]   = ((TH1F*) FindOutput("hDeltaRLeptonsWWLevelNoHt%.1i",qq));
  hDeltaPhiLeptonsWWLevelNoHt[qq] = ((TH1F*) FindOutput("hDeltaPhiLeptonsWWLevelNoHt%.1i",qq));
  hDPhiPtllJetWWLevelNoHt[qq]     = ((TH1F*) FindOutput("hDPhiPtllJetWWLevelNoHt%.1i",qq));
  hSigElNoHt[qq]                  = ((TH1F*) FindOutput("hSigElNoHt%.1i",qq));
  hSigMuNoHt[qq]                  = ((TH1F*) FindOutput("hSigMuNoHt%.1i",qq));

  hPtLepton1WWLevelHtPlus[qq]       = ((TH1F*) FindOutput("hPtLepton1WWLevelHtPlus%.1i",qq));
  hPtLepton2WWLevelHtPlus[qq]       = ((TH1F*) FindOutput("hPtLepton2WWLevelHtPlus%.1i",qq));
  hPtDiLeptonWWLevelHtPlus[qq]      = ((TH1F*) FindOutput("hPtDiLeptonWWLevelHtPlus%.1i",qq));
  hMinvWWLevelHtPlus[qq]            = ((TH1F*) FindOutput("hMinvWWLevelHtPlus%.1i",qq));     // HtPlus%.1i",qq 
  hMtWWLevelHtPlus[qq]              = ((TH1F*) FindOutput("hMtWWLevelHtPlus%.1i",qq));
  hpfMetWWLevelHtPlus[qq]           = ((TH1F*) FindOutput("hpfMetWWLevelHtPlus%.1i",qq));
  hpminMetWWLevelHtPlus[qq]         = ((TH1F*) FindOutput("hpminMetWWLevelHtPlus%.1i",qq));
  hDeltaRLeptonsWWLevelHtPlus[qq]   = ((TH1F*) FindOutput("hDeltaRLeptonsWWLevelHtPlus%.1i",qq));
  hDeltaPhiLeptonsWWLevelHtPlus[qq] = ((TH1F*) FindOutput("hDeltaPhiLeptonsWWLevelHtPlus%.1i",qq));
  hDPhiPtllJetWWLevelHtPlus[qq]     = ((TH1F*) FindOutput("hDPhiPtllJetWWLevelHtPlus%.1i",qq));
  hSigElHtPlus[qq]                  = ((TH1F*) FindOutput("hSigElHtPlus%.1i",qq));
  hSigMuHtPlus[qq]                  = ((TH1F*) FindOutput("hSigMuHtPlus%.1i",qq));

  hHt[qq]                     = ((TH1F*) FindOutput("hHt%.1i",qq));
  hHtAfter[qq]                = ((TH1F*) FindOutput("hHtAfter%.1i",qq));
 }

  // TwoLeptons level histograms                                                                                                                                                   
  //----------------------------------------------------------------------------                                                                                                   
  hPtLepton1TwoLeptonsLevel       = ((TH1F*) FindOutput("hPtLepton1TwoLeptonsLevel"));
  hPtLepton2TwoLeptonsLevel       = ((TH1F*) FindOutput("hPtLepton2TwoLeptonsLevel"));
  hPtDiLeptonTwoLeptonsLevel      = ((TH1F*) FindOutput("hPtDiLeptonTwoLeptonsLevel"));
  hMinvTwoLeptonsLevel            = ((TH1F*) FindOutput("hMinvTwoLeptonsLevel"));
  hMtTwoLeptonsLevel              = ((TH1F*) FindOutput("hMtTwoLeptonsLevel"));
  hpfMetTwoLeptonsLevel           = ((TH1F*) FindOutput("hpfMetTwoLeptonsLevel"));
  hpminMetTwoLeptonsLevel         = ((TH1F*) FindOutput("hpminMetTwoLeptonsLevel"));
  hDeltaRLeptonsTwoLeptonsLevel   = ((TH1F*) FindOutput("hDeltaRLeptonsTwoLeptonsLevel"));
  hDeltaPhiLeptonsTwoLeptonsLevel = ((TH1F*) FindOutput("hDeltaPhiLeptonsTwoLeptonsLevel"));
  hDPhiPtllJetTwoLeptonsLevel     = ((TH1F*) FindOutput("hDPhiPtllJetTwoLeptonsLevel"));

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


float test::Testing(int k)
{
  float my_pt = 105.;

  my_pt += std_vector_lepton_pt->at(k);

  return my_pt;
}

//------------------------------------------------------------------------------
// IsTightLepton
//------------------------------------------------------------------------------
bool test::IsTightLepton(int k)
{
  bool is_tight_lepton = false;

  // Muon tight ID
  if (fabs(std_vector_lepton_id->at(k)) == 13)
    {
      is_tight_lepton = std_vector_lepton_isTightMuon->at(k);
    }
  // Electron cut based medium ID
  else if (fabs(std_vector_lepton_id->at(k)) == 11)
    {
      float aeta = fabs(std_vector_electron_scEta->at(k));

      if (aeta <= 1.479)
	{
	  if (fabs(std_vector_electron_deltaEtaIn->at(k)) < 0.008925 &&
	      fabs(std_vector_electron_deltaPhiIn->at(k)) < 0.035973 &&
	      std_vector_electron_sigmaIetaIeta->at(k)    < 0.009996 &&
	      std_vector_electron_HoE->at(k)              < 0.050537 &&
	      fabs(std_vector_electron_d0->at(k))         < 0.012235 &&
	      fabs(std_vector_electron_dz->at(k))         < 0.042020 &&
	      fabs(std_vector_electron_ooEooP->at(k))     < 0.091942 &&
	      ElectronIsolation(k)                        < 0.107587 &&
	      !std_vector_electron_passConversion->at(k))  // Includes expectedMissingInnerHits
	    {
	      is_tight_lepton = true;
	    }
	}
      else if (aeta > 1.479 && aeta < 2.5)
	{
	  if (fabs(std_vector_electron_deltaEtaIn->at(k)) < 0.007429 &&
	      fabs(std_vector_electron_deltaPhiIn->at(k)) < 0.067879 &&
	      std_vector_electron_sigmaIetaIeta->at(k)    < 0.030135 &&
	      std_vector_electron_HoE->at(k)              < 0.086782 &&
	      fabs(std_vector_electron_d0->at(k))         < 0.036719 &&
	      fabs(std_vector_electron_dz->at(k))         < 0.138142 &&
	      fabs(std_vector_electron_ooEooP->at(k))     < 0.100683 &&
	      ElectronIsolation(k)                        < 0.113254 &&
	      !std_vector_electron_passConversion->at(k))  // Includes expectedMissingInnerHits
	    {
	      is_tight_lepton = true;
	    }
	}
    }

  return is_tight_lepton;
}


//------------------------------------------------------------------------------
// MuonIsolation
//------------------------------------------------------------------------------
float test::MuonIsolation(int k)
{
  float pt = std_vector_lepton_pt->at(k);
  float id = std_vector_lepton_id->at(k);

  float relative_isolation = -999;

  if (fabs(id) != 13) return relative_isolation;

  relative_isolation =
    std_vector_lepton_chargedHadronIso->at(k) +
    max(float(0.0),
	float(std_vector_lepton_photonIso->at(k) +
	      std_vector_lepton_neutralHadronIso->at(k) -
	      0.5*std_vector_lepton_sumPUPt->at(k)));

  relative_isolation /= pt;
  
  return relative_isolation;
}


//------------------------------------------------------------------------------
// ElectronIsolation
//------------------------------------------------------------------------------
float test::ElectronIsolation(int k)
{
  float pt = std_vector_lepton_pt->at(k);
  float id = std_vector_lepton_id->at(k);

  float relative_isolation = -999;

  if (fabs(id) != 11) return relative_isolation;

  relative_isolation =
    std_vector_lepton_chargedHadronIso->at(k) +
    max(float(0.0),
	float(std_vector_lepton_photonIso->at(k) +
	      std_vector_lepton_neutralHadronIso->at(k) -
	      jetRho*std_vector_electron_effectiveArea->at(k)));
  
  relative_isolation /= pt;
  
  return relative_isolation;
}


//------------------------------------------------------------------------------
// IsIsolatedLepton
//------------------------------------------------------------------------------
bool test::IsIsolatedLepton(int k)
{
  float id = std_vector_lepton_id->at(k);

  bool is_isolated_lepton = false;

  if      (fabs(id) == 11) is_isolated_lepton = true;//(ElectronIsolation(k) < 0.15);
  else if (fabs(id) == 13) is_isolated_lepton = (MuonIsolation(k)     < 0.12);
  
  return is_isolated_lepton;
}
