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

float dist(float phi1, float eta1, float phi2, float eta2){
  float dphi = fabs(phi1 - phi2);
  float deta = fabs(eta1 - eta2);
  if (dphi > TMath::Pi())
    dphi = 2*TMath::Pi() - dphi;
  return sqrt(dphi*dphi + deta*deta);
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
  // Create tree and branches
  //------------------------------------------------------------------------------
  
  tree = CreateTree("nt","nt");
  /*
  tree->Branch("pt1",&pt1,"pt1");
  tree->Branch("pt2",&pt2,"pt2");
  tree->Branch("pfType1Met",&pfType1Met,"pfType1Met");
  tree->Branch("trkMet",&trkMet,"trkMet");
  tree->Branch("jetpt1",&jetpt1,"jetpt1");
  tree->Branch("ptll",&ptll,"ptll");
  tree->Branch("dphill",&dphill,"dphill");
  tree->Branch("jetphi1",&jetphi1,"jetphi1");
  tree->Branch("pfType1Metphi",&pfType1Metphi,"pfType1Metphi");
  tree->Branch("mth",&mth,"mth");
  tree->Branch("mll",&mll,"mll");
  tree->Branch("mpmet",&mpmet,"mpmet");
  tree->Branch("metvar",&metvar,"metvar");
  tree->Branch("dphilljet",&dphilljet,"dphilljet");
  tree->Branch("dphillmet",&dphillmet,"dphillmet");
  tree->Branch("njet",&njet,"njet");
  tree->Branch("dphilmet1",&dphilmet1,"dphilmet1");
  tree->Branch("dphilmet2",&dphilmet2,"dphilmet2");
  tree->Branch("bveto_ip",&bveto_ip,"bveto_ip");
  tree->Branch("nbjettche",&nbjettche,"nbjettche");
  tree->Branch("dphiv",&dphiv,"dphiv");
  tree->Branch("dphijet1met",&dphijet1met,"dphijet1met");
  tree->Branch("ratioMet",&ratioMet,"ratioMet");

  tree->Branch("nvtx",&nvtx,"nvtx");
  tree->Branch("nextra",&nextra,"nextra");

  tree->Branch("std_vector_lepton_pt",&std_vector_lepton_pt,"std_vector_lepton_pt");
  tree->Branch("std_vector_jet_pt",&std_vector_jet_pt,"std_vector_jet_pt");
  */

  //------------------------------------------------------------------------------
  // Create histos
  //------------------------------------------------------------------------------
  
  h_n_PV = CreateH1F("h_n_PV","h_n_PV",50,0.,10.);

  // Counting histograms    
  //---------------------------------------------------------------------------- 
  
  hWTrigger     = CreateH1F("hWTrigger",                        "", 10, 0, 10);
  hWMetCut      = CreateH1F("hWMetCut",                         "", 10, 0, 10);
  hWLowMinv     = CreateH1F("hWLowMinv",                        "", 10, 0, 10);
  hWZVeto       = CreateH1F("hWZVeto",                          "", 10, 0, 10);
  hWpMetCut     = CreateH1F("hWpMetCut",                        "", 10, 0, 10);
  hWJetVeto     = CreateH1F("hWJetVeto",                        "", 10, 0, 10);
  hWnJets       = CreateH1F("hWnJets",                          "", 10, 0, 10);
  hWeffnJets    = CreateH1F("hWeffnJets",                       "", 10, 0, 10);
  hWnBtaggedJets     = CreateH1F("hWnBtaggedJets",              "", 10, 0, 10);
  hWeffnBtaggedJets  = CreateH1F("hWeffnBtaggedJets",           "", 10, 0, 10);
  hWnJetsBveto    = CreateH1F("hWnJetsBveto",                   "", 10, 0, 10);
  hWeffnJetsBveto = CreateH1F("hWeffnJetsBveto",                "", 10, 0, 10);
  hNjetsTwoLeptonsLevel    = CreateH1F("hNjetsTwoLeptonsLevel", "", 10, 0, 10);
  hWeffnJetsBvetoAfterHt = CreateH1F("hWeffnJetsBvetoAfterHt",  "", 10, 0, 10);
  
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
  
  hLooseIso = CreateH1F("hLooseIso",  "", 100, 0, 10);

  // WW level histograms     
  //----------------------------------------------------------------------------  
  
  for (Int_t nC=0; nC<4; nC++) {

    hPtLepton1WWLevel[nC]       = CreateH1F(Form("hPtLepton1WWLevel%.1i", nC),        "", 200, 0, 200);
    hPtLepton2WWLevel[nC]       = CreateH1F(Form("hPtLepton2WWLevel%.1i", nC),        "", 200, 0, 200);
    hPtDiLeptonWWLevel[nC]      = CreateH1F(Form("hPtDiLeptonWWLevel%.1i", nC),       "", 200, 0, 200);
    hMinvWWLevel[nC]            = CreateH1F(Form("hMinvWWLevel%.1i", nC),             "", 200, 0, 200);
    hMtWWLevel[nC]              = CreateH1F(Form("hMtWWLevel%.1i", nC),               "", 250, 0, 250);
    hpfMetWWLevel[nC]           = CreateH1F(Form("hpfMetWWLevel%.1i", nC),            "", 150, 0, 150);
    hpminMetWWLevel[nC]         = CreateH1F(Form("hpminMetWWLevel%.1i", nC),          "", 150, 0, 150);
    hDeltaRLeptonsWWLevel[nC]   = CreateH1F(Form("hDeltaRLeptonsWWLevel%.1i", nC),    "",  50, 0,   5);
    hDeltaPhiLeptonsWWLevel[nC] = CreateH1F(Form("hDeltaPhiLeptonsWWLevel%.1i", nC),  "",  32, 0, 3.2);
    hDPhiPtllJetWWLevel[nC]     = CreateH1F(Form("hDPhiPtllJetWWLevel%.1i", nC),      "",  32, 0, 3.2);
    hSigEl[nC]                  = CreateH1F(Form("hSigEl%.1i", nC),                   "", 1000,0.,500);
    hSigMu[nC]                  = CreateH1F(Form("hSigMu%.1i", nC),                   "", 1000,0.,500);

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
  hNjetsPlot1TwoLeptonsLevel      = CreateH1F("hNjetsPlot1TwoLeptonsLevel", "", 10, 0, 10);
  hNjetsPlot2TwoLeptonsLevel      = CreateH1F("hNjetsPlot2TwoLeptonsLevel", "", 10, 0, 10);
  hSigMuNoHtTwoLeptonsLevel       = CreateH1F("hSigMuNoHtTwoLeptonsLevel", "", 10, 0, 10);
  hSigElNoHtTwoLeptonsLevel       = CreateH1F("hSigElNoHtTwoLeptonsLevel", "", 10, 0, 10);
  hDxyTwoLeptonsLevel             = CreateH1F("hDxyTwoLeptonsLevel",       "", 1000, -0.1, 0.1);
  hDzTwoLeptonsLevel              = CreateH1F("hDzTwoLeptonsLevel",        "", 1000, -0.5, 0.5);
}

// The InsideLoop() function is called for each entry in the tree to be processed  
void test::InsideLoop() {

  int lep_size = std_vector_lepton_pt->size();

  //for (int i=0; i<lep_size; i++) printf("testing %f %f\n", Testing(i), std_vector_lepton_pt->at(i)+105.);

  TString TheSample = GetInputParameters()->TheNamedString("theSample");
  TString sameSign  = GetInputParameters()->TheNamedString("SameSign");
  TString _MuonID   = GetInputParameters()->TheNamedString("MuonID");
  
  //Creating the variables we need 
  //--------------------------------------------------------------------------

  Double_t efficiencyW = effW * triggW;
  Double_t totalW      = -999;
  
  efficiencyW = puW * effW * triggW ;
  //totalW      = (1 + 0.5 * (dataset >= 82 && dataset <= 84)) * baseW * efficiencyW * Luminosity;
  totalW = baseW * Luminosity;

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
  
   //building Ht
   Float_t Ht = 0.;
   Ht = std_vector_lepton_pt->at(0) + std_vector_lepton_pt->at(1) + pfType1Met;
   
   if(njet > 10) njet = 10;
   for (int i = 0; i < njet; ++i)
     if(std_vector_jet_pt->at(i) > 0)
       Ht += std_vector_jet_pt->at(i);

   //defining nextra
   nextra = 0;
   for (int i = 0; i < std_vector_lepton_pt -> size(); ++i)
     if(std_vector_lepton_pt -> at(i) > 10)
       ++nextra;
   nextra = nextra -2;
   
   //building dRjet1
   Float_t dRjet1  = 100.;
   TLorentzVector vjet1(0.,0.,0.,0.);
   TLorentzVector vlep1(0.,0.,0.,0.);
   if(std_vector_lepton_pt->at(0) > 0.){
     vlep1.SetPtEtaPhiM(std_vector_lepton_pt->at(0),std_vector_lepton_eta->at(0),std_vector_lepton_phi->at(0),0.);
     for (int i = 0; i < std_vector_jet_pt -> size(); ++i)
       if(std_vector_jet_pt -> at(i) > 0.){
	 vjet1.SetPtEtaPhiM(std_vector_jet_pt->at(i),std_vector_jet_eta->at(i),std_vector_jet_phi->at(i),0.);
	 if( vlep1.DeltaR(vjet1) < dRjet1){
	   dRjet1 = vlep1.DeltaR(vjet1);
	 }
       }
   }

   //building dRjet2
   Float_t dRjet2  = 100.;
   TLorentzVector vjet2(0.,0.,0.,0.);
   TLorentzVector vlep2(0.,0.,0.,0.);
   if(std_vector_lepton_pt->at(1) > 0.){
     vlep2.SetPtEtaPhiM(std_vector_lepton_pt->at(1),std_vector_lepton_phi->at(1),std_vector_lepton_eta->at(1),0.);
     for (int i = 0; i < njet; ++i)
       if(std_vector_jet_pt -> at(i) > 0.){
	 vjet2.SetPtEtaPhiM(std_vector_jet_pt->at(i),std_vector_jet_eta->at(i),std_vector_jet_phi->at(i),0.);
	 if( vlep2.DeltaR(vjet2) < dRjet2){
	   dRjet2 = vlep2.DeltaR(vjet2);
	 }
       }
   }
   
   //Building dphijet1met
   Float_t dphijet1met = 0.;
   if (jetphi1 > 0 && pfType1Metphi > 0){
     dphijet1met = fabs(jetphi1 - pfType1Metphi);
     if (dphijet1met > TMath::Pi()) dphijet1met = 2*TMath::Pi() - dphijet1met;
   }

   //Building RatioMet   
   Float_t ratioMet = 0.;
   if (pfType1Met > 0 && trkMet > 0)
     ratioMet = pfType1Met / sqrt(pfType1Met + trkMet);

// The selection begins here
   //--------------------------------------------------------------------------
   if (std_vector_lepton_pt->at(0) > 20)
     if (std_vector_lepton_pt->at(1) > 10) 
       if ((sameSign == "SS" && ch1*ch2 > 0) || (sameSign == "OS" && ch1*ch2 < 0))
	 if ( (SelectedChannel == -1)                                   || 
	      (channel == SelectedChannel)                          || 
	      (SelectedChannel == 4 && (channel == 2 || channel == 3) ) || 
	      (SelectedChannel == 5 && (channel == 0 || channel == 1) ) 
	      ){
	   /*
	   if (IsTightLepton(0,_MuonID) && !IsTightLepton(1,_MuonID))
	     hLooseIso -> Fill(ElectronIsolation(1), totalW);
	   if (IsTightLepton(1,_MuonID) && !IsTightLepton(0,_MuonID))
	     hLooseIso -> Fill(ElectronIsolation(0), totalW);
	   */
	   
	   hDxyTwoLeptonsLevel ->Fill(std_vector_lepton_BestTrackdxy->at(0),totalW);	       
	   hDzTwoLeptonsLevel  ->Fill(std_vector_lepton_BestTrackdz ->at(0),totalW);	       
	   
	   if (IsIsolatedLepton(0))
	     if (IsIsolatedLepton(1))
	       if (IsTightLepton(0,_MuonID))
		 if (IsTightLepton(1,_MuonID)){
		   
		   
		   //		   tree->Fill();
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
		   hNjetsPlot1TwoLeptonsLevel     ->Fill(dRjet1,    totalW);
		   hNjetsPlot2TwoLeptonsLevel     ->Fill(dRjet2,    totalW);
		   hSigMuNoHtTwoLeptonsLevel      ->Fill(std_vector_lepton_muSIP3D->at(0),totalW);
		   hSigElNoHtTwoLeptonsLevel      ->Fill(std_vector_lepton_elSIP3D->at(0),totalW);
		   
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
				   hHt[2]->Fill(Ht,totalW);			       
				   
				   //b-veto (now not operative)
				   /*if (bveto_mu == 1)*/ {
				     
				     hWSoftMuVeto->Fill(1, totalW);
				     hWeffSoftMuVeto->Fill(1,efficiencyW);
				     
				     hHtAfter[3]->Fill(Ht,totalW);				    
				     
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
				     
				     for (Int_t jetNumber = 0; jetNumber < 3 ; ++jetNumber){
				       if (jetbin >= 3) jetbin = 2;
				       if(jetNumber == jetbin){
					 
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

  hLooseIso = ((TH1F*) FindOutput("hLooseIso"));
  

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
  hNjetsPlot1TwoLeptonsLevel      = ((TH1F*) FindOutput("hNjetsPlot1TwoLeptonsLevel"));
  hNjetsPlot2TwoLeptonsLevel      = ((TH1F*) FindOutput("hNjetsPlot2TwoLeptonsLevel"));
  hSigMuNoHtTwoLeptonsLevel       = ((TH1F*) FindOutput("hSigMuNoHtTwoLeptonsLevel"));
  hSigElNoHtTwoLeptonsLevel       = ((TH1F*) FindOutput("hSigElNoHtTwoLeptonsLevel"));
  hDxyTwoLeptonsLevel             = ((TH1F*) FindOutput("hDxyTwoLeptonsLevel"));
  hDzTwoLeptonsLevel              = ((TH1F*) FindOutput("hDzTwoLeptonsLevel"));
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
//
// https://twiki.cern.ch/twiki/bin/view/CMS/CutBasedElectronIdentificationRun2
// egmGsfElectronIDs:cutBasedElectronID-PHYS14-PU20bx25-V2-standalone-medium
//------------------------------------------------------------------------------
bool test::IsTightLepton(int k, TString _MuonID_)
{
  bool is_tight_lepton = false;

  // Muon ID
  if (fabs(std_vector_lepton_id->at(k)) == 13){
    
    if (_MuonID_ == "MediumID"){
      if (std_vector_lepton_isMediumMuon->at(k) == 1)
	is_tight_lepton = true;
    }
    
    else if(_MuonID_ == "MediumIDTighterIP"){
      if (std_vector_lepton_isMediumMuon->at(k) == 1)
	if (fabs(std_vector_lepton_BestTrackdxy -> at(k)) < 0.02)
	  if (fabs(std_vector_lepton_BestTrackdz -> at(k)) < 0.1)
	    is_tight_lepton = true;
    }
    
    else if (_MuonID_ == "TightID" ){
      if (std_vector_lepton_isTightMuon->at(k) == 1)
	is_tight_lepton = true;
    }
    
    else if(_MuonID_ == "TightIDTighterIP"){
      if (std_vector_lepton_isTightMuon->at(k) == 1)
	if (fabs(std_vector_lepton_BestTrackdxy -> at(k)) < 0.02)
	  if (fabs(std_vector_lepton_BestTrackdz -> at(k)) < 0.1)
	    is_tight_lepton = true;
    }
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

  if      (fabs(id) == 11) is_isolated_lepton = true;  //(ElectronIsolation(k) < 0.15);
  else if (fabs(id) == 13) is_isolated_lepton = (MuonIsolation(k)     < 0.12);
  
  return is_isolated_lepton;
}
