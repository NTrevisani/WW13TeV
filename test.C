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
    //    hNinZevents     [nC] = CreateH1F(Form("hNinZevents%.1i",      nC+1 ), "",   3, 0,   3);
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
  hNJets30TwoLeptonsLevel         = CreateH1F("hNJetsPF30TwoLeptonsLevel",       "",  10, 0,  10);
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

  hEff                          = CreateH1F("hEff",                          "", 50, 0 ,  10);
  hEffHt                        = CreateH1F("hEffHt",                        "", 50, 0 ,  10);

  hEff0bin                      = CreateH1F("hEff0bin",                          "", 50, 0 ,  10);
  hEffHt0bin                    = CreateH1F("hEffHt0bin",                        "", 50, 0 ,  10);

  // Dilepton level differential histograms 
  //----------------------------------------------------------------------------   
  
  //differential in dilepton (eta, pt, ptll, mll, dphill), jets (jetpt1, njet), transverse mass? 2-jet bin?   
  
  Double_t pt1bins[8] = {20,40,60,80,100,125,150,200};
  Double_t eta1bins[9] = {-2.4,-1.8,-1.2,-0.6,0,0.6,1.2,1.8,2.4};
  Double_t mllbins[9] = {20,40,60,80,100,125,150,175,200};
  Double_t dphibins[13] = {0,0.25,0.5,0.75,1,1.25,1.5,1.75,2,2.25,2.5,2.75,3};
  Double_t ptllbins[8] = {30,40,50,60,70,85,120,150};
  
  hPtLepton1DilepLevel_Diff   = CreateH1F("hPtLepton1DilepLevel_Diff",        "", 7, pt1bins);
  hEtaLepton1DilepLevel_Diff  = CreateH1F("hEtaLepton1DilepLevel_Diff",       "", 8, eta1bins);
  hDileptonDilepLevel_Diff    = CreateH1F("hDileptonDilepLevel_Diff",         "", 7, ptllbins);
  hMinvDilepLevel_Diff        = CreateH1F("hMinvDilepLevel_Diff",             "", 8, mllbins);
  hDeltaPhiDilepLevel_Diff    = CreateH1F("hDeltaPhiDilepLevel_Diff",         "", 12, dphibins);
  
  // WW level differential histograms    
  //----------------------------------------------------------------------------
  
  hPtLepton1WWLevel_Diff   = CreateH1F("hPtLepton1WWLevel_Diff",        "", 7, pt1bins);
  hEtaLepton1WWLevel_Diff  = CreateH1F("hEtaLepton1WWLevel_Diff",       "", 8, eta1bins);
  hDileptonWWLevel_Diff    = CreateH1F("hDileptonWWLevel_Diff",         "", 7, ptllbins);
  hMinvWWLevel_Diff        = CreateH1F("hMinvWWLevel_Diff",             "", 8, mllbins);
  hDeltaPhiWWLevel_Diff    = CreateH1F("hDeltaPhiWWLevel_Diff",         "", 12, dphibins);

 
  // Data-driven methods: Z+jets 
  //---------------------------------------------------------------------------- 
  
  /*   
  hNinZevents     [numberDYMVACuts];
  hNoutZevents    [numberDYMVACuts];
  hNinLooseZevents[numberDYMVACuts];
  hMassInZevents  [numberDYMVACuts];
  hMassOutZevents [numberDYMVACuts];
  */

  for (Int_t nC=0; nC<numberDYMVACuts; nC++) {
    hNinZevents     [nC] = CreateH1F(Form("hNinZevents%.1i",      nC+1 ), "",   10, 0,   10);
    hNoutZevents    [nC] = CreateH1F(Form("hNoutZevents%.1i",     nC+1 ), "",   10, 0,   10);
    hNinLooseZevents[nC] = CreateH1F(Form("hNinLooseZevents%.1i", nC+1 ), "",   10, 0,   10);
    hMassInZevents  [nC] = CreateH1F(Form("hMassInZevents%.1i",   nC+1 ), "", 200, 0, 200);
    hMassOutZevents [nC] = CreateH1F(Form("hMassOutZevents%.1i",  nC+1 ), "", 200, 0, 200);
  }
  
  // Data-driven methods: Top                                                                                                                                                      
  //----------------------------------------------------------------------------                                                                                                   
  hTopTaggedEvents            = CreateH1F("hTopTaggedEvents",            "", 10, 0, 10);
  hNTopControlRegion          = CreateH1F("hNTopControlRegion",          "", 10, 0, 10);
  hNTopTaggedTopControlRegion = CreateH1F("hNTopTaggedTopControlRegion", "", 10, 0, 10);
  
  hbTagDisTopTaggedEvents            = CreateH1F("hbTagDisTopTaggedEvents",            "", 300, -10, 20);
  hbTagDisNTopControlRegion          = CreateH1F("hbTagDisNTopControlRegion",          "", 300, -10, 20);
  hbTagDisNTopTaggedTopControlRegion = CreateH1F("hbTagDisNTopTaggedTopControlRegion", "", 300, -10, 20);
  
}


// The InsideLoop() function is called for each entry in the tree to be processed  
void test::InsideLoop() {


  TString TheSample = GetInputParameters()->TheNamedString("theSample");
  
  Double_t efficiencyW = effW * triggW;
  Double_t totalW      = -999;
  
  efficiencyW = puW * effW * triggW ;
  totalW      = (1 + 0.5 * (dataset >= 82 && dataset <= 84)) * baseW * efficiencyW * Luminosity;
  
//   h_n_PV -> Fill(1,efficiencyW);  
  
//   //cout <<mybaseW << endl;
  
//   //if (TheSample.Contains("DYtautau") && dataset != 86) totalW *= (19.4/12.1);
  
  
//   // Help variables
//   //--------------------------------------------------------------------------
//   //   Int_t dphiv = ((njet <= 1 && dphiveto) || (njet > 1 && dphilljetjet < 165.*TMath::DegToRad()));
  
   Int_t dphiv = (njet <= 1 || (njet > 1 && dphilljetjet < 165.*TMath::DegToRad()));
  
   Float_t metvar = (njet <= 1) ? mpmet : pfmet;
  
   Float_t jetbin = njet;
    
   Float_t dyMVA = ( (channel == 2 || channel == 3) || ( (njet!=0 || dymva1>0.88) && (njet!=1 || dymva1>0.84) && ( njet==0 || njet==1 || (pfmet > 45.0)) ) );
  
  Float_t Ht = 0;
  
  if (njet == 0)
    Ht = pt1 + pt2 + pfmet;
  
  
  if (njet == 1)
    Ht = pt1 + pt2 + pfmet + jetpt1;
  
  
  if (njet == 2)
    Ht = pt1 + pt2 + pfmet + jetpt1 + jetpt2;
  
  
  if (njet == 3)
    Ht = pt1 + pt2 + pfmet + jetpt1 + jetpt2 + jetpt3;
  
  
  if (njet >= 4)
    Ht = pt1 + pt2 + pfmet + jetpt1 + jetpt2 + jetpt3 + jetpt4;

   //   Float_t Ht = pt1 + pt2 + pfmet;
   //   for (int i = 0; i < njet; ++i)
   //    Ht += std_vector_jet_pt.at(i);
  
  
  Int_t vbfsel = ((jetpt3 <= 30 || !(jetpt3 > 30 && ((jeteta1-jeteta3 > 0 && jeteta2-jeteta3 < 0) || (jeteta2-jeteta3 > 0 && jeteta1-jeteta3 < 0)))));

    
  // The selection begins here
  //--------------------------------------------------------------------------
  
  if (TheSample != "DY" || mctruth != 2 || (flavorChannel != "EMu" || flavorChannel != "MuE"))
    if (run != 201191)
      if (trigger == 1)
	if (pt2 >= 20) //increase the pt of the leptons to further reduce Wjets 
	  if (pt1 >= 20) //increase the pt of the leptons to further reduce Wjets 
	    if (ch1*ch2 > 0)
	      if ( (SelectedChannel == -1)                                   || 
		   (channel == SelectedChannel)                              || 
		   (SelectedChannel == 4 && (channel == 2 || channel == 3) ) || 
		   (SelectedChannel == 5 && (channel == 0 || channel == 1) ) 
		   )//chooses the right decay channel 
		{
		  
		  // if (dataset == 86 && (flavorChannel == "MuMu" || flavorChannel == "EE")) continue; 
		  // After talking Andrea's just removing this condition!! 
		  /*
		    if ( (SelectedChannel == -1)                                   || 
		    (channel == SelectedChannel)                              || 
		    (SelectedChannel == 4 && (channel == 2 || channel == 3) ) || 
		    (SelectedChannel == 5 && (channel == 0 || channel == 1) ) 
		    ) 
		  */
		  
		  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		  //
		  // Data-driven methods
		  //
		  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		  // if (pfmet > 20 && mpmet > 20 && mll > 12 && ptll > 45 && nextra == 0 && (dphiv || !sameflav)) { OLD CUTS
		  
		  // if( pfmet > 20 && mpmet > 20 && mll > 12&& nextra == 0 && (dphiv || !sameflav) && dyMVA && ptll>30 && (!sameflav || ptll>45) ) { 
		  
		  
		  if( pfmet > 20 && mpmet > 20 && mll > 12 && nextra == 0 && (dphiv || channel == 2 || channel == 3)  && ptll>30 && (channel == 2 || channel == 3 || ptll>45) ) { 
		      
		    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		    //
		    // Z+jets
		    //
		    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		    
		    //if (dphiv && jetbin == jetChannel && bveto_mu && (bveto_ip && (nbjettche == 0 || njet > 3))) { OLD CUTS
		    
		    if( jetbin == jetChannel && bveto_mu && bveto_ip==1 && nbjettche==0 && ptll > 45) { 
		      
		      // Loop over the metvar bins -------------> why only jetbin == 0 and jetbin == 1 ?? ////Z+jets doesn't have more than 1 jet??/////
		      //--------------------------------------------------------------------
		      for (size_t mc=0; mc<numberDYMVACuts; mc ++) {
			
			if ( jetbin == 0 ) { 
			  
			  if ( dymva1 > DYMVACut_0j[mc] && fabs(mll - ZMASS) < 7.5) {
			    hNinLooseZevents[mc]->Fill(1,totalW);
			  }
			  
			  if (dymva1 > DYMVACut_0j[mc] && dymva1 < DYMVACut_0j[mc+1]) {   
			    if (fabs(mll - ZMASS) < 7.5) {
			      hNinZevents[mc]   ->Fill(  1, totalW);
			      hMassInZevents[mc]->Fill(mll, totalW);
			    }
			    else if (fabs(mll - ZMASS) > 15) {  
			      hNoutZevents[mc]   ->Fill(  1, totalW);
			      hMassOutZevents[mc]->Fill(mll, totalW);
			    }
			  }
			}
			
			if ( jetbin == 1 ) { 
			  
			  if ( dymva1 > DYMVACut_1j[mc] && fabs(mll - ZMASS) < 7.5) {
			    hNinLooseZevents[mc]->Fill(1,totalW);
			  }
			  
			  if (dymva1 > DYMVACut_1j[mc] && dymva1 < DYMVACut_1j[mc+1]) {   
			    if (fabs(mll - ZMASS) < 7.5) {
			      hNinZevents[mc]   ->Fill(  1, totalW);
			      hMassInZevents[mc]->Fill(mll, totalW);
			    }
			    else if (fabs(mll - ZMASS) > 15) {  
			      hNoutZevents[mc]   ->Fill(  1, totalW);
			      hMassOutZevents[mc]->Fill(mll, totalW);
			    }
			  }
			}
			
			
		      }
		    }
		    
		    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		    //
		    // Top
		    //
		    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		    
		    //   if (zveto  && metvar > (20 + 25*sameflav)) { OLD CUTS
		    
		    if ( fabs(mll - ZMASS) > 15 && dyMVA ) {
		      
		      // btag_eff denominator
		      if ((jetChannel == 0 && njet == 1 && nbjet == 1) || // jetBProbabilityBJetTags discriminator
			  (jetChannel == 1 && njet == 2 && jettche2 > 2.1) ||
			  (jetChannel == 2)) {
			
			hNTopControlRegion->Fill(1, totalW);
			hbTagDisNTopControlRegion->Fill(jettche2, totalW);
			
			// btag_eff numerator
			if ((jetChannel == 0 && !bveto_nj30) || // bveto_nj30 == tche>2.1 && jetPt[10,30] && softMuonVeto for jetPt<30 
			    (jetChannel == 1 && jettche1 > 2.1) ||
			    (jetChannel == 2)) {
			  
			  hNTopTaggedTopControlRegion->Fill(1, totalW);
			  hbTagDisNTopTaggedTopControlRegion->Fill(jettche2, totalW);
			}
		      }
		    }
		    
		    // Top-tagged events for ttbar estimation
		    //----------------------------------------------------------------------
		    if (fabs(mll - ZMASS) > 15 && dyMVA) {
		      
		      if ((jetChannel == 0 && njet == 0 && !bveto) ||
			  (jetChannel == 1 && njet == 1 && bveto && jettche1 > 2.1) || // bveto == bveto_nj30 == tche>2.1 && jetPt[10,30] && softMuonVeto (no pt requirement)
			  (jetChannel == 2)) {
			
			hTopTaggedEvents->Fill(1, totalW);
			hbTagDisTopTaggedEvents->Fill(jettche2, totalW);
		      }
		    }
		  }
		  
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
		  hNJets30TwoLeptonsLevel        ->Fill(njet,      totalW);
		  hpfMetTwoLeptonsLevel          ->Fill(pfmet,     totalW);
		  hppfMetTwoLeptonsLevel         ->Fill(ppfmet,    totalW);
		  hchMetTwoLeptonsLevel          ->Fill(chmet,     totalW);
		  hpchMetTwoLeptonsLevel         ->Fill(pchmet,    totalW);
		  hpminMetTwoLeptonsLevel        ->Fill(mpmet,     totalW);
		  hDeltaRLeptonsTwoLeptonsLevel  ->Fill(drll,      totalW);
		  hDeltaPhiLeptonsTwoLeptonsLevel->Fill(dphill,    totalW);
		  hDPhiPtllJetTwoLeptonsLevel    ->Fill(dphilljet, totalW);
		  
		  //if ( (pass2012ICHEP1 && !pass2012ICHEP2) || (!pass2012ICHEP1 && pass2012ICHEP2))    
		  h_TwoLeptons_TightFailEvents ->Fill(1, totalW); 
		  
		  //if ( pass2012ICHEP1 && pass2012ICHEP2) 
		  //h_TwoLeptons_TightTightEvents ->Fill(1, totalW); 
		  
		  //if ( pass2012ICHEP1 || pass2012ICHEP2) 
		  //h_TwoLeptons_TightLooseEvents ->Fill(1, totalW); 
		  
		  if (nextra == 0) {
		    
		    hWExtraLepton->Fill(1, totalW);
		    hWeffExtraLepton->Fill(1, efficiencyW);
		    
		    if (pfmet > 20 ) { // removed for differential xsec
		      
		      hWMetCut->Fill(1, totalW);
		      hWeffMetCut->Fill(1, efficiencyW);
		      
		      if (mll > 12) {
			
			hWLowMinv->Fill(1, totalW);
			hWeffLowMinv->Fill(1, efficiencyW);
			
			/*if (fabs(mll - ZMASS) > 15 || channel == 2 || channel == 3)*/ {
			  
			  hWZVeto->Fill(1, totalW); 
			  hWeffZVeto->Fill(1, efficiencyW); 
			  
			  /*if (mpmet > 20)  && dyMVA )*/ { 
			    
			    hWpMetCut->Fill(1, totalW);
			    hWeffpMetCut->Fill(1, efficiencyW);		  
			    
			    if (dphiv || channel == 2 || channel == 3) {
			      
			      hWDeltaPhiJet->Fill(1, totalW);
			      hWeffDeltaPhiJet->Fill(1, efficiencyW);
			      
			      /*if ( ptll>30 && (channel == 2 || channel == 3 || ptll>45) ) */{
				
				hWPtll->Fill(1, totalW);
				hWeffPtll->Fill(1, efficiencyW);
				
				/*if (jetbin == jetChannel)*/ {
				  
				  hWJetVeto->Fill(1, totalW);
				  hWeffJetVeto->Fill(1, efficiencyW);
				  
				  hWnJets->Fill(njet, totalW);
				  hWeffnJets->Fill(njet, efficiencyW);
				  
				  hWnBtaggedJets->Fill(nbjet, totalW);
				  hWeffnBtaggedJets->Fill(nbjet, efficiencyW);
				  
				  for (Int_t jetNumber = 0; jetNumber < 3 ; ++jetNumber){
				    if (jetbin >= 3) jetbin = 2;
				    if(jetNumber == jetbin){
				      hHt[jetNumber]->Fill(Ht,totalW);				    
				    }
				  }
				  
				  hHt[3]->Fill(Ht,totalW);				    
				  
				  //b-veto
				  if (bveto_ip==1 && nbjettche==0) {
				    
				    hWTopTagging->Fill(1, totalW);
				    hWeffTopTagging->Fill(1, efficiencyW);
				      
				    //b-veto
				    /*if (bveto_mu) */{ //--> third soft lepton with pt> 3GeV
					
				      hWSoftMuVeto->Fill(1, totalW);
				      hWeffSoftMuVeto->Fill(1, efficiencyW);
				      
				      hHtAfter[3]->Fill(Ht,totalW);				    
				      
				      hPtLepton1WWLevelNoHt[3]      ->Fill(pt1,       totalW);
				      hPtLepton2WWLevelNoHt[3]      ->Fill(pt2,       totalW);
				      hPtDiLeptonWWLevelNoHt[3]     ->Fill(ptll,      totalW);
				      hMinvWWLevelNoHt[3]           ->Fill(mll,       totalW);
				      hMtWWLevelNoHt[3]             ->Fill(mth,       totalW);
				      hNJets30WWLevelNoHt[3]        ->Fill(jetbin,    totalW);
				      hpfMetWWLevelNoHt[3]          ->Fill(pfmet,     totalW);
				      hppfMetWWLevelNoHt[3]         ->Fill(ppfmet,    totalW);
				      hchMetWWLevelNoHt[3]          ->Fill(chmet,     totalW);
				      hpchMetWWLevelNoHt[3]         ->Fill(pchmet,    totalW);
				      hpminMetWWLevelNoHt[3]        ->Fill(mpmet,     totalW);
				      hDeltaRLeptonsWWLevelNoHt[3]  ->Fill(drll,      totalW);
				      hDeltaPhiLeptonsWWLevelNoHt[3]->Fill(dphill,    totalW);
				      hDPhiPtllJetWWLevelNoHt[3]    ->Fill(dphilljet, totalW);
				      hWnJetsBveto                  ->Fill(njet,      totalW);
				      hWeffnJetsBveto               ->Fill(njet, efficiencyW);
				      hEff                          ->Fill(1,    efficiencyW);
				      
				      //bveto Ht 
				      if(Ht < 250){
					
					hPtLepton1WWLevel[3]      ->Fill(pt1,       totalW);
					hPtLepton2WWLevel[3]      ->Fill(pt2,       totalW);
					hPtDiLeptonWWLevel[3]     ->Fill(ptll,      totalW);
					hMinvWWLevel[3]           ->Fill(mll,       totalW);
					hMtWWLevel[3]             ->Fill(mth,       totalW);
					hNJets30WWLevel[3]        ->Fill(jetbin,    totalW);
					hpfMetWWLevel[3]          ->Fill(pfmet,     totalW);
					hppfMetWWLevel[3]         ->Fill(ppfmet,    totalW);
					hchMetWWLevel[3]          ->Fill(chmet,     totalW);
					hpchMetWWLevel[3]         ->Fill(pchmet,    totalW);
					hpminMetWWLevel[3]        ->Fill(mpmet,     totalW);
					hDeltaRLeptonsWWLevel[3]  ->Fill(drll,      totalW);
					hDeltaPhiLeptonsWWLevel[3]->Fill(dphill,    totalW);
					hDPhiPtllJetWWLevel[3]    ->Fill(dphilljet, totalW);
					hWnJetsBvetoAfterHt       ->Fill(njet,      totalW);
					hWeffnJetsBvetoAfterHt    ->Fill(njet, efficiencyW);					
					hEffHt                    ->Fill(1,    efficiencyW);					
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
					  hpfMetWWLevelNoHt[jetNumber]          ->Fill(pfmet,     totalW);
					  hppfMetWWLevelNoHt[jetNumber]         ->Fill(ppfmet,    totalW);
					  hchMetWWLevelNoHt[jetNumber]          ->Fill(chmet,     totalW);
					  hpchMetWWLevelNoHt[jetNumber]         ->Fill(pchmet,    totalW);
					  hpminMetWWLevelNoHt[jetNumber]        ->Fill(mpmet,     totalW);
					  hDeltaRLeptonsWWLevelNoHt[jetNumber]  ->Fill(drll,      totalW);
					  hDeltaPhiLeptonsWWLevelNoHt[jetNumber]->Fill(dphill,    totalW);
					  hDPhiPtllJetWWLevelNoHt[jetNumber]    ->Fill(dphilljet, totalW);
					  if(jetNumber == 0)
					    hEff0bin                            ->Fill(1,    efficiencyW);
					  
					  //bveto Ht  
					  if(Ht < 250){
					    
					    hPtLepton1WWLevel[jetNumber]      ->Fill(pt1,       totalW);
					    hPtLepton2WWLevel[jetNumber]      ->Fill(pt2,       totalW);
					    hPtDiLeptonWWLevel[jetNumber]     ->Fill(ptll,      totalW);
					    hMinvWWLevel[jetNumber]           ->Fill(mll,       totalW);
					    hMtWWLevel[jetNumber]             ->Fill(mth,       totalW);
					    hNJets30WWLevel[jetNumber]        ->Fill(jetNumber, totalW);
					    hpfMetWWLevel[jetNumber]          ->Fill(pfmet,     totalW);
					    hppfMetWWLevel[jetNumber]         ->Fill(ppfmet,    totalW);
					    hchMetWWLevel[jetNumber]          ->Fill(chmet,     totalW);
					    hpchMetWWLevel[jetNumber]         ->Fill(pchmet,    totalW);
					    hpminMetWWLevel[jetNumber]        ->Fill(mpmet,     totalW);
					    hDeltaRLeptonsWWLevel[jetNumber]  ->Fill(drll,      totalW);
					    hDeltaPhiLeptonsWWLevel[jetNumber]->Fill(dphill,    totalW);
					    hDPhiPtllJetWWLevel[jetNumber]    ->Fill(dphilljet, totalW);
					    if(jetNumber == 0)
					      hEffHt0bin                      ->Fill(1,    efficiencyW);
					  }
					}  
					
				      }
				      hPtLepton1WWLevel_Diff  ->Fill(pt1,       totalW);
				      hEtaLepton1WWLevel_Diff ->Fill(eta1,      totalW);
				      hDileptonWWLevel_Diff   ->Fill(ptll,      totalW);
				      hMinvWWLevel_Diff       ->Fill(mll,       totalW);
				      hDeltaPhiWWLevel_Diff   ->Fill(dphill,    totalW);
				      
				      //if ( (pass2012ICHEP1 && !pass2012ICHEP2) || (!pass2012ICHEP1 && pass2012ICHEP2))    
				      h_WWLevel_TightFailEvents ->Fill(1, totalW); 
				      
				      //if ( (pass2012ICHEP1 && pass2012ICHEP2)) 
				      //h_WWLevel_TightTightEvents ->Fill(1, totalW); 
				      
				      //if ( pass2012ICHEP1 || pass2012ICHEP2) 
				      //h_WWLevel_TightLooseEvents ->Fill(1, totalW); 
				      
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
		  
		  
		}
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
  hNJets30TwoLeptonsLevel         = ((TH1F*) FindOutput("hNJetsPF30TwoLeptonsLevel"));
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

  hEff                          = ((TH1F*) FindOutput("hEff"));
  hEffHt                        = ((TH1F*) FindOutput("hEffHt"));

  hEff0bin                          = ((TH1F*) FindOutput("hEff0bin"));
  hEffHt0bin                        = ((TH1F*) FindOutput("hEffHt0bin"));


  // Dilepton level differential histograms   
  //----------------------------------------------------------------------------   

  //differential in dilepton (eta, pt, ptll, mll, dphill), jets (jetpt1, njet), transverse mass? 2-jet bin?                                                                        

  Double_t pt1bins[8] = {20,40,60,80,100,125,150,200};
  Double_t eta1bins[9] = {-2.4,-1.8,-1.2,-0.6,0,0.6,1.2,1.8,2.4};
  Double_t mllbins[9] = {20,40,60,80,100,125,150,175,200};
  Double_t dphibins[13] = {0,0.25,0.5,0.75,1,1.25,1.5,1.75,2,2.25,2.5,2.75,3};
  Double_t ptllbins[8] = {30,40,50,60,70,85,120,150};

  hPtLepton1DilepLevel_Diff  = ((TH1F*) FindOutput("hPtLepton1DilepLevel_Diff"));
  hEtaLepton1DilepLevel_Diff  = ((TH1F*) FindOutput("hEtaLepton1DilepLevel_Diff"));
  hDileptonDilepLevel_Diff  = ((TH1F*) FindOutput("hDileptonDilepLevel_Diff"));
  hMinvDilepLevel_Diff  = ((TH1F*) FindOutput("hMinvDilepLevel_Diff"));
  hDeltaPhiDilepLevel_Diff  = ((TH1F*) FindOutput("hDeltaPhiDilepLevel_Diff"));


  // WW level differential histograms   
  //---------------------------------------------------------------------------

  hPtLepton1WWLevel_Diff  = ((TH1F*) FindOutput("hPtLepton1WWLevel_Diff"));
  hEtaLepton1WWLevel_Diff  = ((TH1F*) FindOutput("hEtaLepton1WWLevel_Diff"));
  hDileptonWWLevel_Diff  = ((TH1F*) FindOutput("hDileptonWWLevel_Diff"));
  hMinvWWLevel_Diff  = ((TH1F*) FindOutput("hMinvWWLevel_Diff"));
  hDeltaPhiWWLevel_Diff  = ((TH1F*) FindOutput("hDeltaPhiWWLevel_Diff"));




  // Data-driven methods: Z+jets       
  //----------------------------------------------------------------------------  
  /* 
  hNinZevents     [numberDYMVACuts];
  hNoutZevents    [numberDYMVACuts];
  hNinLooseZevents[numberDYMVACuts];
  hMassInZevents  [numberDYMVACuts];
  hMassOutZevents [numberDYMVACuts];
  */

  for (Int_t nC=0; nC<numberDYMVACuts; nC++) {
    hNinZevents     [nC] = ((TH1F*) FindOutput(Form("hNinZevents%.1i",      nC+1 )));
    hNoutZevents    [nC] = ((TH1F*) FindOutput(Form("hNoutZevents%.1i",     nC+1 )));
    hNinLooseZevents[nC] = ((TH1F*) FindOutput(Form("hNinLooseZevents%.1i", nC+1 )));
    hMassInZevents  [nC] = ((TH1F*) FindOutput(Form("hMassInZevents%.1i",   nC+1 )));
    hMassOutZevents [nC] = ((TH1F*) FindOutput(Form("hMassOutZevents%.1i",  nC+1 )));
  }

    // Data-driven methods: Top                                                                                                                                                      
    //----------------------------------------------------------------------------                                                                                                   
    hTopTaggedEvents            = ((TH1F*) FindOutput("hTopTaggedEvents"));
    hNTopControlRegion          = ((TH1F*) FindOutput("hNTopControlRegion"));
    hNTopTaggedTopControlRegion = ((TH1F*) FindOutput("hNTopTaggedTopControlRegion"));

    hbTagDisTopTaggedEvents            = ((TH1F*) FindOutput("hbTagDisTopTaggedEvents"));
    hbTagDisNTopControlRegion          = ((TH1F*) FindOutput("hbTagDisNTopControlRegion"));
    hbTagDisNTopTaggedTopControlRegion = ((TH1F*) FindOutput("hbTagDisNTopTaggedTopControlRegion"));


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
