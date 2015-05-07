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

  h_n_PV = CreateH1F("h_n_PV","h_n_PV",10,-5.,5.);

  // Counting histograms    
  //---------------------------------------------------------------------------- 
                                                               
  hWTrigger     = CreateH1F("hWTrigger",     "", 3, 0, 3);
  hWMetCut      = CreateH1F("hWMetCut",      "", 3, 0, 3);
  hWLowMinv     = CreateH1F("hWLowMinv",     "", 3, 0, 3);
  hWZVeto       = CreateH1F("hWZVeto",       "", 3, 0, 3);
  hWpMetCut     = CreateH1F("hWpMetCut",     "", 3, 0, 3);
  hWJetVeto     = CreateH1F("hWJetVeto",     "", 3, 0, 3);
  hWDeltaPhiJet = CreateH1F("hWDeltaPhiJet", "", 3, 0, 3);
  hWSoftMuVeto  = CreateH1F("hWSoftMuVeto",  "", 3, 0, 3);
  hWExtraLepton = CreateH1F("hWExtraLepton", "", 3, 0, 3);
  hWPtll        = CreateH1F("hWPtll",        "", 3, 0, 3);
  hWTopTagging  = CreateH1F("hWTopTagging",  "", 3, 0, 3);

  hWeffTrigger     = CreateH1F("hWeffTrigger",     "", 3, 0, 3);
  hWeffMetCut      = CreateH1F("hWeffMetCut",      "", 3, 0, 3);
  hWeffLowMinv     = CreateH1F("hWeffLowMinv",     "", 3, 0, 3);
  hWeffZVeto       = CreateH1F("hWeffZVeto",       "", 3, 0, 3);
  hWeffpMetCut     = CreateH1F("hWeffpMetCut",     "", 3, 0, 3);
  hWeffJetVeto     = CreateH1F("hWeffJetVeto",     "", 3, 0, 3);
  hWeffDeltaPhiJet = CreateH1F("hWeffDeltaPhiJet", "", 3, 0, 3);
  hWeffSoftMuVeto  = CreateH1F("hWeffSoftMuVeto",  "", 3, 0, 3);
  hWeffExtraLepton = CreateH1F("hWeffExtraLepton", "", 3, 0, 3);
  hWeffPtll        = CreateH1F("hWeffPtll",        "", 3, 0, 3);
  hWeffTopTagging  = CreateH1F("hWeffTopTagging",  "", 3, 0, 3);

  // WW level histograms     
  //----------------------------------------------------------------------------  

  hPtLepton1WWLevel       = CreateH1F("hPtLepton1WWLevel",       "", 200, 0, 200);
  hPtLepton2WWLevel       = CreateH1F("hPtLepton2WWLevel",       "", 200, 0, 200);
  hPtDiLeptonWWLevel      = CreateH1F("hPtDiLeptonWWLevel",      "", 200, 0, 200);
  hMinvWWLevel            = CreateH1F("hMinvWWLevel",            "", 200, 0, 200);
  hMtWWLevel              = CreateH1F("hMtWWLevel",              "", 250, 0, 250);
  hNJets30WWLevel         = CreateH1F("hNJetsPF30WWLevel",       "",  10, 0,  10);
  hpfMetWWLevel           = CreateH1F("hpfMetWWLevel",           "", 150, 0, 150);
  hppfMetWWLevel          = CreateH1F("hppfMetWWLevel",          "", 150, 0, 150);
  hchMetWWLevel           = CreateH1F("hchMetWWLevel",           "", 150, 0, 150);
  hpchMetWWLevel          = CreateH1F("hpchMetWWLevel",          "", 150, 0, 150);
  hpminMetWWLevel         = CreateH1F("hpminMetWWLevel",         "", 150, 0, 150);
  hDeltaRLeptonsWWLevel   = CreateH1F("hDeltaRLeptonsWWLevel",   "",  50, 0,   5);
  hDeltaPhiLeptonsWWLevel = CreateH1F("hDeltaPhiLeptonsWWLevel", "",  32, 0, 3.2);
  hDPhiPtllJetWWLevel     = CreateH1F("hDPhiPtllJetWWLevel",     "",  32, 0, 3.2);

  h_WWLevel_TightFailEvents = CreateH1F("h_WWLevel_TightFailEvents", "", 3, 0 , 3);
  h_WWLevel_TightTightEvents = CreateH1F("h_WWLevel_TightTightEvents", "", 3, 0 , 3);
  h_WWLevel_TightLooseEvents = CreateH1F("h_WWLevel_TightLooseEvents", "", 3, 0 , 3);


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

  h_TwoLeptons_TightFailEvents = CreateH1F("h_TwoLeptons_TightFailEvents", "", 3, 0 , 3);
  h_TwoLeptons_TightTightEvents = CreateH1F("h_TwoLeptons_TightTightEvents", "", 3, 0 , 3);
  h_TwoLeptons_TightLooseEvents = CreateH1F("h_TwoLeptons_TightLooseEvents", "", 3, 0 , 3);


  // Dilepton level differential histograms 
  //----------------------------------------------------------------------------   

  //differential in dilepton (eta, pt, ptll, mll, dphill), jets (jetpt1, njet), transverse mass? 2-jet bin?   

  Double_t pt1bins[8] = {20,40,60,80,100,125,150,200};
  Double_t eta1bins[9] = {-2.4,-1.8,-1.2,-0.6,0,0.6,1.2,1.8,2.4};
  Double_t mllbins[9] = {20,40,60,80,100,125,150,175,200};
  Double_t dphibins[13] = {0,0.25,0.5,0.75,1,1.25,1.5,1.75,2,2.25,2.5,2.75,3};
  Double_t ptllbins[8] = {30,40,50,60,70,85,120,150};

  hPtLepton1DilepLevel_Diff  = CreateH1F("hPtLepton1DilepLevel_Diff",       "", 7, pt1bins);
  hEtaLepton1DilepLevel_Diff  = CreateH1F("hEtaLepton1DilepLevel_Diff",       "", 8, eta1bins);
  hDileptonDilepLevel_Diff  = CreateH1F("hDileptonDilepLevel_Diff",       "", 7, ptllbins);
  hMinvDilepLevel_Diff  = CreateH1F("hMinvDilepLevel_Diff",       "", 8, mllbins);
  hDeltaPhiDilepLevel_Diff  = CreateH1F("hDeltaPhiDilepLevel_Diff",       "", 12, dphibins);


  // WW level differential histograms    
  //----------------------------------------------------------------------------

  hPtLepton1WWLevel_Diff  = CreateH1F("hPtLepton1WWLevel_Diff",       "",7, pt1bins);
  hEtaLepton1WWLevel_Diff  = CreateH1F("hEtaLepton1WWLevel_Diff",       "", 8, eta1bins);
  hDileptonWWLevel_Diff  = CreateH1F("hDileptonWWLevel_Diff",       "", 7, ptllbins);
  hMinvWWLevel_Diff  = CreateH1F("hMinvWWLevel_Diff",       "", 8, mllbins);
  hDeltaPhiWWLevel_Diff  = CreateH1F("hDeltaPhiWWLevel_Diff",       "", 12, dphibins);




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
    hNinZevents     [nC] = CreateH1F(Form("hNinZevents%.1i",      nC+1 ), "",   3, 0,   3);
    hNoutZevents    [nC] = CreateH1F(Form("hNoutZevents%.1i",     nC+1 ), "",   3, 0,   3);
    hNinLooseZevents[nC] = CreateH1F(Form("hNinLooseZevents%.1i", nC+1 ), "",   3, 0,   3);
    hMassInZevents  [nC] = CreateH1F(Form("hMassInZevents%.1i",   nC+1 ), "", 200, 0, 200);
    hMassOutZevents [nC] = CreateH1F(Form("hMassOutZevents%.1i",  nC+1 ), "", 200, 0, 200);
  
}

    // Data-driven methods: Top                                                                                                                                                      
    //----------------------------------------------------------------------------                                                                                                   
    hTopTaggedEvents            = CreateH1F("hTopTaggedEvents",            "", 3, 0, 3);
    hNTopControlRegion          = CreateH1F("hNTopControlRegion",          "", 3, 0, 3);
    hNTopTaggedTopControlRegion = CreateH1F("hNTopTaggedTopControlRegion", "", 3, 0, 3);

    hbTagDisTopTaggedEvents            = CreateH1F("hbTagDisTopTaggedEvents",            "", 300, -10, 20);
    hbTagDisNTopControlRegion          = CreateH1F("hbTagDisNTopControlRegion",          "", 300, -10, 20);
    hbTagDisNTopTaggedTopControlRegion = CreateH1F("hbTagDisNTopTaggedTopControlRegion", "", 300, -10, 20);


  }



void test::InsideLoop() {
 
  TString TheSample = GetInputParameters()->TheNamedString("theSample");
  h_n_PV -> Fill(zveto);

 // The InsideLoop() function is called for each entry in the tree to be processed  




    Double_t efficiencyW = effW * triggW;
    Double_t totalW      = -999;

    
    Double_t mybaseW = 1;

    Float_t nllW = 1;


    if ( TheSample.Contains("_nnll") ) {
      if ( TheSample.Contains("_pow_") ) mybaseW = 5984.0/999864; 
      if ( TheSample.Contains("_mad_") ) mybaseW = 5984.0/1933235;
      if ( TheSample.Contains("_mcnlo_") ) mybaseW = 5984.0/539594;
    }
    

    //5812.3/539594; // mcnlo 
    //5812.3/1933235; // madgraph (1933232)
    // 5812.3/999864; // powheg (999860)

    if (TheSample.Contains("Data"))
      totalW = 1.0;

    else if (TheSample.Contains("WJetsFakes")) { 
      totalW = fakeW;
    
    
      
      if ( dataset >= 82 && dataset <= 86 ) 
	totalW  = fakeW * (1 + 0.5 * (dataset >= 82 && dataset <= 84)) * baseW * puW * effW * triggW * Luminosity;
    
    }


    else if (TheSample.Contains("_nnll") ) {
      efficiencyW = puW * effW * triggW * nllW;
      totalW      = (1 + 0.5 * (dataset >= 82 && dataset <= 84)) * mybaseW * efficiencyW * Luminosity;
    } 

    else {
      efficiencyW = puW * effW * triggW ;
      totalW      = (1 + 0.5 * (dataset >= 82 && dataset <= 84)) * baseW * efficiencyW * Luminosity;
    }

    //cout <<mybaseW << endl;

    //if (TheSample.Contains("DYtautau") && dataset != 86) totalW *= (19.4/12.1);
    
    
    // Help variables
    //--------------------------------------------------------------------------
    //   Int_t dphiv = ((njet <= 1 && dphiveto) || (njet > 1 && dphilljetjet < 165.*TMath::DegToRad()));

    Int_t dphiv = (njet <= 1 || (njet > 1 && dphilljetjet < 165.*TMath::DegToRad()));

    Float_t metvar = (njet <= 1) ? mpmet : pfmet;

    Float_t jetbin = njet;

    Float_t dyMVA = ( !sameflav || ( (njet!=0 || dymva1>0.88) && (njet!=1 || dymva1>0.84) && ( njet==0 || njet==1 || (pfmet > 45.0)) ) );

	    
    if (njet == 3) jetbin = 2;  


    Int_t vbfsel = ((jetpt3 <= 30 || !(jetpt3 > 30 && ((jeteta1-jeteta3 > 0 && jeteta2-jeteta3 < 0) || (jeteta2-jeteta3 > 0 && jeteta1-jeteta3 < 0)))));

    
    // The selection begins here
    //--------------------------------------------------------------------------
 
    if (TheSample != "DY" || mctruth != 2 || (flavorChannel != "EMu" || flavorChannel != "MuE"))
      if (run != 201191)
	if (trigger == 1)
	  if (pt2 >= 20) //increase the pt of the leptons to further reduce Wjets 
	    if (pt1 >= 20) //increase the pt of the leptons to further reduce Wjets 
	      if (ch1*ch2 < 0)
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


      if( pfmet > 20 && mpmet > 20 && mll > 12 && nextra == 0 && (dphiv || !sameflav)  && ptll>30 && (!sameflav || ptll>45) ) { 



	//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	//
	// Z+jets
	//
	//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

	//if (dphiv && jetbin == jetChannel && bveto_mu && (bveto_ip && (nbjettche == 0 || njet > 3))) { OLD CUTS

	if( jetbin == jetChannel && bveto_mu && bveto_ip==1 && nbjettche==0 && ptll > 45) { 


	  // Loop over the metvar bins
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

	if ( zveto && dyMVA ) {

	  
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
	if (zveto && dyMVA) {

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
	  
	  if ( pass2012ICHEP1 && pass2012ICHEP2) 
	    h_TwoLeptons_TightTightEvents ->Fill(1, totalW); 
	  

	  if ( pass2012ICHEP1 || pass2012ICHEP2) 
	    h_TwoLeptons_TightLooseEvents ->Fill(1, totalW); 
	


    if (nextra == 0) {
			
	    hWExtraLepton->Fill(1, totalW);
	    hWeffExtraLepton->Fill(1, efficiencyW);
      
	    if (pfmet > 20 ) { // removed for differential xsec
	   
	      hWMetCut->Fill(1, totalW);
	      hWeffMetCut->Fill(1, efficiencyW);
	
	      if (mll > 12) {
	    
		hWLowMinv->Fill(1, totalW);
		hWeffLowMinv->Fill(1, efficiencyW);

		if (zveto==1 || !sameflav) {
		  
		  hWZVeto->Fill(1, totalW); 
		  hWeffZVeto->Fill(1, efficiencyW); 

		  if (mpmet > 20  && dyMVA ) { 

		
		    hWpMetCut->Fill(1, totalW);
		    hWeffpMetCut->Fill(1, efficiencyW);		  
		    

		    if (dphiv || !sameflav) {


		      hWDeltaPhiJet->Fill(1, totalW);
		      hWeffDeltaPhiJet->Fill(1, efficiencyW);


		      if (bveto_mu) { //--> third soft lepton with pt> 3GeV 
		      
			hWSoftMuVeto->Fill(1, totalW); 
			hWeffSoftMuVeto->Fill(1, efficiencyW);
		     

			 
			if ( ptll>30 && (!sameflav || ptll>45) ) {
			    
			  hWPtll->Fill(1, totalW);
			  hWeffPtll->Fill(1, efficiencyW);
			  

			  if (jetbin == jetChannel) {
			    
			    hWJetVeto->Fill(1, totalW);
			    hWeffJetVeto->Fill(1, efficiencyW);
			    

			    if (bveto_ip==1 &&  nbjettche==0) {
			      
			      hWTopTagging->Fill(1, totalW);
			      hWeffTopTagging->Fill(1, efficiencyW);

		       
			
			  
			  hPtLepton1WWLevel      ->Fill(pt1,       totalW);
			  hPtLepton2WWLevel      ->Fill(pt2,       totalW);
			  hPtDiLeptonWWLevel     ->Fill(ptll,      totalW);
			  hMinvWWLevel           ->Fill(mll,       totalW);
			  hMtWWLevel             ->Fill(mth,       totalW);
			  hNJets30WWLevel        ->Fill(njet,      totalW);
			  hpfMetWWLevel          ->Fill(pfmet,     totalW);
			  hppfMetWWLevel         ->Fill(ppfmet,    totalW);
			  hchMetWWLevel          ->Fill(chmet,     totalW);
			  hpchMetWWLevel         ->Fill(pchmet,    totalW);
			  hpminMetWWLevel        ->Fill(mpmet,     totalW);
			  hDeltaRLeptonsWWLevel  ->Fill(drll,      totalW);
			  hDeltaPhiLeptonsWWLevel->Fill(dphill,    totalW);
			  hDPhiPtllJetWWLevel    ->Fill(dphilljet, totalW);

			  hPtLepton1WWLevel_Diff  ->Fill(pt1,       totalW);
			  hEtaLepton1WWLevel_Diff ->Fill(eta1,      totalW);
			  hDileptonWWLevel_Diff   ->Fill(ptll,      totalW);
			  hMinvWWLevel_Diff       ->Fill(mll,       totalW);
			  hDeltaPhiWWLevel_Diff   ->Fill(dphill,    totalW);


			  //if ( (pass2012ICHEP1 && !pass2012ICHEP2) || (!pass2012ICHEP1 && pass2012ICHEP2))    
			    h_WWLevel_TightFailEvents ->Fill(1, totalW); 
			  
			   if ( (pass2012ICHEP1 && pass2012ICHEP2)) 
			     h_WWLevel_TightTightEvents ->Fill(1, totalW); 
			  
			    if ( pass2012ICHEP1 || pass2012ICHEP2) 
			      h_WWLevel_TightLooseEvents ->Fill(1, totalW); 
	


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

   // Counting histograms                                                                  //----------------------------------------------------------------------------       

  hWTrigger     = ((TH1F*) FindOutput("hWTrigger"));
  hWMetCut      = ((TH1F*) FindOutput("hWMetCut"));
  hWLowMinv     = ((TH1F*) FindOutput("hWLowMinv"));
  hWZVeto       = ((TH1F*) FindOutput("hWZVeto"));
  hWpMetCut     = ((TH1F*) FindOutput("hWpMetCut"));
  hWJetVeto     = ((TH1F*) FindOutput("hWJetVeto"));
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
  hPtLepton1WWLevel       = ((TH1F*) FindOutput("hPtLepton1WWLevel"));
  hPtLepton2WWLevel       = ((TH1F*) FindOutput("hPtLepton2WWLevel"));
  hPtDiLeptonWWLevel      = ((TH1F*) FindOutput("hPtDiLeptonWWLevel"));
  hMinvWWLevel            = ((TH1F*) FindOutput("hMinvWWLevel"));
  hMtWWLevel              = ((TH1F*) FindOutput("hMtWWLevel"));
  hNJets30WWLevel         = ((TH1F*) FindOutput("hNJetsPF30WWLevel"));
  hpfMetWWLevel           = ((TH1F*) FindOutput("hpfMetWWLevel"));
  hppfMetWWLevel          = ((TH1F*) FindOutput("hppfMetWWLevel"));
  hchMetWWLevel           = ((TH1F*) FindOutput("hchMetWWLevel"));
  hpchMetWWLevel          = ((TH1F*) FindOutput("hpchMetWWLevel"));
  hpminMetWWLevel         = ((TH1F*) FindOutput("hpminMetWWLevel"));
  hDeltaRLeptonsWWLevel   = ((TH1F*) FindOutput("hDeltaRLeptonsWWLevel"));
  hDeltaPhiLeptonsWWLevel = ((TH1F*) FindOutput("hDeltaPhiLeptonsWWLevel"));
  hDPhiPtllJetWWLevel     = ((TH1F*) FindOutput("hDPhiPtllJetWWLevel"));

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

  h_TwoLeptons_TightFailEvents = ((TH1F*) FindOutput("h_TwoLeptons_TightFailEvents"));
  h_TwoLeptons_TightTightEvents = ((TH1F*) FindOutput("h_TwoLeptons_TightTightEvents"));
  h_TwoLeptons_TightLooseEvents = ((TH1F*) FindOutput("h_TwoLeptons_TightLooseEvents"));




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
