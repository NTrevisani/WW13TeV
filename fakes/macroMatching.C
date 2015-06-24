//root -l 'macroMatching.C("hNjetsTwoLeptonsLevel")'

//#include "leptonClass.h"

//#include "latinoTree.h"
#include "TROOT.h"
#include "TCanvas.h"
#include "THStack.h"
#include "TH1F.h"
#include "TLatex.h"
#include "TSystem.h"
#include "TTree.h"
#include "TFile.h"
#include "TGraph.h"
#include "TList.h"
#include "TKey.h"
#include "TObject.h"
#include "TStyle.h"
//#include "latinoTree.h"
#include <vector>

TStyle *tdrStyle = new TStyle("tdrStyle","Style for P-TDR");
TH1F* hlepNumber = new TH1F("hlepNumber","hlepNumber", 20,0.,20.);
TH1F* hjetNumber = new TH1F("hjetNumber","hjetNumber", 20,0.,20.);

/*
//number of tight electrons with pt > 20GeV
int nel(std::vector<float> *lepPtvec,std::vector<float> *lepIdvec){
  int cont = 0;
  for(int q = 0; q < lepPtvec -> size(); ++q)
    if(fabs(lepIdvec) -> at(q) == 11)
      if(lepPtvec -> at(q) > 20)
	++cont;
  return cont;
}

//number of tight muons with pt > 20GeV
int nmu(std::vector<float> *lepPtvec,std::vector<float> *lepIdvec){
  int cont = 0;
  for(int q = 0; q < lepPtvec -> size(); ++q)
    if(fabs(lepIdvec) -> at(q) == 13)
      if(lepPtvec -> at(q) > 20)
	++cont;
  return cont;
}
*/
//number of tight leptons with pt > 20GeV
int nlep(std::vector<float> *lepPtvec,std::vector<int> *lepIdvec){
  int cont = 0;
  for(int q = 0; q < lepPtvec -> size(); ++q)
    if(fabs(lepIdvec -> at(q)) == 13 || fabs(lepIdvec -> at(q)) == 11)
      if(lepPtvec -> at(q) > 20)
	++cont;
  return cont;
}

int njet(std::vector<float> *jetPtvec){
  int cont = 0;
  for(int q = 0; q < jetPtvec -> size(); ++q)
    if(jetPtvec->at(q) > 30)
      ++cont;
  return cont;
}
/*
//number of prompt electrons with pt > 20GeV which pass the tight lepton requests
int npromptel(std::vector<float> *lepPtvec,std::vector<float> *lepIdvec){
  int cont = 0;
  for(int q = 0; q < lepMumId -> size(); ++q)
    if(fabs(lepIdvec) -> at(q) == 11)
      if(lepPtvec -> at(q) > 20)
	if(fabs(lepMumId -> at(q)) == 24 || lepMumId -> at(q) == 23)
	  if(lepMumStatus -> at(q) > 20 && lepMumStatus -> at(q) < 30)
	    ++cont;
  return cont;
}

//number of prompt muons with pt > 20GeV which pass the tight lepton requests
int npromptmu(std::vector<float> *lepPtvec,std::vector<float> *lepIdvec){
  int cont = 0;
  for(int q = 0; q < lepMumId -> size(); ++q)
    if(fabs(lepIdvec) -> at(q) == 13)
      if(lepPtvec -> at(q) > 20)
	if(fabs(lepMumId -> at(q)) == 24 || lepMumId -> at(q) == 23)
	  if(lepMumStatus -> at(q) > 20 && lepMumStatus -> at(q) < 30)
	    ++cont;
  return cont;
}

//number of prompt leptons with pt > 20GeV which pass the tight lepton requests
int nprompt(std::vector<float> *lepPtvec,std::vector<int> *lepIdvec){
  int cont = 0;
  for(int q = 0; q < lepMumId -> size(); ++q)
    if(fabs(lepIdvec -> at(q)) == 13 || fabs(lepIdvec -> at(q)) == 11)
      if(lepPtvec -> at(q) > 20)
	if(fabs(lepMumId -> at(q)) == 24 || lepMumId -> at(q) == 23)
	  if(lepMumStatus -> at(q) > 20 && lepMumStatus -> at(q) < 30)
	    ++cont;
  return cont;
}
*/

float dist(float phi1, float eta1, float phi2, float eta2){
  float dphi = fabs(phi1 - phi2);
  float deta = fabs(eta1 - eta2);
  if (dphi > TMath::Pi()) 
    dphi = 2*TMath::Pi() - dphi;
  return sqrt(dphi*dphi + deta*deta);
}

int macroMatching(TString plotName = "hWnJetsBvetoAfterHt"){

  //gInterpreter->LoadMacro("latinoTree.C+");

  //leptonClass lc;// = new leptonClass();

  TFile *f = new TFile("/gpfs/csic_projects/cms/trevisanin/newLatino/25ns/latino_WJetsToLNu.root");
  TTree *t = (TTree*) f -> Get("latino");

  //  latinoTree lt(t);

  std::vector<float> *lepton_pt  = 0;      t -> SetBranchAddress("std_vector_lepton_pt",&lepton_pt);
  std::vector<float> *lepton_eta = 0;      t -> SetBranchAddress("std_vector_lepton_eta",&lepton_eta);
  std::vector<float> *lepton_phi = 0;      t -> SetBranchAddress("std_vector_lepton_phi",&lepton_phi);
  std::vector<int>   *lepton_id  = 0;      t -> SetBranchAddress("std_vector_lepton_id",&lepton_id);

  std::vector<float> *leptonGen_pt  = 0;   t -> SetBranchAddress("std_vector_leptonGen_pt",&leptonGen_pt);
  std::vector<float> *leptonGen_eta = 0;   t -> SetBranchAddress("std_vector_leptonGen_eta",&leptonGen_eta);
  std::vector<float> *leptonGen_phi = 0;   t -> SetBranchAddress("std_vector_leptonGen_phi",&leptonGen_phi);
  std::vector<float> *leptonGen_mid = 0;   t -> SetBranchAddress("std_vector_leptonGen_mpid",&leptonGen_mid);
  std::vector<int>   *leptonGen_mst = 0;   t -> SetBranchAddress("std_vector_leptonGen_mstatus",&leptonGen_mst);

  std::vector<float> *jet_pt     = 0;      t -> SetBranchAddress("std_vector_jet_pt",&jet_pt);
  std::vector<float> *jet_eta    = 0;      t -> SetBranchAddress("std_vector_jet_eta",&jet_eta);
  std::vector<float> *jet_phi    = 0;      t -> SetBranchAddress("std_vector_jet_phi",&jet_phi);

  std::vector<float> *jetGen_pt  = 0;      t -> SetBranchAddress("std_vector_jetGen_pt",&jetGen_pt);
  std::vector<float> *jetGen_eta = 0;      t -> SetBranchAddress("std_vector_jetGen_eta",&jetGen_eta);
  std::vector<float> *jetGen_phi = 0;      t -> SetBranchAddress("std_vector_jetGen_phi",&jetGen_phi);

  std::vector<int>   *hadronFl   = 0;      t -> SetBranchAddress("std_vector_jet_HadronFlavour",&hadronFl);
  std::vector<int>   *partonFl   = 0;      t -> SetBranchAddress("std_vector_jet_PartonFlavour",&partonFl);

  TFile* out = new TFile("fakeNtu.root","recreate");
  out->cd();
  TTree* nt = new TTree("nt","nt");
  nt->SetDirectory(0);

  Float_t jetPt;        nt->Branch("jetPt",&jetPt,"jetPt");
  Float_t jetEta;       nt->Branch("jetEta",&jetEta,"jetEta");
  Float_t jetPhi;       nt->Branch("jetPhi",&jetPhi,"jetPhi");
  Float_t leptonPt;     nt->Branch("leptonPt",&leptonPt,"leptonPt");
  Float_t leptonPhi;    nt->Branch("leptonPhi",&leptonPhi,"leptonPhi");
  Float_t leptonEta;    nt->Branch("leptonEta",&leptonEta,"leptonEta");
  Float_t distance;     nt->Branch("distance",&distance,"distance");
  Int_t   ID;           nt->Branch("ID",&ID,"ID");
  Bool_t  sameSign;     nt->Branch("sameSign",&sameSign,"sameSign");
  Float_t jetGenPt;     nt->Branch("jetGenPt",&jetGenPt,"jetGenPt");
  Float_t jetGenEta;    nt->Branch("jetGenEta",&jetGenEta,"jetGenEta");
  Float_t jetGenPhi;    nt->Branch("jetGenPhi",&jetGenPhi,"jetGenPhi");
  Float_t leptonGenPt;  nt->Branch("leptonGenPt",&leptonGenPt,"leptonGenPt");
  Float_t leptonGenPhi; nt->Branch("leptonGenPhi",&leptonGenPhi,"leptonGenPhi");
  Float_t leptonGenEta; nt->Branch("leptonGenEta",&leptonGenEta,"leptonGenEta");
  Int_t   jetPoint;     nt->Branch("jetPoint",&jetPoint,"jetPoint");

  std::vector<int> point;
  std::vector<float> dist;
  std::vector<int> particleId;
  float distance_ = 0.;  

  for(int iEntry = 0; iEntry < t -> GetEntries(); ++iEntry){
    t->GetEntry(iEntry);
    if(iEntry == 0) cout<<"Total events: "<<t->GetEntries()<<endl;
    if(iEntry % 10000 == 0) cout<<"Reading Entry "<<iEntry<<endl;
    
    int lepNumber = nlep(lepton_pt, lepton_id);
    hlepNumber -> Fill(lepNumber);
    hjetNumber -> Fill(njet(jet_pt));
    
    if (lepton_pt->at(0) > 20 && lepton_pt->at(1) > 20){
      //first lepton is prompt, second lepton is non-prompt (looking at mother id and status)
      if(fabs(leptonGen_mid->at(0)) == 24 || fabs(leptonGen_mid->at(0)) == 23){
	//if(fabs(leptonGen_mid->at(1)) != 24 && fabs(leptonGen_mid->at(1)) != 23){
	//if(leptonGen_mst -> at(0) > 20 && leptonGen_mst -> at(0) < 30){
	if(leptonGen_mst -> at(1) < 20 || leptonGen_mst -> at(1) > 30){
	  //then the second lepton (1) is a fake -> Fill ntuple with its coordinate
	  //	  leptonPt  = lepton_pt ->at(1);
	  //	  leptonEta = lepton_eta->at(1);
	  //	  leptonPhi = lepton_phi->at(1);
	  //find the closest jet to the lepton
	  distance_ = 100000.;
	  int jetRef = -9999;
	  for(int i = 0; i < jet_pt->size(); ++i){
	    if(jet_pt->at(i) > 0){
	      float dumpDist = dist(leptonPhi, leptonEta, jet_phi->at(i), jet_eta->at(i));
	      if( dumpDist < distance_){
		distance_ = dumpDist;
		jetRef = i;
	      }
	    }
	  }
	  if (jetRef > 0){
	    jetPt  = jet_pt ->at(jetRef);
	    jetEta = jet_eta->at(jetRef);
	    jetPhi = jet_phi->at(jetRef);
	    leptonPt  = lepton_pt ->at(1);
	    leptonEta = lepton_eta->at(1);
	    leptonPhi = lepton_phi->at(1);	  
	    distance = distance_;
	    nt->Fill();
	  }
	  else
	    cout<<"No Matching!!";
	}
      }
      
      //first lepton is non-prompt, second lepton is prompt (looking at mother id and status)
      else if(fabs(leptonGen_mid->at(1)) == 24 || fabs(leptonGen_mid->at(1)) == 23){
	//if( fabs(leptonGen_mid->at(0)) != 24 && fabs(leptonGen_mid->at(0)) != 23){
	//if(leptonGen_mst -> at(1) > 20 && leptonGen_mst -> at(1) < 30)
	if(leptonGen_mst -> at(0) < 20 || leptonGen_mst -> at(0) > 30){
	  //then the first lepton (0) is a fake -> Fill ntuple with its coordinate
	  //	  leptonPt  = lepton_pt ->at(0);
	  //	  leptonEta = lepton_eta->at(0);
	  //	  leptonPhi = lepton_phi->at(0);
	  cout<<"ciao1"<<endl;
	  //find the closest jet to the lepton                                                                                               
	  distance_ = 10000.;
	  int jetRef = -9999;
	  for(int i = 0; i < jet_pt->size(); ++i){
	    if(jet_pt->at(i) > 0){
	      float dumpDist = dist(leptonPhi, leptonEta, jet_phi->at(i), jet_eta->at(i));
	      if( dumpDist < distance_){
		distance_ = dumpDist;
		jetRef = i;
	      }
	    }
	  }
	  if (jetRef > 0){
	    jetPt  = jet_pt ->at(jetRef);
	    jetEta = jet_eta->at(jetRef);
	    jetPhi = jet_phi->at(jetRef);
	    leptonPt  = lepton_pt ->at(0);
	    leptonEta = lepton_eta->at(0);
	    leptonPhi = lepton_phi->at(0);
	    distance = distance_;
	    nt->Fill();
	  }
	  else
	    cout<<"No Matching!!";
	}
      }
      
      else
	cout<<"This event has no prompt leptons!! :O"<<endl;
    }
  }
  nt->Write();
  hlepNumber->Write();
  hjetNumber->Write();
  out->Close();
  
  return 0;
}
