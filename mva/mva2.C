//c++ `root-config --cflags --glibs` -o mva mva.cpp

#include<vector>
#include<iostream>
#include<cmath>

#include "TGraphErrors.h"
#include "TLegend.h"
#include "THStack.h"
#include "TChain.h"
#include "TString.h"
#include "TH1.h"
#include "TH2.h"
#include "TCanvas.h"
#include "TNtuple.h"
#include "TLorentzVector.h"
#include "TMVA/Types.h"
#include "TMVA/Factory.h"
#include "TMVA/Reader.h"
#include "TFile.h"

using namespace std ;

Float_t fullpmet;
Float_t trkpmet;
Float_t ratioMet;
Float_t ptll;
Float_t mth;
Float_t jetpt1;
Float_t ptWW;
Float_t dphilljet;
Float_t dphillmet;
Float_t dphijet1met;
Float_t nvtx;

Float_t baseW;       

int readDataset (TString datasetBaseName)
{
  TFile *f = new TFile(datasetBaseName,"read");
  TTree* ch = (TTree*) f -> Get("nt");
  Float_t htotalW;
  ch -> SetBranchAddress("baseW",&baseW);
  ch -> GetEntry(10);
  return htotalW;
}

void test_train(TString signalName = "WW",
		TString bkgName = "DY")
{
  TFile *outFile = new TFile("myAnalysisFile.root","RECREATE");
  
  TMVA::Factory *factory = new TMVA::Factory(signalName, outFile,"");
  
  TString directory = "../rootFiles/SF/MediumIDTighterIP/";
  //signalName = directory + signalName;
  
  //defining WW signal
  TFile *MySignalFile = new TFile("../rootFiles/SF/MediumIDTighterIP/WW.root","READ");
  TTree* sigTree = (TTree*)MySignalFile->Get("nt");
  factory->AddSignalTree(sigTree,1);
  
  //defining DY background
  TFile *MyBkgFile = new TFile("../rootFiles/SF/MediumIDTighterIP/DY.root","READ");
  TTree* bkgTree = (TTree*)MyBkgFile->Get("nt");
  factory->AddBackgroundTree(bkgTree,1);

  factory->SetWeightExpression("baseW");

  //************************************ FACTORY  
  
  factory->AddVariable("fullpmet");
  factory->AddVariable("trkpmet");
  factory->AddVariable("ratioMet");
  factory->AddVariable("ptll");
  factory->AddVariable("mth");
  factory->AddVariable("jetpt1");
  factory->AddVariable("ptWW");
  factory->AddVariable("dphilljet");
  factory->AddVariable("dphillmet");
  factory->AddVariable("dphijet1met");
  factory->AddVariable("nvtx");

  factory->PrepareTrainingAndTestTree("",500,500,500,500);
  cout<<"I've prepared trees"<<endl;
  //factory->BookMethod(TMVA::Types::kFisher, "Fisher","");
  factory->BookMethod(TMVA::Types::kBDT, "BDT","");
  
  cout<<"I've booked method"<<endl;
  factory->TrainAllMethods();
  factory->TestAllMethods();
  cout<<"I've tested all methods"<<endl;
  factory->EvaluateAllMethods();
  cout<<"I've evaluated all methods"<<endl;
  
}

void allBranches(TTree* inTree){
  
  inTree->SetBranchAddress("fullpmet",&fullpmet);
  inTree->SetBranchAddress("trkpmet",&trkpmet);
  inTree->SetBranchAddress("ratioMet",&ratioMet);
  inTree->SetBranchAddress("ptll",&ptll);
  inTree->SetBranchAddress("mth",&mth);
  inTree->SetBranchAddress("jetpt1",&jetpt1);
  inTree->SetBranchAddress("ptWW",&ptWW);
  inTree->SetBranchAddress("dphilljet",&dphilljet);
  inTree->SetBranchAddress("dphillmet",&dphillmet);
  inTree->SetBranchAddress("dphijet1met",&dphijet1met);
  inTree->SetBranchAddress("nvtx",&nvtx);

  inTree->SetBranchAddress("baseW",&baseW);        
}

void createOutput(TTree* outTree){

  outTree->Branch("fullpmet",&fullpmet);
  outTree->Branch("trkpmet",&trkpmet);
  outTree->Branch("ratioMet",&ratioMet);
  outTree->Branch("ptll",&ptll);
  outTree->Branch("mth",&mth);
  outTree->Branch("jetpt1",&jetpt1);
  outTree->Branch("ptWW",&ptWW);
  outTree->Branch("dphilljet",&dphilljet);
  outTree->Branch("dphillmet",&dphillmet);
  outTree->Branch("dphijet1met",&dphijet1met);
  outTree->Branch("nvtx",&nvtx);

  outTree->Branch("baseW",&baseW);        

}

void read(TString sampleName = "Dark10")
{
  //calling the reader of the MVA analysis
  TMVA::Reader* reader = new TMVA::Reader("");
  
  reader->AddVariable("fullpmet",&fullpmet);
  reader->AddVariable("trkpmet",&trkpmet);
  reader->AddVariable("ratioMet",&ratioMet);
  reader->AddVariable("ptll",&ptll);
  reader->AddVariable("mth",&mth);
  reader->AddVariable("jetpt1",&jetpt1);
  reader->AddVariable("ptWW",&ptWW);
  reader->AddVariable("dphilljet",&dphilljet);
  reader->AddVariable("dphillmet",&dphillmet);
  reader->AddVariable("dphijet1met",&dphijet1met);
  reader->AddVariable("nvtx",&nvtx);
  
  //reader->BookMVA("Fisher", "weights/MVAnalysis_Fisher.weights.xml");
  reader->BookMVA("BDT", "weights/" + sampleName + "_BDT.weights.xml");
  
  //calling Dark Matter File and Tree
  TFile *MyDarkMatterFile = new TFile("../rootFiles/AllJet/OF/" + sampleName + ".root","READ");
  TTree* DarkMatterTree = (TTree*)MyDarkMatterFile->Get("nt");

  DarkMatterTree->SetBranchAddress("fullpmet",&fullpmet);
  DarkMatterTree->SetBranchAddress("trkpmet",&trkpmet);
  DarkMatterTree->SetBranchAddress("ratioMet",&ratioMet);
  DarkMatterTree->SetBranchAddress("ptll",&ptll);
  DarkMatterTree->SetBranchAddress("mth",&mth);
  DarkMatterTree->SetBranchAddress("jetpt1",&jetpt1);
  DarkMatterTree->SetBranchAddress("ptWW",&ptWW);
  DarkMatterTree->SetBranchAddress("dphilljet",&dphilljet);
  DarkMatterTree->SetBranchAddress("dphillmet",&dphillmet);
  DarkMatterTree->SetBranchAddress("dphijet1met",&dphijet1met);
  DarkMatterTree->SetBranchAddress("nvtx",&nvtx);

  //calling ZH File and Tree
  TFile *MyZHFile = new TFile("../rootFiles/AllJet/OF/ZH.root","READ");
  TTree* ZHTree = (TTree*)MyZHFile->Get("nt");

  ZHTree->SetBranchAddress("fullpmet",&fullpmet);
  ZHTree->SetBranchAddress("trkpmet",&trkpmet);
  ZHTree->SetBranchAddress("ratioMet",&ratioMet);
  ZHTree->SetBranchAddress("ptll",&ptll);
  ZHTree->SetBranchAddress("mth",&mth);
  ZHTree->SetBranchAddress("jetpt1",&jetpt1);
  ZHTree->SetBranchAddress("ptWW",&ptWW);
  ZHTree->SetBranchAddress("dphilljet",&dphilljet);
  ZHTree->SetBranchAddress("dphillmet",&dphillmet);
  ZHTree->SetBranchAddress("dphijet1met",&dphijet1met);
  ZHTree->SetBranchAddress("nvtx",&nvtx);

  //calling HWW File and Tree
  TFile *MyHWWFile = new TFile("../rootFiles/AllJet/OF/HWW.root","READ");
  TTree* HWWTree = (TTree*)MyHWWFile->Get("nt");

  HWWTree->SetBranchAddress("fullpmet",&fullpmet);
  HWWTree->SetBranchAddress("trkpmet",&trkpmet);
  HWWTree->SetBranchAddress("ratioMet",&ratioMet);
  HWWTree->SetBranchAddress("ptll",&ptll);
  HWWTree->SetBranchAddress("mth",&mth);
  HWWTree->SetBranchAddress("jetpt1",&jetpt1);
  HWWTree->SetBranchAddress("ptWW",&ptWW);
  HWWTree->SetBranchAddress("dphilljet",&dphilljet);
  HWWTree->SetBranchAddress("dphillmet",&dphillmet);
  HWWTree->SetBranchAddress("dphijet1met",&dphijet1met);
  HWWTree->SetBranchAddress("nvtx",&nvtx);

  //calling WW File and Tree 
  TFile *MyWWFile = new TFile("../rootFiles/AllJet/OF/WW.root","READ");
  TTree* WWTree = (TTree*)MyWWFile->Get("nt");

  WWTree->SetBranchAddress("fullpmet",&fullpmet);
  WWTree->SetBranchAddress("trkpmet",&trkpmet);
  WWTree->SetBranchAddress("ratioMet",&ratioMet);
  WWTree->SetBranchAddress("ptll",&ptll);
  WWTree->SetBranchAddress("mth",&mth);
  WWTree->SetBranchAddress("jetpt1",&jetpt1);
  WWTree->SetBranchAddress("ptWW",&ptWW);
  WWTree->SetBranchAddress("dphilljet",&dphilljet);
  WWTree->SetBranchAddress("dphillmet",&dphillmet);
  WWTree->SetBranchAddress("dphijet1met",&dphijet1met);
  WWTree->SetBranchAddress("nvtx",&nvtx);

  //creating output Dark Matter Tree
  TTree *outDarkMatterTree = new TTree ("Dark","Dark");
  outDarkMatterTree -> SetDirectory(0);
  createOutput(outDarkMatterTree);
  
  //creating output ZH Tree
  TTree *outZHTree = new TTree ("ZH","ZH");
  outZHTree -> SetDirectory(0);
  createOutput(outZHTree);
  
  //creating output HWW Tree
  TTree *outHWWTree = new TTree ("HWW","HWW");
  outHWWTree -> SetDirectory(0);
  createOutput(outHWWTree);
  
  //creating output WW Tree
  TTree *outWWTree = new TTree ("WW","WW");
  outWWTree -> SetDirectory(0);
  createOutput(outWWTree);

  Float_t cutpt1        = 20;
  Float_t cutpt2        = 20;
  Float_t cutptll       = 0;
  Float_t cutmll        = 0;
  Float_t cutmth        = 0;
  Float_t cutpfType1Met = 0;
  Float_t cutdrll       = 0;
  Float_t cutdphill     = 0;
  Float_t cutdphilljet  = 0;
  Float_t cutdphillmet  = 0;
  Float_t cuttrkMet     = 0;
  Float_t cutMt1        = 0;
  Float_t cutMt2        = 0;
  Float_t cutmpmet      = 0;
  Float_t cutMc         = 0;
  Float_t cutptWW       = 0;
  Float_t cutHt         = 0;
  Float_t cutValue      = 0.15;

  Float_t value         = 0;

  //applying selections on Dark Matter sample
  Int_t contSig = 0;
  for (int i = 0; i < DarkMatterTree->GetEntries(); ++i){
    if (i == 0) cout<<"Dark Matter Entries : "<<DarkMatterTree->GetEntries()<<endl;
    DarkMatterTree->GetEntry(i);

    if (pt1        < cutpt1)        continue;
    if (pt2        < cutpt2)        continue;
    if (ptll       < cutptll)       continue;
    if (mll        < cutmll)        continue;
    if (mth        < cutmth)        continue;
    if (pfType1Met < cutpfType1Met) continue;
    if (drll       < cutdrll)       continue;
    if (dphill     < cutdphill)     continue;
    if (dphilljet  < cutdphilljet)  continue;
    if (dphillmet  < cutdphillmet)  continue;
    if (trkMet     < cuttrkMet)     continue;
    if (Mt1        < cutMt1)        continue;
    if (Mt2        < cutMt2)        continue;
    if (mpmet      < cutmpmet)      continue;
    if (Mc         < cutMc)         continue;
    if (ptWW       < cutptWW)       continue;
    if (Ht         < cutHt)         continue;

    value = reader->EvaluateMVA("BDT");
    if(value       < cutValue)      continue;

    ++contSig;
    outDarkMatterTree->Fill();
  }
  cout<<"Signal survived: "<<contSig<<endl;  

  //applying selections on ZH sample
  Int_t contZH = 0;
  for (int i = 0; i < ZHTree->GetEntries(); ++i){
    if (i == 0) cout<<"ZH Entries : "<<ZHTree->GetEntries()<<endl;
    ZHTree->GetEntry(i);

    if (pt1        < cutpt1)        continue;
    if (pt2        < cutpt2)        continue;
    if (ptll       < cutptll)       continue;
    if (mll        < cutmll)        continue;
    if (mth        < cutmth)        continue;
    if (pfType1Met < cutpfType1Met) continue;
    if (drll       < cutdrll)       continue;
    if (dphill     < cutdphill)     continue;
    if (dphilljet  < cutdphilljet)  continue;
    if (dphillmet  < cutdphillmet)  continue;
    if (trkMet     < cuttrkMet)     continue;
    if (Mt1        < cutMt1)        continue;
    if (Mt2        < cutMt2)        continue;
    if (mpmet      < cutmpmet)      continue;
    if (Mc         < cutMc)         continue;
    if (ptWW       < cutptWW)       continue;
    if (Ht         < cutHt)         continue;

    value = reader->EvaluateMVA("BDT");
    if(value       < cutValue)      continue;

    ++contZH;
    outZHTree->Fill();
  }
    
  cout<<"ZH survived: "<<contZH<<endl;  

  //applying selections on HWW sample
  Int_t contHWW = 0;
  for (int i = 0; i < HWWTree->GetEntries(); ++i){
    if (i == 0) cout<<"HWW Entries : "<<HWWTree->GetEntries()<<endl;
    HWWTree->GetEntry(i);

    if (pt1        < cutpt1)        continue;
    if (pt2        < cutpt2)        continue;
    if (ptll       < cutptll)       continue;
    if (mll        < cutmll)        continue;
    if (mth        < cutmth)        continue;
    if (pfType1Met < cutpfType1Met) continue;
    if (drll       < cutdrll)       continue;
    if (dphill     < cutdphill)     continue;
    if (dphilljet  < cutdphilljet)  continue;
    if (dphillmet  < cutdphillmet)  continue;
    if (trkMet     < cuttrkMet)     continue;
    if (Mt1        < cutMt1)        continue;
    if (Mt2        < cutMt2)        continue;
    if (mpmet      < cutmpmet)      continue;
    if (Mc         < cutMc)         continue;
    if (ptWW       < cutptWW)       continue;
    if (Ht         < cutHt)         continue;

    value = reader->EvaluateMVA("BDT");
    if(value       < cutValue)      continue;

    ++contHWW;
    outHWWTree->Fill();
  }

  cout<<"HWW survived: "<<contHWW<<endl;
  
  //applying selections on WW sample
  Int_t contWW = 0;
  for (int i = 0; i < WWTree->GetEntries(); ++i){
    if (i == 0) cout<<"WW Entries : "<<WWTree->GetEntries()<<endl;
    WWTree->GetEntry(i);

    if (pt1        < cutpt1)        continue;
    if (pt2        < cutpt2)        continue;
    if (ptll       < cutptll)       continue;
    if (mll        < cutmll)        continue;
    if (mth        < cutmth)        continue;
    if (pfType1Met < cutpfType1Met) continue;
    if (drll       < cutdrll)       continue;
    if (dphill     < cutdphill)     continue;
    if (dphilljet  < cutdphilljet)  continue;
    if (dphillmet  < cutdphillmet)  continue;
    if (trkMet     < cuttrkMet)     continue;
    if (Mt1        < cutMt1)        continue;
    if (Mt2        < cutMt2)        continue;
    if (mpmet      < cutmpmet)      continue;
    if (Mc         < cutMc)         continue;
    if (ptWW       < cutptWW)       continue;
    if (Ht         < cutHt)         continue;

    value = reader->EvaluateMVA("BDT");
    if(value             < cutValue)           continue;
  
    ++contWW;
    outWWTree->Fill();
  }
  
  cout<<"WW survived: "<<contWW<<endl;

  //saving trees
  TFile *outMVA = new TFile("outMVA" + sampleName + ".root","RECREATE");
  outMVA -> cd();
  outDarkMatterTree -> Write();
  outZHTree -> Write();
  outHWWTree -> Write();
  outWWTree -> Write();
  outMVA -> Close();
}

void draw(){
  
TGraphErrors *gll = new TGraphErrors();
TGraphErrors *gjj = new TGraphErrors();

 for (int i = 0 ; i < 15 ; ++i){
   float content = h_lep_m_signal->GetBinContent(i) + h_lep_m_background->GetBinContent(i) + h_lep_m_wjets->GetBinContent(i);
   gll -> SetPoint(i,i*100-50.,content);
   float errorb = h_lep_m_background->GetBinContent(i)*0.08*h_lep_m_background->GetBinContent(i);
   float errors = h_lep_m_signal->GetBinContent(i)*0.08*h_lep_m_signal->GetBinContent(i);
   float errorw = 0.40 * h_lep_m_wjets->GetBinContent(i)*h_lep_m_wjets->GetBinContent(i);
   float pesi = h_lep_m_background->GetBinContent(i) + h_lep_m_signal->GetBinContent(i) + h_lep_m_wjets->GetBinContent(i);
   float error = sqrt(errorb*errorb + errorw*errorw + errors*errors) / pesi;
   gll -> SetPointError(i,100/2.,error);
 }

for (int i = 0 ; i < 10 ; ++i){
  float content = h_jet_m_signal->GetBinContent(i) + h_jet_m_background->GetBinContent(i) + h_jet_m_wjets->GetBinContent(i);
  gjj -> SetPoint(i,i*500-250.,content);	
  float errorb = h_jet_m_background->GetBinContent(i)*0.08*h_jet_m_background->GetBinContent(i);
  float errors = h_jet_m_signal->GetBinContent(i)*0.08*h_jet_m_signal->GetBinContent(i);
  float errorw = 0.40 * h_jet_m_wjets->GetBinContent(i)*h_jet_m_wjets->GetBinContent(i);
  float pesi = h_jet_m_background->GetBinContent(i) + h_jet_m_signal->GetBinContent(i) + h_jet_m_wjets->GetBinContent(i);	
  float error = sqrt(errorb*errorb + errorw*errorw + errors*errors) / pesi;
  gjj -> SetPointError(i,500/2.,error);
 }
 
 gll->SetFillStyle(3002);
 gjj->SetFillStyle(3002);
 
 gll -> SetFillColor(kBlack);
 gjj -> SetFillColor(kBlack);
 
 
 cout<<"signal tot = "<<scontTot<<endl;
 cout<<"signal mjj = "<<scontMjj<<endl;
 cout<<"signal deta = "<<scontDeta<<endl;
 cout<<"signal mll = "<<scontMll<<endl;
 cout<<"signal dy = "<<scontDY<<endl;
 cout<<"signal met = "<<scontMET<<"\n"<<endl;
 
 cout<<"background tot = "<<contTot<<endl;
 cout<<"background mjj = "<<contMjj<<endl;
 cout<<"background deta = "<<contDeta<<endl;
 cout<<"background mll = "<<contMll<<endl;
 cout<<"background dy = "<<contDY<<endl;
 cout<<"background met = "<<contMET<<"\n"<<endl;
 
 cout<<"wjets tot = "<<wcontTot<<endl;
 cout<<"wjets mjj = "<<wcontMjj<<endl;
 cout<<"wjets deta = "<<wcontDeta<<endl;
 cout<<"wjets mll = "<<wcontMll<<endl;
 cout<<"wjets dy = "<<wcontDY<<endl;
 cout<<"wjets met = "<<wcontMET<<"\n"<<endl;
 
 TH1F *axis = new TH1F ("mll Signal & Backgrounds","mll Signal & Backgrounds",15,0.,1500.);
 axis -> GetXaxis() -> SetRangeUser(0.,1500.);
 axis -> GetYaxis() -> SetRangeUser(0.,120.);
 axis -> GetXaxis () -> SetTitle("Leptons Invariant Mass [GeV]");
 axis -> GetYaxis () -> SetTitle("Counts per Bin");
 axis -> GetYaxis () -> SetTitleOffset(1.3);
 axis -> SetStats(0);
 
 TH1F *axis2 = new TH1F ("mjj Signal & Backgrounds","mjj Signal & Backgrounds",10,0.,5000.);
 axis2 -> GetXaxis() -> SetRangeUser(0.,5000.);
 axis2 -> GetYaxis() -> SetRangeUser(0.,105.);
 axis2 -> GetXaxis () -> SetTitle("Tag Jets Invariant Mass [GeV]");
 axis2 -> GetYaxis () -> SetTitle("Counts per Bin");
 axis2 -> GetYaxis () -> SetTitleOffset(1.3);
 axis2 -> SetStats(0);
 
 THStack* hs = new THStack("hs","mll Signal & Background");
 
 h_lep_m_signal -> SetFillColor(kGreen+1);
 h_lep_m_background -> SetFillColor(kRed+1);
 h_lep_m_wjets -> SetFillColor(kBlue+1);
 
 hs -> Add ( h_lep_m_background );
 hs -> Add ( h_lep_m_wjets );
 hs -> Add ( h_lep_m_signal );
 
 TCanvas *c1 = new TCanvas("c1","c1",600.,600.);
 c1 -> cd();
 axis -> Draw();
 //c1 -> DrawFrame(0.,0.,1500.,40.);
 hs -> Draw("same");
 hs -> GetXaxis () -> SetTitle("Leptons Invariant Mass [GeV]");
 hs -> GetYaxis () -> SetRangeUser(0.,40.);
 TLegend *leg = new TLegend(0.54,0.3,0.89,0.5);
 leg->SetHeader("Luminosity = 300fb^{-1} @ 13TeV");
 leg->AddEntry(h_lep_m_signal,"Signal","f");
 leg->AddEntry(h_lep_m_background,"QCD Background","f");
 leg->AddEntry(h_lep_m_wjets,"W + Jets Background","f");
 leg->SetFillColor(kWhite);
 leg->SetLineColor(kWhite);
 leg->Draw();
 gll->Draw("same,2");
 c1 -> Print("stack.C");
 c1 -> Print("stack.pdf","pdf");
 c1 -> Print("stack.png","png");
 
 TCanvas *c2 = new TCanvas("c2","c2",600.,600.);
 c2 -> cd();
 hs -> Draw("nostack");
 hs -> GetXaxis () -> SetTitle("Leptons Invariant Mass [GeV]");
 hs -> GetYaxis () -> SetRangeUser(0.,20.);
 c2 -> Print("nostack.C");
 c2 -> Print("nostack.pdf","pdf");
 
 THStack* hs2 = new THStack("hs2","mjj Signal & Background");

 h_jet_m_signal -> SetFillColor(kGreen+1);
 h_jet_m_background -> SetFillColor(kRed+1);
 h_jet_m_wjets -> SetFillColor(kBlue+1);
 
 hs2 -> Add ( h_jet_m_background );
 hs2 -> Add ( h_jet_m_wjets );
 hs2 -> Add ( h_jet_m_signal );
 
  TCanvas *c12 = new TCanvas("c12","c12",600.,600.);
 c12 -> cd();
 axis2 -> Draw();
 //c12 -> DrawFrame(0.,0.,5000.,35.);
 hs2 -> Draw("same");
 hs2 -> GetXaxis () -> SetTitle("Tag Jets Invariant Mass [GeV]");
 hs2 -> GetYaxis () -> SetRangeUser(0.,35.);
 TLegend *leg2 = new TLegend(0.54,0.3,0.89,0.5);
 leg2->SetHeader("Luminosity = 300fb^{-1} @ 13TeV");
 leg2->AddEntry(h_jet_m_signal,"Signal","f");
 leg2->AddEntry(h_jet_m_background,"QCD Background","f");
 leg2->AddEntry(h_jet_m_wjets,"W + Jets Background","f");
 leg2->SetFillColor(kWhite);
 leg2->SetLineColor(kWhite);
 leg2->Draw();
 gjj->Draw("same,2");
 c12 -> Print("stack_mjj.C");
 c12 -> Print("stack_mjj.pdf","pdf");
 c12 -> Print("stack_mjj.png","png");
 
 TCanvas *c22 = new TCanvas("c22","c22",600.,600.);
 c22 -> cd();
 hs2 -> Draw("nostack");
 hs2 -> GetXaxis () -> SetTitle("Tag Jet Invariant Mass [GeV]");
 hs2 -> GetYaxis () -> SetRangeUser(0.,20.);
 c22 -> Print("nostack_mjj.C");
 c22 -> Print("nostack_mjj.pdf","pdf");

 THStack* hs3 = new THStack("hs3","BDT Signal & Background");
 
 h_BDT_signal -> SetFillColor(kGreen+1);
 h_BDT_background -> SetFillColor(kRed+1);
 h_BDT_wjets -> SetFillColor(kBlue+1);
 
 hs3 -> Add ( h_BDT_background );
 hs3 -> Add ( h_BDT_wjets );
 hs3 -> Add ( h_BDT_signal );
 
 TCanvas *c13 = new TCanvas("c13","c13",600.,600.);
 c13 -> cd();
 hs3 -> Draw();
 hs3 -> GetXaxis () -> SetTitle("BDT discriminant");
 hs3 -> GetYaxis () -> SetRangeUser(0.,60.);
 TLegend *leg3 = new TLegend(0.54,0.3,0.89,0.5);
 leg3->SetHeader("Luminosity = 300fb^{-1} @ 13TeV");
 leg3->AddEntry(h_BDT_signal,"Signal","f");
 leg3->AddEntry(h_BDT_background,"QCD Background","f");
 leg3->AddEntry(h_BDT_wjets,"W + Jets Background","f");
 leg3->SetFillColor(kWhite);
 leg3->SetLineColor(kWhite);
 leg3->Draw();
 c13 -> Print("stack_BDT.C");
 c13 -> Print("stack_BDT.pdf","pdf");
 
 TCanvas *c23 = new TCanvas("c23","c23",600.,600.);
 c23 -> cd();
 hs3 -> Draw("nostack");
 hs3 -> GetXaxis () -> SetTitle("BDT discriminant");
 hs3 -> GetYaxis () -> SetRangeUser(0.,20.);
 c23 -> Print("nostack_BDT.C");
 c23 -> Print("nostack_BDT.pdf","pdf");

 THStack* hs4 = new THStack("hs4","BDT Signal & Backgrounds Before Cuts");

 h_BDT_signalNoCuts -> SetFillColor(kGreen+1);
 h_BDT_backgroundNoCuts -> SetFillColor(kRed+1);
 h_BDT_wjetsNoCuts -> SetFillColor(kBlue+1);

 hs4 -> Add ( h_BDT_backgroundNoCuts );
 hs4 -> Add ( h_BDT_wjetsNoCuts );
 hs4 -> Add ( h_BDT_signalNoCuts );

 TCanvas *c14 = new TCanvas("c14","c14",600.,600.);
 c14 -> cd();
 hs4 -> Draw();
 hs4 -> GetXaxis () -> SetTitle("BDT discriminant");
 hs4 -> GetYaxis () -> SetRangeUser(0.,60.);
 TLegend *leg4 = new TLegend(0.54,0.3,0.89,0.5);
 leg4->SetHeader("Luminosity = 300fb^{-1} @ 13TeV");
 leg4->AddEntry(h_BDT_signal,"Signal","f");
 leg4->AddEntry(h_BDT_background,"QCD Background","f");
 leg4->AddEntry(h_BDT_wjets,"W + Jets Background","f");
 leg4->SetFillColor(kWhite);
 leg4->SetLineColor(kWhite);
 leg4->Draw();
 c14 -> Print("stack_BDT_nocuts.C");
 c14 -> Print("stack_BDT_nocuts.pdf","pdf");
 
 TCanvas *c24 = new TCanvas("c24","c24",600.,600.);
 c24 -> cd();
 hs4 -> Draw("nostack");
 hs4 -> GetXaxis () -> SetTitle("BDT discriminant");
 hs4 -> GetYaxis () -> SetRangeUser(0.,20.);
 c24 -> Print("nostack_BDT_nocuts.C");
 c24 -> Print("nostack_BDT_nocuts.pdf","pdf");
 
 //TFile *outMVA = new TFile("outMVA_EWK_BDT_norm.root","RECREATE");
 //TFile *outMVA = new TFile("outMVA_EWK_Fisher.root","RECREATE");
 //TFile *outMVA = new TFile("outMVA_QCD_Fisher.root","RECREATE");
 
}


void mva2(){
  test_train();
  //read();
  //draw();
}

