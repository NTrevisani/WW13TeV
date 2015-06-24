#include "TH2F.h"
#include "TTree.h"
#include "TPad.h"
#include "TTree.h"
#include "TStyle.h"

TStyle *tdrStyle = new TStyle("tdrStyle","Style for P-TDR");
int contGood = 0;
int contBad  = 0;

void macroPlot(){

  TFile *f = new TFile("fakeNtu.root");
  TTree *t = (TTree*) f -> Get("nt");

  Float_t leptonPt; t->SetBranchAddress("leptonPt",&leptonPt);
  Float_t jetPt;    t->SetBranchAddress("jetPt",&jetPt);
  Float_t jetPhi;   t->SetBranchAddress("jetPhi",&jetPhi);
  Float_t jetEta;   t->SetBranchAddress("jetEta",&jetEta);
  Float_t distance; t->SetBranchAddress("distance",&distance);

  TH2F *plot  = new TH2F("plot" ,"plot" ,100 ,0.,5. ,100,0.,5.);
  TH1F *hdist = new TH1F("hdist","hdist",100,0.,10);

  for(int iEntry = 0; iEntry < t -> GetEntries(); ++iEntry){
    t -> GetEntry(iEntry);
    if(iEntry % 100000 == 0) cout<<"Plotting Entry "<<iEntry<<endl;
    hdist->Fill(distance);
    if(distance < 0.5){ 
      plot->Fill(leptonPt/jetPt,fabs(jetEta));
      ++contGood;
    }
    else ++contBad;
}

  cout<<"contGood = "<<contGood<<endl;
  cout<<"contBad = "<<contBad<<endl;

  TCanvas *c1 = new TCanvas("c1","c1",600,600);
  c1->cd();
  TPad* pad1 = new TPad("pad1", "pad1", 0., 0., 1.0, 1.0);
  pad1->SetLeftMargin(0.14);
  pad1->SetRightMargin(0.14);
  pad1->SetBottomMargin(0.14);
  pad1->Draw();
  pad1->cd();

  tdrStyle->SetTitleFontSize(0.12);
  tdrStyle->SetTitleH(1); // Set the height of the title box
  tdrStyle->SetTitleW(1); // Set the width of the title box

  plot->SetTitle("Fake Rate");
  plot->SetStats(0);
  plot->GetXaxis()->SetTitle("p_{T}^{lep} / p_{T}^{jet}");
  plot->GetXaxis()->SetTitleSize(0.05);
  plot->GetYaxis()->SetTitle("|#eta_{jet}|");
  plot->GetYaxis()->SetTitleSize(0.05);

  plot->Draw("colz");
  c1->Print("FakeRatio13TeV.pdf","pdf");
  hdist->Draw();
}

