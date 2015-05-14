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

TStyle *tdrStyle = new TStyle("tdrStyle","Style for P-TDR");

TH1F *hss = new TH1F("hss","hss",10,0.,10.);
TH1F *hos = new TH1F("hos","hos",10,0.,10.);

int macroSSOS(TString plotName = "hWnJetsBvetoAfterHt"){

  //Defining OS files
  TFile *QCDos = new TFile("rootFiles/AllJet/OF/QCD.root","read");
  TFile *Topos = new TFile("rootFiles/AllJet/OF/Top.root","read");
  TFile *TTJetsos = new TFile("rootFiles/AllJet/OF/TTJets.root","read");
  TFile *WJetsos = new TFile("rootFiles/AllJet/OF/WJets.root","read");

  //Defining SS files
  TFile *QCDss = new TFile("rootFilesSS/AllJet/OF/QCD.root","read");
  TFile *Topss = new TFile("rootFilesSS/AllJet/OF/Top.root","read");
  TFile *TTJetsss = new TFile("rootFilesSS/AllJet/OF/TTJets.root","read");
  TFile *WJetsss = new TFile("rootFilesSS/AllJet/OF/WJets.root","read");

  //Filling plots
  /*
  hss -> Add((TH1F*)QCDss->Get(plotName));
  hos -> Add((TH1F*)QCDos->Get(plotName));
  
  hss -> Add((TH1F*)Topss->Get(plotName));
  hos -> Add((TH1F*)Topos->Get(plotName));

  hss -> Add((TH1F*)TTJetsss->Get(plotName));
  hos -> Add((TH1F*)TTJetsos->Get(plotName));
*/  
  hss -> Add((TH1F*)WJetsss->Get(plotName));
  hos -> Add((TH1F*)WJetsos->Get(plotName));

  TCanvas *c1 = new TCanvas("c1","c1",600,600);
  c1->cd();
  TPad* pad1 = new TPad("pad1", "pad1", 0., 0., 1.0, 1.0);
  pad1->SetLeftMargin(0.14);
  pad1->SetBottomMargin(0.18);
  pad1->Draw();
  pad1->cd();

  tdrStyle->SetTitleFontSize(0.12);
  tdrStyle->SetTitleH(1); // Set the height of the title box
  tdrStyle->SetTitleW(1); // Set the width of the title box

  //hos->SetTitle("Number of Jets After B-Veto and Ht Cut \n in W + Jets MC Sample");
  
  hos->GetXaxis()->SetTitle("Number of Jets");
  
  hos->GetXaxis()->SetTitleSize(0.07);
  hos->GetXaxis()->SetTitleOffset(0.9);
  hos->GetXaxis()->SetLabelSize(0.05);
  hos->GetYaxis()->SetTitleSize(0.05);
  hos->GetYaxis()->SetLabelSize(0.05);
  hos->SetLineWidth(3);
  hss->SetLineWidth(3);
  hss->SetStats(0);
  hos->SetStats(0);
  hos->SetLineColor(kBlue);
  hss->SetLineColor(kRed);

  hos->Draw("");
  hss->Draw("same");

  TLegend* leg = new TLegend(0.50,0.79,0.70,0.89);
  leg->AddEntry(hss,"Same Sign","l");
  leg->AddEntry(hos,"Opposite Sign","l");
  leg->SetTextSize(0.05);
  leg->SetFillColor(kWhite);
  leg->SetLineColor(kWhite);
  leg->Draw();
  c1->Print("FakeRatio13TeV.pdf","pdf");

  for(int i = 0; i < hos->GetNbinsX(); ++i){
    if(hss->GetBinContent(i) != 0){
      cout<<hos->GetBinContent(i)/hss->GetBinContent(i)<<endl;
    }
  }
}
