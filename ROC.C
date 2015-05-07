#include "TROOT.h"
#include "TCanvas.h"
#include "THStack.h"
#include "TH1F.h"
#include "TLatex.h"
#include "TSystem.h"
#include "TTree.h"
#include "TFile.h"
#include "TGraph.h"


int ROC(){
  //signal files
  TFile *fsig1 = new TFile("rootFiles/AllJet/OF/WWTo2L2Nu_pow_nnll.root","read");
  TFile *fsig2 = new TFile("rootFiles/AllJet/OF/ggWWto2L.root","read");

  //all the background files
  TFile *fbkg1 = new TFile("rootFiles/AllJet/OF/TTbar.root","read");
  TFile *fbkg2 = new TFile("rootFiles/AllJet/OF/TW.root","read");
  TFile *fbkg3 = new TFile("rootFiles/AllJet/OF/WZ.root","read");
  TFile *fbkg4 = new TFile("rootFiles/AllJet/OF/ZZ.root","read");
  TFile *fbkg5 = new TFile("rootFiles/AllJet/OF/WgammaNoStar.root","read");
  TFile *fbkg6 = new TFile("rootFiles/AllJet/OF/WgammaStar.root","read");
  TFile *fbkg7 = new TFile("rootFiles/AllJet/OF/WJetsFakes_Total.root","read");
  TFile *fbkg8 = new TFile("rootFiles/AllJet/OF/DY.root","read");
  TFile *fbkg9 = new TFile("rootFiles/AllJet/OF/DYtautau.root","read");
  TFile *fbkg10 = new TFile("rootFiles/AllJet/OF/Zgamma.root","read");
  TFile *fbkg11 = new TFile("rootFiles/AllJet/OF/VVV.root","read");
  TFile *fbkg12 = new TFile("rootFiles/AllJet/OF/HWW125.root","read");
  TFile *fbkg13 = new TFile("rootFiles/AllJet/OF/VBF.root","read");

  //defining graphs
  TGraph *g = new TGraph();
  TGraph *g2 = new TGraph();
  TGraph *g3 = new TGraph();
  TGraph *g4 = new TGraph();

  //defining signal histograms
  TH1F* hsig1 = (TH1F*)fsig1 -> Get("hHtAfter3");
  TH1F* hsig2 = (TH1F*)fsig2 -> Get("hHtAfter3");

  //defining background histograms
  TH1F *hbkg[13];

  hbkg[0]  = (TH1F*)fbkg13 -> Get("hHtAfter3");
  hbkg[1]  = (TH1F*)fbkg1 -> Get("hHtAfter3");
  hbkg[2]  = (TH1F*)fbkg2 -> Get("hHtAfter3");
  hbkg[3]  = (TH1F*)fbkg3 -> Get("hHtAfter3");
  hbkg[4]  = (TH1F*)fbkg4 -> Get("hHtAfter3");
  hbkg[5]  = (TH1F*)fbkg5 -> Get("hHtAfter3");
  hbkg[6]  = (TH1F*)fbkg6 -> Get("hHtAfter3");
  hbkg[7]  = (TH1F*)fbkg7 -> Get("hHtAfter3");
  hbkg[8]  = (TH1F*)fbkg8 -> Get("hHtAfter3");
  hbkg[9]  = (TH1F*)fbkg9 -> Get("hHtAfter3");
  hbkg[10]  = (TH1F*)fbkg10 -> Get("hHtAfter3");
  hbkg[11]  = (TH1F*)fbkg11 -> Get("hHtAfter3");
  hbkg[12]  = (TH1F*)fbkg12 -> Get("hHtAfter3");
  
  //total number of signal events (no cuts applied)
  Float_t totSig = hsig1->Integral(0.,3000.) + hsig2->Integral(0.,3000.);

  //total number of background events (no cuts applied)
  Float_t totBkg = 0;
  for(int q = 0; q < 13; ++q)
    totBkg = totBkg + hbkg[q] ->Integral(0.,3000.);
  //totBkg = hbkg[1] ->Integral(0.,3000.);

  //defining variables
  Float_t effSig = 0.;
  Float_t effBkg = 0.;

  //main loop: scanning efficiencies
  for (Int_t i = 0; i < 3000; ++i){

    //calculating signal efficiency
    Float_t sum = 0.;
    sum = hsig1->Integral(0.,i) + hsig2->Integral(0.,i);
    effSig = sum / totSig;
    
    //calculating background efficiency (and rejection)
    Float_t bkgInt = 0;
    for(int pp = 0; pp < 13; ++pp)
      bkgInt = bkgInt + hbkg[pp]->Integral(0.,i);
    //bkgInt = hbkg[1]->Integral(0.,i);
    effBkg = bkgInt / totBkg; 

    //filling graphs
    g ->SetPoint(i, effSig, 1. - effBkg);
    Float_t den = (bkgInt + hsig1->Integral(0.,i) + hsig2->Integral(0.,i));
    if (den != 0)
      g2->SetPoint(i, i, (hsig1->Integral(0.,i) + hsig2->Integral(0.,i)) / sqrt(den));
    
    g3->SetPoint(i, i, effSig);
    g4->SetPoint(i, i, 1. - effBkg);
  }

  g->SetTitle("ROC");
  g->GetXaxis()->SetTitle("Signal Efficiency");
  g->GetYaxis()->SetTitle("Background Rejection");
  g->SetLineWidth(5);  
  g->GetXaxis()->SetTitleSize(0.05);
  g->GetXaxis()->SetTitleOffset(1.2);
  g->GetXaxis()->SetLabelSize(0.05);
  g->GetYaxis()->SetTitleOffset(1.4);
  g->GetYaxis()->SetTitleSize(0.05);
  g->GetYaxis()->SetLabelSize(0.05);

  TGraph *p = new TGraph();
  p->SetPoint(0.,0.84,0.55);
  p->SetMarkerColor(kRed);
  p->SetMarkerStyle(20);
  p->SetMarkerSize(2);
  
  g2->SetTitle("Significance");
  g2->GetXaxis()->SetTitle("Ht [GeV]");
  g2->SetLineWidth(5);
  g2->GetXaxis()->SetRangeUser(0.,500.);
  g2->GetYaxis()->SetRangeUser(0.,70.);
  g2->GetXaxis()->SetTitleSize(0.05);
  g2->GetXaxis()->SetTitleOffset(1.2);
  g2->GetXaxis()->SetLabelSize(0.05);
  g2->GetYaxis()->SetTitle("S / #sqrt{S + B}");
  g2->GetYaxis()->SetTitleOffset(1.4);
  g2->GetYaxis()->SetTitleSize(0.05);
  g2->GetYaxis()->SetLabelSize(0.05);

  TLine *l = new TLine(250.,0.,250.,70.);
  l->SetLineColor(kViolet);
  l->SetLineWidth(5);

  g3->SetLineColor(kRed);
  g3->SetLineWidth(5);
  g3->SetTitle("Signal Efficiency and Background Rejection");
  g3->GetXaxis()->SetTitle("Ht [GeV]");
  g3->GetYaxis()->SetTitle("Signal Efficiency");
  g3->GetXaxis()->SetRangeUser(0.,500.);  
  g3->GetYaxis()->SetRangeUser(0.,1.);
  g3->GetXaxis()->SetLabelOffset(0.02);
  g3->GetXaxis()->SetTitleSize(0.05);
  g3->GetXaxis()->SetTitleOffset(1.3);
  g3->GetXaxis()->SetLabelSize(0.05);
  g3->GetYaxis()->SetTitleOffset(1.4);
  g3->GetYaxis()->SetTitleSize(0.05);
  g3->GetYaxis()->SetLabelSize(0.05);
  g4->SetLineColor(kBlack);
  g4->SetLineWidth(5);

  TLegend* leg = new TLegend(0.60,0.50,0.75,0.70);
  leg->AddEntry(g3,"Sig Eff","l");
  leg->AddEntry(g4,"Bkg Rej","l");
  leg->SetTextSize(0.04);
  leg->SetFillColor(kWhite);
  leg->SetLineColor(kWhite);

  TLine *ll = new TLine(250.,0.,250.,1.);
  ll->SetLineColor(kBlue);
  ll->SetLineWidth(5);

  TCanvas *c1 = new TCanvas("ROC","ROC",600.,600.);
  c1->cd();

  TPad* pad1 = new TPad("pad1", "pad1", 0., 0., 1.0, 1.0);
  pad1->SetLeftMargin(0.15);
  pad1->SetBottomMargin(0.15);
  pad1->Draw();
  pad1->cd();

  g->Draw("AL");
  p->Draw("Psame");

  c1->Print("RocAll.pdf","pdf");  

  TCanvas *c2 = new TCanvas("Significance","Significance",600.,600.);
  c2->cd();

  TPad* pad2 = new TPad("pad2", "pad2", 0., 0., 1.0, 1.0);
  pad2->SetLeftMargin(0.15);
  pad2->SetBottomMargin(0.15);
  pad2->Draw();
  pad2->cd();

  g2->Draw("AL");
  l->Draw("same");
  //  g4->Draw("same");

  c2->Print("SignificanceAll.pdf","pdf");  

  TCanvas *c3 = new TCanvas("Efficiencies","Efficiencies",600.,600.);
  c3->cd();

  TPad* pad3 = new TPad("pad3", "pad3", 0., 0., 1.0, 1.0);
  pad3->SetLeftMargin(0.15);
  pad3->SetRightMargin(0.15);
  pad3->SetBottomMargin(0.15);
  pad3->Draw();
  pad3->cd();

  g3->Draw("AL");
  g4->Draw("Lsame");
  ll->Draw("same");
  pad3->Update();
  TGaxis *axis = new TGaxis(pad3->GetUxmax(),pad3->GetUymin(),pad3->GetUxmax(),pad3->GetUymax(),0,1,510,"+L");
  axis->SetLabelSize(g3->GetXaxis()->GetLabelSize());
  axis->SetLabelFont(g3->GetXaxis()->GetLabelFont());
  axis->SetTitleSize(g3->GetXaxis()->GetTitleSize());
  axis->SetTitleFont(g3->GetXaxis()->GetTitleFont());
  axis->SetTitleOffset(1.4);
  axis->SetTitle("Background Rejection");
  axis->Draw("same");
  leg->Draw("same");

  c3->Print("EfficienciesAll.pdf","pdf");
  /*  hsig2->Draw();
  hsig1->Draw("same");
  hbkg[1]->Draw("same");
  */
  
  return 0;
}
