//run typing:  root -l -b -q 'macroHisto.C("png","logon","normoff")'

#include "TFile.h"
#include "TH1.h"
#include "TCanvas.h"
#include "TROOT.h"
#include "TTree.h"
#include "TChain.h
#include "TString.h"
#include "TLegend.h"
#include "TObjArray.h"
#include <algorithm>
#include "TCut.h"
#include "TStyle.h"
#include <iostream.h>
#include "THStack.h"
#include "TGraphErrors.h"
#include "TLatex.h"

const int nProcesses = 4;

enum {iWW, iTT, iWJets,iData};

TFile *input[nProcesses];
TH1F  *histo[nProcesses];

TGraphErrors* ratio = new TGraphErrors();

TString process[nProcesses];

process[iWW]     = "WW50";
process[iTT]     = "TTbar50";
process[iWJets]  = "WJets50";
process[iData]   = "WW50";

Color_t color[nProcesses];

color[iWW]     = kAzure - 9;
color[iWJets]  = kGray  + 1;
color[iTT]     = kYellow;
color[iData]   = kBlack;

TGraphErrors *errors = new TGraphErrors();

//drawing instructions
void drawPlots(TString variable,
	       Float_t left, 
	       Float_t right, 
	       Int_t   nrebin,
	       TString units = "",
	       TString format,
	       TString drawLog,
	       TString norm,
	       TString Channel,
	       TString StackMode,
	       TString muonID,
	       TString DataMode = "dataon" 
	       ){

  Float_t rangeY = 0.;

  TString path = "rootFiles/" + Channel + "/" + muonID + "/";

  for (int ip = 0; ip < nProcesses; ++ip){
    TString fileName = path + process[ip] + ".root";
    input[ip] = new TFile(fileName, "read");

    if (input[ip] -> GetListOfKeys() -> Contains(variable) == 0){
      cout<<"I cannot find the histogram named "<<variable<<". I'll skip it."<<endl;
      return;
    }

    histo[ip]  = (TH1F*) input[ip] -> Get(variable);
    histo[ip] -> Rebin(nrebin);
    histo[ip] -> GetXaxis() -> SetRangeUser(left,right);

    //histograms normalization
    if (norm == "normon")
      histo[ip] -> Scale(1./histo[ip] -> Integral());
    
    //Y axis range
    if( histo[ip] -> GetMaximum() > rangeY )
      rangeY = histo[ip] -> GetMaximum();
  }
  
  TCanvas *c1 = new TCanvas("variable","variable",600,800);
  c1->cd();
  
  TPad* pad1 = new TPad("pad1", "pad1", 0., 0., 1.0, 1.0);
  pad1->SetLeftMargin(0.20);
  pad1->SetBottomMargin(0.18);
  pad1->Draw();
  pad1->cd();

  histo[0]->SetTitleSize(4.);

  Float_t rangeMin = 0.0001;

  if (drawLog == "logon"){
    histo[0] -> GetYaxis() -> SetRangeUser(rangeMin, rangeY / rangeMin);
    pad1->SetLogy();
  }  
  else if (drawLog == "logoff"){
    histo[0] -> GetYaxis() -> SetRangeUser(0.,2.*rangeY);
  }

  if(units != "[]")
    histo[0]->GetXaxis()->SetTitle(units);
  histo[0]->GetXaxis()->SetTitleSize(0.07);
  histo[0]->GetXaxis()->SetTitleOffset(0.9);
  histo[0]->GetXaxis()->SetLabelSize(0.045);
  histo[0]->GetYaxis()->SetTitleSize(0.05);
  histo[0]->GetYaxis()->SetLabelSize(0.045);  
  histo[0]->GetXaxis()->SetNdivisions(408);
  histo[0]->GetYaxis()->SetNdivisions(408);
  histo[0]->GetYaxis()->SetTitleOffset(2.0);
  histo[0]->GetYaxis()->SetTitle(Form("entries / %.1f", histo[0]->GetBinWidth(0)));

  TLegend* leg = new TLegend(0.25,0.70,0.75,0.89);
  Float_t maxYaxis = 0.;
  
  for (int ip = 0; ip < nProcesses-1; ++ip){
    histo[ip]->SetLineWidth(3);
    histo[ip]->SetStats(0);
    histo[ip]->SetLineColor(color[ip]);
    if( ip == 0 ) histo[ip] -> Draw();
    else histo[ip] -> Draw("same");
    leg->AddEntry(histo[nProcesses -1 - ip],process[nProcesses -1 - ip],"f");
    maxYaxis += histo[ip] -> GetMaximum();
  }
  leg->AddEntry(histo[nProcesses -1],"Data","lep");
  leg->SetTextSize(0.03);
  leg->SetFillColor(kWhite);
  leg->SetLineColor(kWhite);
  leg->Draw();

  if (StackMode == "stackoff")
    c1->Print(variable + "." + format, format);

  //building the tstack
  for (int ip = 0; ip < nProcesses-1; ++ip){
    histo[ip]    -> SetFillColor(color[ip]);
    histo[ip]    -> SetFillStyle(1001);
  }

  TString title  = variable;
  TString njets  = "Inclusive";

  if (title.Contains("Level0")) njets = "0  Jet";
  if (title.Contains("Level1")) njets = "1  Jet";
  if (title.Contains("Level2")) njets = "2+ Jet";

  /*
  if(title.Contains("stack"))
    title.Remove(title.Length() - 5);
  if(title.Contains("TwoLeptonsLevel"))
    title.Remove(title.Length() - 15);
  if(title.Contains("WWLevel"))
    title.Remove(title.Length() - 8);
  cout<<variable<<","<<title<<endl;
  */
  THStack* hstack = new THStack("","");//title, title);
  
  //use this histogram for stack Y axis range
  TH1F haxis("haxis","haxis",histo[0]->GetNbinsX(),histo[0]->GetXaxis()->GetXmin(),histo[0]->GetXaxis()->GetXmax());
  for(int i = 0; i < nProcesses-1; ++i){
  haxis.Add(histo[i]);
  }
  
  Float_t large = haxis.GetBinWidth(0) / 2.;
  //building error graph
  for(int e = 0; e < haxis.GetNbinsX(); ++e){
    errors -> SetPoint(e, haxis.GetXaxis()->GetBinCenter(e), haxis.GetBinContent(e));
    errors -> SetPointError(e, large, haxis.GetBinError(e));
  }
  float maxYaxisStack = 0.;
  Int_t maxBin = haxis.GetMaximumBin();
  maxYaxisStack = haxis.GetBinContent(maxBin) + errors -> GetErrorY(maxBin) / 2;
  
  errors -> SetMarkerStyle(8);
  errors -> SetFillStyle(3005);	
  errors -> SetFillColor(kBlack);
  
  //Y-axis draw options
  if (drawLog == "logon"){
    hstack -> SetMinimum(rangeMin);
    hstack -> SetMaximum(maxYaxisStack / rangeMin);// * rangeMin);
  }

  else if (drawLog == "logoff"){
    hstack -> SetMinimum(0.);
    hstack -> SetMaximum(maxYaxisStack*1.5);
  }

  //actually building the stack
  for (int w = 0; w < nProcesses-1; ++w){
    hstack -> Add(histo[w]);      
  }

  TCanvas *c2 = new TCanvas("stack","stack",600,800);
  c2->cd();
  
  if (DataMode == "nodata")
    TPad* pad2 = new TPad("pad2", "pad2", 0.0, 0.0, 1.0, 1.0);
  
  else if (DataMode == "dataon") {
    TPad* pad2 = new TPad("pad2", "pad2", 0.0, 0.160, 1.0, 1.0);
    TPad* pad3 = new TPad("pad3", "pad3", 0.0, 0.000, 1.0, 0.3);
  }

  c2->Update();
  c2->cd();
  pad2->SetLeftMargin(0.20);
  pad2->SetBottomMargin(0.18);
  pad2->SetTitle(title);
  pad2->Draw();
  pad2->cd();

  if (drawLog == "logon")
    pad2->SetLogy();
  hstack->Draw("hist");

  hstack -> GetXaxis() -> SetRangeUser(left,right);
  hstack -> GetYaxis() -> SetTitle(Form("entries / %.1f", histo[0]->GetBinWidth(0)));
  hstack -> GetYaxis() -> SetTitleOffset(2.0);
  
  hstack -> GetXaxis() -> SetTitleSize(0.07);
  hstack -> GetXaxis() -> SetTitleOffset(0.9);
  hstack -> GetXaxis() -> SetLabelSize(0.045);
  hstack -> GetXaxis() -> SetLabelOffset(0.045);
  hstack -> GetYaxis() -> SetTitleSize(0.05);
  hstack -> GetYaxis() -> SetLabelSize(0.045);  
  hstack -> GetXaxis() -> SetTitle(units);
  hstack -> GetXaxis() -> SetNdivisions(408);
  hstack -> GetYaxis() -> SetNdivisions(408);

  histo[3]->SetMarkerStyle(kFullCircle);
  histo[3]->SetMarkerColor(kBlack);

  histo[iData] -> Draw("ep");
  hstack -> Draw("hist");
  histo[iData] -> Draw("epsame");
  
  errors -> Draw("same,2");

  if(DataMode == "dataon"){
    for(int p = 0; p < histo[0]->GetNbinsX(); ++p){
      if(haxis.GetBinContent(p) != 0){
	ratio->SetPoint(p, histo[0]->GetBinCenter(p), histo[iData]->GetBinContent(p) / haxis.GetBinContent(p));
	ratio->SetPointError(p, 0.5*histo[0]->GetBinWidth(p), 0.1);//histo[iData]->GetBinContent(p) / haxis.GetBinContent(p));
      }
      else 
	ratio->SetPoint(p, histo[0]->GetBinCenter(p), 0.);
    }
    c2->Update();
    c2->cd();
    //ratio->GetXaxis()->SetRangeUser(left,right);
    ratio->GetXaxis()->SetLimits(left,right);
    ratio->GetYaxis()->SetRangeUser(0.,2.);
    ratio->GetXaxis()->SetTitleSize(0.15);
    ratio->GetXaxis()->SetTitleOffset(1.2);
    ratio->GetXaxis()->SetLabelSize(0.10);
    ratio->GetYaxis()->SetTitleSize(0.14);
    ratio->GetYaxis()->SetLabelSize(0.10);  
    ratio->GetXaxis()->SetTitle(units);
    ratio->GetXaxis()->SetNdivisions(306);
    ratio->GetYaxis()->SetNdivisions(306);
    ratio->SetMarkerStyle(8);
    pad3->SetLeftMargin(0.20);
    pad3->SetBottomMargin(0.45);
    pad3->SetTitle(title);
    pad3->Draw();
    pad3->cd();
    ratio->Draw("AP");
  }
  /*  
      TLegend* leg2 = new TLegend(0.20,0.75,0.70,0.89);
      for (int ip = 0; ip < nProcesses; ++ip)
      leg2->AddEntry(histo[ip],process[ip],"l");
      leg2->SetTextSize(0.03);
      leg2->SetFillColor(kWhite);
      leg2->SetLineColor(kWhite);
  */

  pad2->Update();
  pad2->cd();

  leg -> Draw();
  DrawTLatex(0.88, 0.860, 0.04, "CMS preliminary");
  DrawTLatex(0.88, 0.830, 0.03, "L = 5 fb^{-1}");
  DrawTLatex(0.88, 0.780, 0.05, Channel);
  DrawTLatex(0.88, 0.740, 0.04, njets);
  
  pad2->Update();
  
  if (StackMode == "stackon")
    c2 -> Print(variable + "stack." + format, format);
  
  for(int ip = 0; ip < nProcesses-1; ++ip)
    delete histo[ip];
  delete c1;
  delete c2;
  delete hstack;
  
  if (drawLog == "logoff" && norm == "normoff"){
    if (StackMode == "stackoff")
      gSystem->Exec("mv " + variable + "." + format + " " + "distributions/" + Channel + "/" + muonID + "/" + format + "/");
    else
      gSystem->Exec("mv " + variable + "stack." + format + " " + "distributions/" + Channel + "/" + muonID + "/" + format + "/");
  }
  else if(drawLog == "logon" && norm == "normoff"){
    if (StackMode == "stackoff")
      gSystem->Exec("mv " + variable + "." + format + " " + "distributions/" + Channel + "/" + muonID + "/" + format + "Log/");
    else
      gSystem->Exec("mv " + variable + "stack." + format + " " + "distributions/" + Channel + "/" + muonID + "/" + format + "Log/");
  }
  else if(drawLog == "logoff" && norm == "normon"){
    if (StackMode == "stackoff")
      gSystem->Exec("mv " + variable + "." + format + " " + "distributions/" + Channel + "/" + muonID + "/" + format + "Norm/");
    else
      gSystem->Exec("mv " + variable + "stack." + format + " " + "distributions/" + Channel + "/" + muonID + "/" + format + "Norm/");
  }
  else if(drawLog == "logon" && norm == "normon"){
    if (StackMode == "stackoff")
      gSystem->Exec("mv " + variable + "." + format + " " + "distributions/" + Channel + "/" + muonID + "/" + format + "NormLog/");
    else
      gSystem->Exec("mv " + variable + "stack." + format + " " + "distributions/" + Channel + "/" + muonID + "/" + format + "NormLog/");
  }
  
}

//main function
void macroHisto(TString printMode = "", 
		TString logMode   = "", 
		TString normMode  = "",
		TString channel   = "",
		TString stackMode = "",
		TString MuonID    = ""
		){
  
  if(printMode == "" || logMode == "" || normMode == "" || channel == "" || stackMode == "" || MuonID == ""){
    cout<<"****************************************************************************************************"<<endl;
    cout<<"Please choose a set of options."<<endl;
    cout<<"root -l -b -q 'macroHisto.C(printMode,logMode,normMode,channel)', where:"<<endl;
    cout<<"printMode = 'C', 'png' or 'pdf'"<<endl;
    cout<<"logMode = 'logon' or 'logoff'"<<endl;
    cout<<"normMode = 'normon' or 'normoff'"<<endl;
    cout<<"channel = 'All' or 'OF' or 'SF' or 'MuMu' or 'EE' or 'EMu' or 'MuE'"<<endl;
    cout<<"stackMode = 'stackon' or 'stackoff'"<<endl;
    cout<<"MuonID = 'MediumID' or 'MediumIDTighterIP' or 'TightID' or 'TightIDTighterIP'"<<endl;
    cout<<"I suggest, for example:"<<endl;
    cout<<"root -l -b -q 'macroHisto.C(\"pdf\",\"logoff\",\"normoff\",\"OF\",\"stackon\",\"TightIDTighterIP\")'"<<endl;
    cout<<"****************************************************************************************************"<<endl;
    return;
  }

  if (printMode == "C" || printMode == "png" || printMode == "pdf"){
    gSystem->Exec("mkdir distributions");
    gSystem->Exec("mkdir distributions/" + channel);
    gSystem->Exec("mkdir distributions/" + channel + "/" + MuonID);
    if (logMode == "logoff" && normMode == "normoff")
      gSystem->Exec("mkdir distributions/"  + channel + "/" + MuonID + "/" + printMode);
    else if (logMode == "logoff" && normMode == "normon")
      gSystem->Exec("mkdir distributions/"  + channel + "/" + MuonID + "/"  + printMode + "Norm");
    else if (logMode == "logon" && normMode == "normoff")
      gSystem->Exec("mkdir distributions/"  + channel + "/" + MuonID + "/"  + printMode + "Log");
    else if (logMode == "logon" && normMode == "normon")
      gSystem->Exec("mkdir distributions/"  + channel + "/" + MuonID + "/"  + printMode + "NormLog");
  }
  else{
    cout<<"please print a valid plot format: 'C', 'png' or 'pdf'"<<endl;
    return;
  }

  if(logMode != "logon" && logMode != "logoff"){
    cout<<"I'd like to know if you want the log mode to be 'logon' or 'logoff'"<<endl;
    return;
  }

  if(normMode != "normon" && normMode != "normoff"){
    cout<<"I'd like to know if you want the normalization mode to be 'normon' or 'normoff'"<<endl;
    return;
  }

  if (channel != "All" && channel !=  "OF" && channel !=  "SF" && channel !=  "MuMu" && channel !=  "EE" && channel !=  "EMu" && channel !=  "MuE"){
    cout<<"Please select a valid flavour channel to look at: 'All' or 'OF' or 'SF' or 'MuMu' or 'EE' or 'EMu' or 'MuE'"<<endl;
    return;
  }

  if (stackMode != "stackon" && stackMode != "stackoff"){
    cout<<"Please select a valid print mode: 'stackon' or 'stackoff'"<<endl;
    return;
  }

  if (MuonID != "MediumID" && MuonID != "TightID" && MuonID != "TightIDTighterIP" && MuonID != "MediumIDTighterIP"){
    cout<<"Please select the muon ID plots you want to draw"<<endl;
    return;
}

  Int_t cont = 0;
  TString var;
  Float_t leftBound = 0;
  Float_t rightBound = 0;
  TString units = "";
  Int_t nbin = 10;
  ifstream inFile("input.txt");
  std::string line;

  while (getline(inFile,line)){
    
    std::ofstream outFile("out.tmp",std::ios::out);
    outFile<<line<<endl;
    outFile.close();
    std::ifstream input_;
    input_.open("out.tmp",std::ios::in);
    input_ >> var >> leftBound >> rightBound >> nbin >> units;
    input_.close();
    drawPlots(var, leftBound, rightBound, nbin, units, printMode, logMode, normMode, channel, stackMode, MuonID);
  } 
  
  inFile.close();
  gSystem->Exec("rm out.tmp");
  //  gSystem->Exec("rm outHisto.tmp");
}

//------------------------------------------------------------------------------
// DrawTLatex
//------------------------------------------------------------------------------
void DrawTLatex(Double_t x, Double_t y, Double_t tsize, const char* text)
{
  TLatex* tl = new TLatex(x, y, text);

  tl->SetNDC();
  tl->SetTextAlign(   32);
  tl->SetTextFont (   42);
  tl->SetTextSize (tsize);

  tl->Draw("same");
}