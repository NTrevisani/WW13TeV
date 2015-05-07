{
  gROOT->LoadMacro("PUWeight.C+");
  PUWeight* pu = new PUWeight(12103.3, Summer12_53X, "2012");
  //PUWeight* pu = new PUWeight(2143.3, Summer11ITSmear, "2011A");
  //PUWeight* pu = new PUWeight(731.1, Fall11, "2011B");
  TCanvas* c = new TCanvas;
  c->Divide(1,2);
  c->cd(1);
  gPad->SetGridx();
  gPad->SetGridy();
  TH1D* hdata = pu.GetDataHisto();
  hdata->SetTitle("PU Profiles");
  hdata->Draw();
  hdata->GetXaxis()->SetTitle("N_{PU}");

  TH1D* hmc = pu.GetMCHisto();
  hmc->SetLineColor(kRed);
  hmc->SetTitle("PU^{MC}");
  hmc->Draw("SAMES");

  TLegend *legend = new TLegend(0.78,0.50,0.98,0.75);
  legend->AddEntry(hdata, "Data", "l");
  legend->AddEntry(hmc, "MC", "l");
  legend->Draw();



  c->cd(2);
  gPad->SetGridx();
  gPad->SetGridy();
  TH1D* hweights = pu.GetWeightsHisto();
  hweights->SetTitle("PU^{WEIGHTS}");
  hweights->Draw();
  hweights->GetXaxis()->SetTitle("N_{PU}");
  
  TLine l(0,1,hweights->GetXaxis()->GetXmax(),1);
  l.SetLineColor(kRed);
  l.Draw();
}
