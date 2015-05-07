const UInt_t nchannels = 7;

TString channel     [nchannels] = {"All", "SF", "OF", "MuMu",   "EE", "EMu",  "MuE" };
TString channelLabel[nchannels] = {"all", "SF", "DF", "#mu#mu", "ee", "e#mu", "#mue"};


// NNLO WW cross-section
//------------------------------------------------------------------------------
const Double_t nloXS      = 59.84;         // [pb]
const Double_t nloXSplus  =  0.022*nloXS;  // [pb]
const Double_t nloXSminus =  0.019*nloXS;  // [pb]


TGraphErrors* gStat = NULL;
TGraphErrors* gSyst = NULL;
TGraphErrors* gLumi = NULL;

UInt_t   mychannels;
Double_t absoluteMin;
Double_t absoluteMax;


enum {All, SF, OF, MuMu, EE, EMu, MuE};

enum {SILENT, NOTE, DEBUG};


//------------------------------------------------------------------------------
// runXS
//------------------------------------------------------------------------------
//
// ichannel
//     -1 --> do all channels
//
//------------------------------------------------------------------------------
void runXS(Int_t    njet          = 0,
	   Int_t    ichannel      = All,   //was 'All'
	   Double_t luminosity    = 19468,
	   TString  cutLevel      = "TopTagging",
	   Bool_t   useNM1        = false,
	   Int_t    printLevel    = DEBUG,
	   Bool_t   useDataDriven = true,  // false => DY and Top from MC
	   Bool_t   drawTheXS     = true, // was true!!!
	   Double_t yMin          = -999,
	   Double_t yMax          = -999)
{

  cout << ichannel << endl;
  if (ichannel >= Int_t(nchannels)) return;

  TString directory = "rootFiles";//"rootfiles.2014Selection.4L.NewMC"; // synch with Guillelmo?

  gInterpreter->LoadMacro("DY.C+");
  gInterpreter->LoadMacro("Top.C+");
  gInterpreter->LoadMacro("XS.C+");
 
  mychannels = (ichannel > -1) ? 1 : nchannels;

  gStat = new TGraphErrors(mychannels);
  gSyst = new TGraphErrors(mychannels);
  gLumi = new TGraphErrors(mychannels);


  // Loop
  //----------------------------------------------------------------------------
  absoluteMin =  999;
  absoluteMax = -999;

   Double_t xsValue, xsStat, xsSyst, xsLumi;

  for (UInt_t i=0; i<mychannels; i++) {

    UInt_t thechannel = (ichannel > -1) ? ichannel : i;

    cout<<"numero = "<<i<<endl;
    XS(xsValue,
       xsStat,
       xsSyst,
       xsLumi,
       luminosity,
       njet,
       channel[ichannel],
       cutLevel,
       useNM1,
       directory,
       useDataDriven,
       printLevel);

  cout<<"sono qui!!"<<endl;

    gStat->SetPoint(i, i, xsValue);
    gSyst->SetPoint(i, i, xsValue);
    gLumi->SetPoint(i, i, xsValue);

    Double_t xsError = xsStat;

    gStat->SetPointError(i, 0.5, xsError);

    xsError = sqrt(xsError*xsError + xsSyst*xsSyst);

    gSyst->SetPointError(i, 0.5, xsError);

    xsError = sqrt(xsError*xsError + xsLumi*xsLumi);

    gLumi->SetPointError(i, 0.5, xsError);

    if (absoluteMin > (xsValue - xsError)) absoluteMin = xsValue - xsError; 
    if (absoluteMax < (xsValue + xsError)) absoluteMax = xsValue + xsError; 
  }

  if (absoluteMin > (nloXS - nloXSminus)) absoluteMin = nloXS - nloXSminus;
  if (absoluteMax < (nloXS + nloXSplus))  absoluteMax = nloXS + nloXSplus;


  // Optional drawing
  //----------------------------------------------------------------------------
  if (drawTheXS) DrawXS(ichannel, yMin, yMax);
}


//------------------------------------------------------------------------------
// DrawXS
//------------------------------------------------------------------------------
void DrawXS(Int_t    ichannel,
	    Double_t yMin,
	    Double_t yMax)
{
  gLumi->SetLineColor(kGreen+1);
  gSyst->SetLineColor(kBlue);
  gStat->SetMarkerStyle(kFullCircle);

  gLumi->SetLineWidth(2);
  gSyst->SetLineWidth(2);
  gStat->SetLineWidth(2);

 TCanvas* canvas = new TCanvas("wwxs", "wwxs");
  
  TMultiGraph* mg = new TMultiGraph();

  mg->SetTitle("CMS preliminary");

  mg->Add(gLumi);
  mg->Add(gSyst);
  mg->Add(gStat);

  mg->Draw("ap");

  mg->GetYaxis()->SetTitle("#sigma_{WW} [pb]");
  mg->GetYaxis()->SetTitleOffset(1.6);

  (yMin > 0) ? mg->SetMinimum(yMin) : mg->SetMinimum(absoluteMin -  5);
  (yMax > 0) ? mg->SetMaximum(yMax) : mg->SetMaximum(absoluteMax + 30);

  canvas->Update();


  // x-axis label
  //----------------------------------------------------------------------------
  TAxis* xaxis = mg->GetXaxis();

  xaxis->SetLabelSize(0.07);

  for (UInt_t i=0; i<mychannels; i++) {

    UInt_t thechannel = (ichannel > -1) ? ichannel : i;
    xaxis->SetBinLabel(xaxis->FindBin(i), channelLabel[thechannel].Data());
  }

  xaxis->CenterLabels();
  xaxis->LabelsOption("h");


  // NLO WW cross-section
  //----------------------------------------------------------------------------
  TH1F* nlo = new TH1F("nlo", "nlo", mychannels, canvas->GetUxmin(), canvas->GetUxmax());

  for(UInt_t i=0; i<=mychannels; i++) {
    nlo->SetBinContent(i, nloXS + 0.5 * (nloXSplus - nloXSminus));
    nlo->SetBinError  (i,         0.5 * (nloXSplus + nloXSminus));
  }

  nlo->SetLineColor  (0);
  nlo->SetFillColor  (kRed+1);
  nlo->SetFillStyle  (3354);
  nlo->SetMarkerColor(0);

  TLine* line = new TLine(canvas->GetUxmin(), nloXS, canvas->GetUxmax(), nloXS);

  line->SetLineColor(kRed+1);
  line->SetLineWidth(2);


  // Legend
  //----------------------------------------------------------------------------
  TLegend* legend = new TLegend(0.22, 0.65, 0.49, 0.89);

  legend->SetFillColor(0);
  legend->SetTextFont(42);
  legend->SetTextSize(0.035);

  legend->AddEntry(nlo,   " #sigma_{WW}^{NLO}", "f");
  legend->AddEntry(gStat, " #sigma_{WW} #pm stat. ", "lp");
  legend->AddEntry(gSyst, " #sigma_{WW} #pm (stat. #oplus syst.)", "l");
  legend->AddEntry(gLumi, " #sigma_{WW} #pm (stat. #oplus syst. #oplus lumi.)", "l");


  // Put everything together
  //----------------------------------------------------------------------------
  nlo ->Draw("e2,same");
  line->Draw("same");
  mg  ->Draw("p,same");
  legend->Draw("same");


  // And save it
  //----------------------------------------------------------------------------
  canvas->Update();
  canvas->GetFrame()->DrawClone();

  TString suffix = (ichannel > -1) ? channel[ichannel] : "";
  
  gSystem->mkdir("png", kTRUE);

  canvas->SaveAs("png/wwxs" + suffix + ".png");
}
