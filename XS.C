#include "TFile.h"
#include "TH1F.h"
#include "TInterpreter.h"

#include "DY.C"
#include "Top.C"


const Double_t nlo8tev       = 59.8;  // [pb]
const Double_t nlo8tevPlus   = 0.022 * nlo8tev / 1e2;
const Double_t nlo8tevMinus  = 0.019 * nlo8tev / 1e2;
const Double_t BR_WW_to_lnln = (3 * 0.108) * (3 * 0.108);

const Double_t ggWW_xs = nlo8tev * BR_WW_to_lnln * 0.03;
const Double_t qqWW_xs = nlo8tev * BR_WW_to_lnln * 0.97;


// Madgraph # events == 1933235
// MCNLO #events == 539595
// Powgeh # events == 999864

const Double_t NTotalqqWW = 999864;
const Double_t NTotalggWW =  109987;

// Fidutial Xsec 
//const Double_t NTotalqqWW = 172114;
//const Double_t NTotalggWW = 21410;

const UInt_t Nsyst = 36+18;
              

// Relative errors directly from data cards (in %) 
//------------------------------------------------------------------------------

enum backgrounds {
   WW,
   ggWW,
   Vg,
   Wjets,
   VVV,
   DYtt,
   VV,
   ggH,
   VgS
};
 

enum myChannels {
  jet0_SF=0*9,
  jet0_OF=1*9,
  jet1_SF=2*9,
  jet1_OF=3*9,
  jet0=4*9,
  jet1=5*9
};



//[*] WW / ggWW/  Vg/ Wjets/ VVV/ DY/ VV/ ggH/ VgS

Double_t Systematics [Nsyst] = { 6.9,  30.4,  31.3,  27.4,  50.6,   2.6,   5.9,  32.5,  34.1,   // jet0_SF
				 6.5,  30.4,  30.0,  28.3,  50.4,  30.7,   6.2,  28.7,  30.3,   // jet0_OF
				 13.7,  32.1,  33.3,  23.9,  50.3,   2.0,   6.4,  41.0,  58.3,  // jet1_SF 
				 13.3,  31.8,  30.1,  26.4,  50.3,  30.4,   5.6,  35.5,  31.2,  // jet1_OF
				 9.5,	43.0,  43.4,  39.4,  71.4,  30.8,   8.6,  43.4,	45.6,   // jet0
				 19.1,	45.2,  44.9,  35.6,  71.1,  30.5,   8.5,  54.2,	66.1 }; // jet1



// [*] Using SystematicsqqWW as 12.5,  32.7,  30.1,  26.4,  50.3,  32.0,   5.8,  31.5,  31.2placeholder


TString analysisLevel;


//------------------------------------------------------------------------------
// XS
//------------------------------------------------------------------------------
void XS(Double_t &xsValue,
	Double_t &xsStat,
	Double_t &xsSyst,
	Double_t &xsLumi,
	Double_t  luminosity,
	Int_t     njet,
	TString   channel,
	TString   cutLevel,
	Bool_t    useNM1,
	TString   directory,
	Bool_t    useDataDriven,
	Int_t     printLevel)
{


  Int_t defineChannel = -1; 

  if ( njet == 0 && channel == "SF")  defineChannel  = jet0_SF;
  if ( njet == 0 && channel == "OF")  defineChannel  = jet0_OF;
  if ( njet == 1 && channel == "SF")  defineChannel  = jet1_SF;
  if ( njet == 1 && channel == "OF")  defineChannel  = jet1_OF;
  if ( njet == 0 && channel == "All")  defineChannel  = jet0;
  if ( njet == 1 && channel == "All")  defineChannel  = jet1;

  analysisLevel = cutLevel;

  if (useNM1) analysisLevel += "_NM1";

  TString dyChannel = "SF";

  if (channel.Contains("All") && printLevel > 0) {
    printf("\n");
    printf("   WWFull xs = %.2f pb;   WWto2L xs = %.2f pb\n", nlo8tev, nlo8tev*BR_WW_to_lnln);
    printf(" qqWWFull xs = %.2f pb; qqWWto2L xs = %.2f pb\n", nlo8tev*0.97, nlo8tev*0.97*BR_WW_to_lnln);
  }

  if (channel.Contains("EE"))   dyChannel = "EE";
  if (channel.Contains("MuMu")) dyChannel = "MuMu";

  Double_t dyScaleFactor       = -999;
  Double_t dyScaleFactorForTop = -999;
  Double_t topScaleFactor      = -999;

  Double_t NDY [3];
  Double_t NTop[3];


  // Needed to always have the SF DY scale factor in the top estimation
  //----------------------------------------------------------------------------
  DY(NDY[0],
     NDY[1],
     NDY[2],
     dyScaleFactorForTop,
     njet,
     "SF",  // "SF" --> dyChannel (individual topScaleFactor)
     directory,
     useDataDriven,
     0,
     false);

    Top(NTop[0],
      NTop[1],
      NTop[2],
      topScaleFactor,
      dyScaleFactorForTop, /// should be 1 in case of OF 0nly??? 
      njet,
      "All",  // "All" --> channel (individual topScaleFactor)
      directory,
      useDataDriven,
      printLevel);

  
    DY(NDY[0],
     NDY[1],
     NDY[2],
     dyScaleFactor,
     njet,
     dyChannel,
     directory,
     useDataDriven,
     printLevel,
     false);
  

  // Input files
  //----------------------------------------------------------------------------
  TString path = Form("%s/%djet/%s/", directory.Data(), njet, channel.Data());

  TFile* inputWW       = new TFile(path + "WW_pow.root"); //"WW_pow.root");WW_pow_nnll.root");
  TFile* inputggWW     = new TFile(path + "ggWWto2L.root");  
  TFile* inputqqWW     = new TFile(path + "WWTo2L2Nu_pow_nnll.root");  
  TFile* inputTT       = new TFile(path + "TTbar.root");  
  TFile* inputTW       = new TFile(path + "TW.root");  
  TFile* inputWj       = new TFile(path + "WJetsFakes_Total.root");
  TFile* inputWZ       = new TFile(path + "WZ.root");  
  TFile* inputZZ       = new TFile(path + "ZZ.root");  
  TFile* inputDY       = new TFile(path + "DY.root");  
  TFile* inputDYtautau = new TFile(path + "DYtautau.root");  
  TFile* inputWg       = new TFile(path + "Wgamma.root");  
  TFile* inputWgS      = new TFile(path + "WgammaStar.root"); 
  TFile* inputWgNoS    = new TFile(path + "WgammaNoStar.root");  
  TFile* inputH125     = new TFile(path + "HWW125.root"); 
  TFile* inputZgamma   = new TFile(path + "Zgamma.root");  
  TFile* inputVVV      = new TFile(path + "VVV.root");  

  TFile* inputData     = new TFile(path + "DataRun2012_Total.root");

  TFile* inputDYOF = new TFile(Form("%s/%djet/OF/DY.root", directory.Data(), njet));


  //----------------------------------------------------------------------------
  //
  // Estimate WW cross-section
  //
  //----------------------------------------------------------------------------
  TH1F* hNWW       = (TH1F*) inputWW      ->Get("hW" + analysisLevel);
  TH1F* hNggWW     = (TH1F*) inputggWW    ->Get("hW" + analysisLevel);
  TH1F* hNqqWW     = (TH1F*) inputqqWW    ->Get("hW" + analysisLevel);
  TH1F* hNTT       = (TH1F*) inputTT      ->Get("hW" + analysisLevel);
  TH1F* hNTW       = (TH1F*) inputTW      ->Get("hW" + analysisLevel);
  TH1F* hNWj       = (TH1F*) inputWj      ->Get("hW" + analysisLevel);
  TH1F* hNWZ       = (TH1F*) inputWZ      ->Get("hW" + analysisLevel);
  TH1F* hNZZ       = (TH1F*) inputZZ      ->Get("hW" + analysisLevel);
  TH1F* hNDY       = (TH1F*) inputDY      ->Get("hW" + analysisLevel);
  TH1F* hNDYtautau = (TH1F*) inputDYtautau->Get("hW" + analysisLevel);
  TH1F* hNWg       = (TH1F*) inputWg      ->Get("hW" + analysisLevel);
  TH1F* hNWgS      = (TH1F*) inputWgS     ->Get("hW" + analysisLevel);
  TH1F* hNWgNoS    = (TH1F*) inputWgNoS   ->Get("hW" + analysisLevel);
  TH1F* hNH125     = (TH1F*) inputH125    ->Get("hW" + analysisLevel);
  TH1F* hNZgamma   = (TH1F*) inputZgamma  ->Get("hW" + analysisLevel);
  TH1F* hNVVV      = (TH1F*) inputVVV     ->Get("hW" + analysisLevel);
  TH1F* hNData     = (TH1F*) inputData    ->Get("hW" + analysisLevel);
  TH1F* hNDYOF     = (TH1F*) inputDYOF    ->Get("hW" + analysisLevel);


  // Yields
  //----------------------------------------------------------------------------
  Double_t NggWW[]     = {hNggWW    ->GetBinContent(2), hNggWW    ->GetBinError(2), 0.00};
  Double_t NqqWW[]     = {hNqqWW    ->GetBinContent(2), hNqqWW    ->GetBinError(2), 0.00};
  Double_t NWW[]       = {hNWW      ->GetBinContent(2), hNWW      ->GetBinError(2), 0.00};
  Double_t NWj[]       = {hNWj      ->GetBinContent(2), hNWj      ->GetBinError(2), 0.36*hNWj->GetBinContent(2)};
  Double_t NWZ[]       = {hNWZ      ->GetBinContent(2), hNWZ      ->GetBinError(2), 0.00};
  Double_t NZZ[]       = {hNZZ      ->GetBinContent(2), hNZZ      ->GetBinError(2), 0.00};
  Double_t NDYtautau[] = {hNDYtautau->GetBinContent(2), hNDYtautau->GetBinError(2), 0.50*hNDYtautau->GetBinContent(2)}; // [*]
  Double_t NWg[]       = {hNWg      ->GetBinContent(2), hNWg      ->GetBinError(2), 0.00};
  Double_t NWgS[]      = {hNWgS     ->GetBinContent(2), hNWgS     ->GetBinError(2), 0.00};
  Double_t NWgNoS[]    = {hNWgNoS   ->GetBinContent(2), hNWgNoS   ->GetBinError(2), 0.00};
  Double_t NH125[]     = {hNH125    ->GetBinContent(2), hNH125    ->GetBinError(2), 0.00};
  Double_t NZgamma[]   = {hNZgamma  ->GetBinContent(2), hNZgamma  ->GetBinError(2), 0.00};
  Double_t NVVV[]      = {hNVVV     ->GetBinContent(2), hNVVV     ->GetBinError(2), 0.00};
  Double_t NData[]     = {hNData    ->GetBinContent(2), hNData    ->GetBinError(2)};


  // [*] Andrea: A normalization uncertainty of 50% has been added on final estimation of dytautau.
   

  // Top yields for the non-inclusive channels
  //----------------------------------------------------------------------------
  if (!channel.Contains("All")) {

    Double_t NTT[] = {hNTT->GetBinContent(2), hNTT->GetBinError(2), 0.0};
    Double_t NTW[] = {hNTW->GetBinContent(2), hNTW->GetBinError(2), 0.0};

    NTop[2] = NTop[2] * topScaleFactor * (NTT[0] + NTW[0]) / NTop[0];  // syst.
    NTop[0] = topScaleFactor * (NTT[0] + NTW[0]);                      // yield
    NTop[1] = topScaleFactor * sqrt(NTT[1]*NTT[1] + NTW[1]*NTW[1]);    // stat.
    }


  // DY yields for the {OF, EMu, MuE} channels
  //----------------------------------------------------------------------------
  if (channel.Contains("OF")  ||
      channel.Contains("EMu") ||
      channel.Contains("MuE")) {

    NDY[2] = NDY[2] * hNDY->GetBinContent(2) / NDY[0];  // syst.
    NDY[0] = hNDY->GetBinContent(2);                    // yield  // ---> use MC because is OF 
    NDY[1] = hNDY->GetBinError(2);                      // stat.
  }


  // DY yield for the All channel
  //----------------------------------------------------------------------------
  if (channel.Contains("All")) {

    NDY[2] = NDY[2] * (NDY[0] + hNDYOF->GetBinContent(2)) / NDY[0];                // syst.
    NDY[0] = NDY[0] + hNDYOF->GetBinContent(2);                                    // yield
    NDY[1] = sqrt(NDY[1]*NDY[1] + hNDYOF->GetBinError(2)*hNDYOF->GetBinError(2));  // stat.
  }


  // Add relative systematic uncertainties
  //----------------------------------------------------------------------------
  

  NggWW[2] = Systematics [ggWW+defineChannel]   ;  // Relative systematic uncertainty
  NqqWW[2] = Systematics [WW+defineChannel]   ;  // Relative systematic uncertainty

  cout << WW+defineChannel << endl;

  NWW[2]  = (NggWW[2]*NggWW[0] * NggWW[2]*NggWW[0]);
  NWW[2] += (NqqWW[2]*NqqWW[0] * NqqWW[2]*NqqWW[0]);
  NWW[2]  = sqrt(NWW[2]) / NWW[0];  // Relative systematic uncertainty

  NWZ    [2] =  Systematics [VV+defineChannel] * NWZ    [0] / 1e2;  // Absolute systematic uncertainty
  NZZ    [2] =  Systematics [VV+defineChannel] * NZZ    [0] / 1e2;  // Absolute systematic uncertainty

  NWgNoS [2] =  Systematics [Vg+defineChannel];  // Relative systematic uncertainty
  NWgS   [2] =  Systematics [VgS+defineChannel]; // Relative systematic uncertainty

  NWg    [2] = (NWgNoS [2]*NWgNoS [0] * NWgNoS [2]*NWgNoS [0]) / 1e4;  
  NWg    [2] += (NWgS [2]*NWgS [0] * NWgS [2]*NWgS [0]) / 1e4;
  NWg    [2] = sqrt(NWg [2]);// Absolute systematic uncertainty

  NH125  [2] =  Systematics [ggH+defineChannel]* NH125  [0] / 1e2;  // Absolute systematic uncertainty
  NZgamma[2] =  Systematics [Vg+defineChannel] * NZgamma[0] / 1e2;  // Absolute systematic uncertainty
  NVVV   [2] =  Systematics [VVV+defineChannel]* NVVV   [0] / 1e2;  // Absolute systematic uncertainty




  Double_t Background = NWj[0] + NWZ[0] + NZZ[0] + NWg[0] + NTop[0] + NDY[0] + NDYtautau[0] + NH125[0]  + NZgamma[0] + NVVV[0];
  
  Double_t statErrorB = sqrt(NWj      [1]*NWj      [1] +
			     NWZ      [1]*NWZ      [1] +
			     NZZ      [1]*NZZ      [1] +
			     NWg      [1]*NWg      [1] +
			     NTop     [1]*NTop     [1] +
			     NDY      [1]*NDY      [1] +
			     NDYtautau[1]*NDYtautau[1] +
			     NH125    [1]*NH125    [1] +
			     NZgamma  [1]*NZgamma  [1] +
			     NVVV     [1]*NVVV     [1]);

  Double_t systErrorB = sqrt(NWj      [2]*NWj      [2] +
			     NWZ      [2]*NWZ      [2] +
			     NZZ      [2]*NZZ      [2] +
			     NWg      [2]*NWg      [2] +
			     NTop     [2]*NTop     [2] +
			     NDY      [2]*NDY      [2] +
			     NDYtautau[2]*NDYtautau[2] +
			     NH125    [2]*NH125    [2] +
			     NZgamma  [2]*NZgamma  [2] +
			     NVVV     [2]*NVVV     [2]);

  Double_t totalErrorB = sqrt(statErrorB*statErrorB + systErrorB*systErrorB);


  //----------------------------------------------------------------------------
  //
  // Estimate WW efficiency
  //
  //----------------------------------------------------------------------------
  TH1F* hgg = (TH1F*) inputggWW->Get("hWeff" + analysisLevel);
  TH1F* hqq = (TH1F*) inputqqWW->Get("hWeff" + analysisLevel);

  Double_t NGenggWW = NTotalggWW;
  Double_t NGenqqWW = NTotalqqWW;

  Double_t f_gg = ggWW_xs / (ggWW_xs + qqWW_xs);
  Double_t f_qq = 1 - f_gg;


  Printf("%5.3f, %5.3f" ,hqq->GetBinContent(2),  NGenqqWW);

  Double_t ggWW_efficiency = ratioValue(hgg->GetBinContent(2), NGenggWW);
  Double_t qqWW_efficiency = ratioValue(hqq->GetBinContent(2), NGenqqWW);
  Double_t   WW_efficiency = f_gg*ggWW_efficiency + f_qq*qqWW_efficiency;

  Double_t ggWW_efficiencyErr = ratioError(hgg->GetBinContent(2), NGenggWW, hgg->GetBinError(2), sqrt(NGenggWW));
  Double_t qqWW_efficiencyErr = ratioError(hqq->GetBinContent(2), NGenqqWW, hqq->GetBinError(2), sqrt(NGenqqWW));

  // Systematic component
  Double_t WW_efficiencyErr = (f_gg*NggWW[2]*NggWW[0] + f_qq*NqqWW[2]*NqqWW[0]) / NWW[0];
  WW_efficiencyErr *= WW_efficiency / 1e2;
  WW_efficiencyErr *= WW_efficiencyErr;

  // Statistical component (negligible)
  WW_efficiencyErr += (f_gg*ggWW_efficiencyErr)*(f_gg*ggWW_efficiencyErr);
  WW_efficiencyErr += (f_qq*qqWW_efficiencyErr)*(f_qq*qqWW_efficiencyErr);

  // Absolute efficiency uncertainty
  WW_efficiencyErr = sqrt(WW_efficiencyErr);


  if (printLevel > 0) {
    printf("\n signal efficiencies\n");
    printf(" -------------------------------------------------\n");
    printf("       ggWW efficiency = (%6.3f +- %5.3f)%s (stat.)\n", 1e2*ggWW_efficiency, 1e2*ggWW_efficiencyErr, "%");
    printf("       qqWW efficiency = (%6.3f +- %5.3f)%s (stat.)\n", 1e2*qqWW_efficiency, 1e2*qqWW_efficiencyErr, "%");
    printf("         WW efficiency = (%6.3f +- %5.3f)%s (total)\n", 1e2*  WW_efficiency, 1e2*  WW_efficiencyErr, "%");
  }


  // Estimate WW cross-section
  //----------------------------------------------------------------------------
  Double_t xs = (NData[0] - Background) / (luminosity * WW_efficiency);

  xs /= BR_WW_to_lnln;

  //cout << "datos!! " << (NData[0] - Background) << endl;

  // Relative errors
  //----------------------------------------------------------------------------
  Double_t errxsStats = 1e2 * sqrt(NData[0]) / (NData[0] - Background);
  Double_t errxsBkg   = 1e2 * totalErrorB    / (NData[0] - Background);
  Double_t errxsEff   = 1e2 * WW_efficiencyErr / WW_efficiency;
  Double_t errxsSyst  = sqrt(errxsBkg*errxsBkg + errxsEff*errxsEff);
  Double_t errxsLumi  = 2.6;


  if (printLevel > 0) {
    printf("\n [PAS]");
    printf(" The total uncertainty on the background estimation is about %.0f%s,",
	   1e2 * totalErrorB / Background, "%");
    printf(" [PAS]\n");
  }


  // Save the result
  //----------------------------------------------------------------------------
  xsValue = xs;
  xsStat  = xs * errxsStats / 1e2;
  xsSyst  = xs * errxsSyst  / 1e2;
  xsLumi  = xs * errxsLumi  / 1e2;


  // Print
  //----------------------------------------------------------------------------
  Double_t statErrorVV = sqrt(NWZ[1]*NWZ[1] + NZZ[1]*NZZ[1]);
  Double_t systErrorVV = sqrt(NWZ[2]*NWZ[2] + NZZ[2]*NZZ[2]);

  Double_t statErrorDYAll = sqrt(NDY[1]*NDY[1] + NDYtautau[1]*NDYtautau[1]);
  Double_t systErrorDYAll = sqrt(NDY[2]*NDY[2] + NDYtautau[2]*NDYtautau[2]);

  Double_t statErrorSPlusB = sqrt(NWW[1]*NWW[1]                           + statErrorB*statErrorB);
  Double_t systErrorSPlusB = sqrt((NWW[2]*NWW[0]/1e2)*(NWW[2]*NWW[0]/1e2) + systErrorB*systErrorB);


  if (printLevel > 0) {
    printf("\n \\hline \n");
    printf("[%s] yields \\\\ \n", channel.Data());
    printf(" \\hline \n");
    printf("       sample              &  yield $\\pm$ stat. $\\pm$ syst.\\\\ \n ");
    printf(" \\hline \n");
    printf("       gg to WW            & %6.1f $\\pm$ %5.1f $\\pm$ %5.1f \\\\ \n", NggWW[0],            NggWW[1],        NggWW[2]*NggWW[0]/1e2);
    printf("       qq to WW            & %6.1f $\\pm$ %5.1f $\\pm$ %5.1f \\\\ \n", NqqWW[0],            NqqWW[1],        NqqWW[2]*NqqWW[0]/1e2);
    printf(" \\hline \n");
    printf("       tt+tW               & %6.1f $\\pm$ %5.1f $\\pm$ %5.1f \\\\ \n", NTop[0],             NTop[1],         NTop[2]);
    printf("       W+jets              & %6.1f $\\pm$ %5.1f $\\pm$ %5.1f \\\\ \n", NWj[0],              NWj[1],          NWj[2]);
    printf("       WZ+ZZ               & %6.1f $\\pm$ %5.1f $\\pm$ %5.1f \\\\ \n", NWZ[0]+NZZ[0],       statErrorVV,     systErrorVV);
    printf("       Z/g*                & %6.1f $\\pm$ %5.1f $\\pm$ %5.1f \\\\ \n",   NDY[0]+NDYtautau[0], statErrorDYAll,  systErrorDYAll);
    //printf(" (%.1f Z/g* + %.1f Z/g* -> tautau)\n",                           NDY[0],              NDYtautau[0]);
    printf("       Wg+Wg*              & %6.1f $\\pm$ %5.1f $\\pm$ %5.1f \\\\ \n", NWg[0],              NWg[1],          NWg[2]);
    printf("       Zgamma              & %6.1f $\\pm$ %5.1f $\\pm$ %5.1f \\\\ \n", NZgamma[0],          NZgamma[1],      NZgamma[2]);
    printf("       VVV              & %6.1f $\\pm$ %5.1f $\\pm$ %5.1f \\\\ \n", NVVV[0],          NVVV[1],      NVVV[2]);
    printf("       H125                & %6.1f $\\pm$ %5.1f $\\pm$ %5.1f \\\\ \n", NH125[0],            NH125[1],        NH125[2]);
    printf(" \\hline \n");
    printf("       total background    & %6.1f $\\pm$ %5.1f $\\pm$ %5.1f \\\\ \n", Background,          statErrorB,      systErrorB);
    printf("       signal              & %6.1f $\\pm$ %5.1f $\\pm$ %5.1f \\\\ \n", NWW[0],              NWW[1],          NWW[2]*NWW[0]/1e2);
    printf(" \\hline \n");
    printf("       signal + background & %6.1f $\\pm$ %5.1f $\\pm$ %5.1f \\\\ \n", NWW[0]+Background,   statErrorSPlusB, systErrorSPlusB);
    printf(" \\hline \n");
    printf("       data                & %4.0f \\\\ \n",                           NData[0]);
    printf(" \\hline \n");
    printf("\n [%s] WW cross-section\n", channel.Data());
    printf(" -------------------------------------------------\n");
  }

  printf("       sigmaWW(%s,%s) = %5.2f +- %5.2f (stat.) +- %5.2f (syst.) +- %5.2f (lumi.) pb\n",
    	 channel.Data(), analysisLevel.Data(), xsValue, xsStat, xsSyst, xsLumi);


  // Difference wrt. the theoretical value
  //----------------------------------------------------------------------------
  if (channel.Contains("All") && printLevel > 0) {

    Double_t deltaXS    = xsValue - nlo8tev;
    Double_t deltaXSErr = xsStat*xsStat + xsSyst*xsSyst + xsLumi*xsLumi + nlo8tevPlus*nlo8tevPlus;

    deltaXSErr = sqrt(deltaXSErr);

    printf("\n       sigma(NLO) = %.2f + %.2f - %.2f pb\n", nlo8tev, nlo8tevPlus, nlo8tevMinus);
    printf("\n       sigmaWW(%s) - sigma(NLO) = %.2f +- %.2f pb", channel.Data(), deltaXS, deltaXSErr);
    printf(" = (%.0f +- %.0f)\%s of the theoretical value\n\n", 1e2*deltaXS/nlo8tev, 1e2*deltaXSErr/nlo8tev, "%");
  }


  // WW/Z cross-section ratio
  //----------------------------------------------------------------------------
  if (channel.Contains("All") && printLevel > 0) {

    Double_t wwzRatio = xsValue / (3.342e04);

    printf("\n [%s] WW/Z cross-section ratio\n", channel.Data());
    printf(" -------------------------------------------------\n");
    printf("       sigmaWW/sigmaZ = (%.2f +- %.2f) x 10^{-3}\n\n",
	   1e3*wwzRatio,
	   1e3*wwzRatio*9.4/1e2);
  }
}
