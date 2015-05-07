////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////

/////            DUMMY            ////

////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////





#ifndef MUONFR_NEWJETETDEF_H
#define MUONFR_NEWJETETDEF_H 1

#include "PAFAnalysis.h"

#include <TH1F.h>
#include <TMatrix.h>
#include <TH2F.h>
#include "TCounterUI.h"
#include <TLorentzVector.h>
#include "Riostream.h"  
#include "PUWeight.h"

const Double_t ZMASS = 91.1876;

const UInt_t numberMetCuts = 5;
const UInt_t numberDYMVACuts = 5;

Double_t MetCut[] = {20, 25, 30, 45, 1000};

Double_t DYMVACut_0j[numberDYMVACuts] = {-0.9, -0.86, -0.6, 0.88, 1000};

Double_t DYMVACut_1j[numberDYMVACuts] = {-0.9, -0.86, -0.6, 0.84, 1000};

Bool_t runAtOviedo = false;
Bool_t runAtIfca   = !runAtOviedo;

Float_t SelectedChannel;

class test: public PAFAnalysis{

 public:
   test(TTree *tree=0);
   virtual ~test() {}

 protected:
   virtual void              Initialise();
   virtual void              InsideLoop();
   virtual void              SetDataMembersAtTermination();
   virtual void              Summary(); 

   TH1F* h_n_PV; 


   // Counting histograms                                                                  //----------------------------------------------------------------------------       

   TH1F* hWTrigger;
   TH1F* hWMetCut;
   TH1F* hWLowMinv;
   TH1F* hWZVeto;
   TH1F* hWpMetCut;
   TH1F* hWJetVeto;
   TH1F* hWnJetsBeforeBtag;
   TH1F* hWeffnJetsBeforeBtag;
   TH1F* hWnJets;
   TH1F* hWeffnJets;
   TH1F* hWnBtaggedJets;
   TH1F* hWeffnBtaggedJets;
   TH1F* hWnJetsBveto;
   TH1F* hWeffnJetsBveto;
   TH1F* hWnJetsBvetoAfterHt;
   TH1F* hWeffnJetsBvetoAfterHt;

   TH1F* hWDeltaPhiJet;
   TH1F* hWSoftMuVeto;
   TH1F* hWExtraLepton;
   TH1F* hWPtll;
   TH1F* hWTopTagging;

   TH1F* hWeffTrigger;
   TH1F* hWeffMetCut;
   TH1F* hWeffLowMinv;
   TH1F* hWeffZVeto;
   TH1F* hWeffpMetCut;
   TH1F* hWeffJetVeto;
   TH1F* hWeffDeltaPhiJet;
   TH1F* hWeffSoftMuVeto;
   TH1F* hWeffExtraLepton;
   TH1F* hWeffPtll;
   TH1F* hWeffTopTagging;

   // WW level histograms                                                         
   //---------------------------------------------------------------------------- 
   
   TH1F* hPtLepton1WWLevel[4];
   TH1F* hPtLepton2WWLevel[4];
   TH1F* hPtDiLeptonWWLevel[4];
   TH1F* hMinvWWLevel[4];
   TH1F* hMtWWLevel[4];
   TH1F* hNJets30WWLevel[4];
   TH1F* hpfMetWWLevel[4];
   TH1F* hppfMetWWLevel[4];
   TH1F* hchMetWWLevel[4];
   TH1F* hpchMetWWLevel[4];
   TH1F* hpminMetWWLevel[4];
   TH1F* hDeltaRLeptonsWWLevel[4];
   TH1F* hDeltaPhiLeptonsWWLevel[4];
   TH1F* hDPhiPtllJetWWLevel[4];
   TH1F* hDataEvents[4];
   TH1F* hBackgroundEvents[4];

   TH1F* hPtLepton1WWLevelNoHt[4];
   TH1F* hPtLepton2WWLevelNoHt[4];
   TH1F* hPtDiLeptonWWLevelNoHt[4];
   TH1F* hMinvWWLevelNoHt[4];
   TH1F* hMtWWLevelNoHt[4];
   TH1F* hNJets30WWLevelNoHt[4];
   TH1F* hpfMetWWLevelNoHt[4];
   TH1F* hppfMetWWLevelNoHt[4];
   TH1F* hchMetWWLevelNoHt[4];
   TH1F* hpchMetWWLevelNoHt[4];
   TH1F* hpminMetWWLevelNoHt[4];
   TH1F* hDeltaRLeptonsWWLevelNoHt[4];
   TH1F* hDeltaPhiLeptonsWWLevelNoHt[4];
   TH1F* hDPhiPtllJetWWLevelNoHt[4];
   TH1F* hDataEventsNoHt[4];
   TH1F* hBackgroundEventsNoHt[4];

   TH1F* hHt[4];
   TH1F* hHtAfter[4];

   TH1F* h_WWLevel_TightFailEvents;
   TH1F* h_WWLevel_TightTightEvents;
   TH1F* h_WWLevel_TightLooseEvents;


   // TwoLeptons level histograms                                            
   //---------------------------------------------------------------------------- 
   
   TH1F* hPtLepton1TwoLeptonsLevel;
   TH1F* hPtLepton2TwoLeptonsLevel ;
   TH1F* hPtDiLeptonTwoLeptonsLevel;
   TH1F* hMinvTwoLeptonsLevel;
   TH1F* hMtTwoLeptonsLevel;
   TH1F* hNJets30TwoLeptonsLevel;
   TH1F* hpfMetTwoLeptonsLevel;
   TH1F* hppfMetTwoLeptonsLevel;
   TH1F* hchMetTwoLeptonsLevel;
   TH1F* hpchMetTwoLeptonsLevel;
   TH1F* hpminMetTwoLeptonsLevel;
   TH1F* hDeltaRLeptonsTwoLeptonsLevel;
   TH1F* hDeltaPhiLeptonsTwoLeptonsLevel;
   TH1F* hDPhiPtllJetTwoLeptonsLevel;


   TH1F* h_TwoLeptons_TightFailEvents;
   TH1F* h_TwoLeptons_TightTightEvents;
   TH1F* h_TwoLeptons_TightLooseEvents;

   TH1F* hEff;
   TH1F* hEffHt;

   TH1F* hEff0bin;
   TH1F* hEffHt0bin;


   // Dilepton level differential histograms  
   //---------------------------------------------------------------------------

   //differential in dilepton (eta, pt, ptll, mll, dphill), jets (jetpt1, njet), transverse mass? 2-jet bin?     
   /*
   Double_t pt1bins[8] = {20,40,60,80,100,125,150,200};
   Double_t eta1bins[9] = {-2.4,-1.8,-1.2,-0.6,0,0.6,1.2,1.8,2.4};
   Double_t mllbins[9] = {20,40,60,80,100,125,150,175,200};
   Double_t dphibins[13] = {0,0.25,0.5,0.75,1,1.25,1.5,1.75,2,2.25,2.5,2.75,3};
   Double_t ptllbins[8] = {30,40,50,60,70,85,120,150};
   */

   TH1F* hPtLepton1DilepLevel_Diff;
   TH1F* hEtaLepton1DilepLevel_Diff;
   TH1F* hDileptonDilepLevel_Diff;
   TH1F* hMinvDilepLevel_Diff;
   TH1F* hDeltaPhiDilepLevel_Diff;


   // WW level differential histograms
   //----------------------------------------------------------------------------

   TH1F* hPtLepton1WWLevel_Diff;
   TH1F* hEtaLepton1WWLevel_Diff;
   TH1F* hDileptonWWLevel_Diff;
   TH1F* hMinvWWLevel_Diff;
   TH1F* hDeltaPhiWWLevel_Diff;


   // Data-driven methods: Z+jets 
   //--------------------------------------------------------------------------

   TH1F* hNinZevents     [numberDYMVACuts];
   TH1F* hNoutZevents    [numberDYMVACuts];
   TH1F* hNinLooseZevents[numberDYMVACuts];
   TH1F* hMassInZevents  [numberDYMVACuts];
   TH1F* hMassOutZevents [numberDYMVACuts];

   /*
   for (Int_t nC=0; nC<numberDYMVACuts; nC++) {
     hNinZevents     [nC] = new TH1F(Form("hNinZevents%.1i",      nC+1 ), "",   3, 0,   3);
     hNoutZevents    [nC] = new TH1F(Form("hNoutZevents%.1i",     nC+1 ), "",   3, 0,   3);
     hNinLooseZevents[nC] = new TH1F(Form("hNinLooseZevents%.1i", nC+1 ), "",   3, 0,   3);
     hMassInZevents  [nC] = new TH1F(Form("hMassInZevents%.1i",   nC+1 ), "", 200, 0, 200);
     hMassOutZevents [nC] = new TH1F(Form("hMassOutZevents%.1i",  nC+1 ), "", 200, 0, 200);

   */

   // Data-driven methods: Top                                                                                                                                                      
   //----------------------------------------------------------------------------                                                                                                   
   TH1F* hTopTaggedEvents;
   TH1F* hNTopControlRegion;
   TH1F* hNTopTaggedTopControlRegion;

   TH1F* hbTagDisTopTaggedEvents;
   TH1F* hbTagDisNTopControlRegion;
   TH1F* hbTagDisNTopTaggedTopControlRegion;


 public:
 
   // My Declarations:
   // Define global variables

   PUWeight* fPUWeight;
 

   bool G_Debug_DefineAnalysisVariables;

 
   //* Histograms 
  
// * Input parameters
TString Signal; // Type of Signal
int NEvents; // Total number of events in the sample before skim
double Luminosity; // Total luminosity
double XSection; // Process cross section
bool IsDATA; // True if is Data, False in case MC
int WhichRun; // 1 in case of RunI samples. 2 In case of RunII samples.;
 TString TheSample; //path to the input files
 TString flavorChannel; //selected decay channel 
 int jetChannel; //number of jets in the event
   ClassDef(test,0);
};
#endif

