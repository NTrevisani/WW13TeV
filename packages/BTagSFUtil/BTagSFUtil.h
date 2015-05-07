/*************************************************************

Class Usage:

This class should only be used for upgrading and downgrading 
if a single operating point is used in an analysis. 

Based on the BTagSFUtil from Michael Segala

*************************************************************/

#include <Riostream.h>
#include "TRandom3.h"
#include "TMath.h"

class BTagSFUtil{

 public:
    
  BTagSFUtil(TString BTagAlgorithm, int seed=0 );
  ~BTagSFUtil();
    
  
  bool IsTagged(float JetDiscriminant, int JetFlavor, float JetPt, float JetEta, int SystematicFlag, TString DataPeriod = "AB");

 private:

  float ScaleFactorB(float JetPt, int SystematicFlag);
  float ScaleFactorLight(float JetPt, float JetEta, int SystematicFlag, TString DataPeriod);
  float ScaleFactorJet(int JetFlavor, float JetPt, float JetEta, int SystematicFlag, TString DataPeriod);

  float JetTagEfficiency(int JetFlavor, float JetPt, float JetEta);
  float TagEfficiencyB(float JetPt, float JetEta);
  float TagEfficiencyC(float JetPt, float JetEta);
  float TagEfficiencyLight(float JetPt, float JetEta);
  
  TRandom3* rand_;

  TString TaggerName;
  float TaggerCut;
  
};

