#include "BTagSFUtil.h"

BTagSFUtil::BTagSFUtil(TString BTagAlgorithm, int seed ) {

  rand_ = new TRandom3(seed);

  TaggerName = BTagAlgorithm;

  TaggerCut = -1.;
  //if (TaggerName=="SSVHEM") TaggerCut = 1.74;
  //if (TaggerName=="SSVHPT") TaggerCut = 2.0;
  //if (TaggerName== "TCHEM") TaggerCut = 3.3;
  //if (TaggerName== "TCHEL") TaggerCut = 1.7;
  if (TaggerName== "TCHPT") TaggerCut = 3.41;
  //if (TaggerName== "TCHPM") TaggerCut = 1.93;
  if (TaggerName==   "JPL") TaggerCut = 0.275;
  if (TaggerName==   "JPM") TaggerCut = 0.545;
  if (TaggerName==   "JPT") TaggerCut = 0.790;
  //if (TaggerName==  "JBPL") TaggerCut = 1.33;
  //if (TaggerName==  "JPBM") TaggerCut = 2.55;
  //if (TaggerName==  "JBPT") TaggerCut = 3.74;
  if (TaggerName==  "CSVL") TaggerCut = 0.244;
  if (TaggerName==  "CSVM") TaggerCut = 0.679;
  if (TaggerName==  "CSVT") TaggerCut = 0.898;

  if (TaggerCut<0.) 
    cout << "BTagSFUtil: " << BTagAlgorithm << " not a supported b-tagging algorithm" << endl;

}

BTagSFUtil::~BTagSFUtil() {

  delete rand_;

}

float BTagSFUtil::ScaleFactorB(float JetPt, int SystematicFlag) {
  
  float SFb = -1.;
 
  if (JetPt<20.) { cout << "SFb is not available for jet pt<20 GeV" << endl; return -1.; }

  float x = JetPt;
  if (JetPt>800.) x = 800.;

  const int nBTagBins = 16;
  float ptmin[nBTagBins] = {20, 30, 40, 50, 60, 70, 80, 100, 120, 160, 210, 260, 320, 400, 500, 600};
 
  int JetPtBin = -1;
  for (int ptbtagbin = 0; ptbtagbin<nBTagBins; ptbtagbin++) 
    if (x>=ptmin[ptbtagbin]) JetPtBin++;

  if (TaggerName=="CSVL") {

    //Tagger: CSVL within 20 < pt < 800 GeV, abs(eta) < 2.4, x = pt
    SFb = 0.981149*((1.+(-0.000713295*x))/(1.+(-0.000703264*x)));
    float SFb_error[nBTagBins] = {
      0.0484285,
      0.0126178,
      0.0120027,
      0.0141137,
      0.0145441,
      0.0131145,
      0.0168479,
      0.0160836,
      0.0126209,
      0.0136017,
      0.019182,
      0.0198805,
      0.0386531,
      0.0392831,
      0.0481008,
      0.0474291 };

    if (SystematicFlag!=0) {
      
      if (JetPt>800.) SFb += SystematicFlag*2.*SFb_error[JetPtBin];
      else SFb += SystematicFlag*SFb_error[JetPtBin];

    }
    
  } else if (TaggerName=="CSVM") {

    //Tagger: CSVM within 20 < pt < 800 GeV, abs(eta) < 2.4, x = pt
      SFb = 0.726981*((1.+(0.253238*x))/(1.+(0.188389*x)));
    float SFb_error[nBTagBins] = {
      0.0554504,
      0.0209663,
      0.0207019,
      0.0230073,
      0.0208719,
      0.0200453,
      0.0264232,
      0.0240102,
      0.0229375,
      0.0184615,
      0.0216242,
      0.0248119,
      0.0465748,
      0.0474666,
      0.0718173,
      0.0717567 };

    if (SystematicFlag!=0) {
      
      if (JetPt>800.) SFb += SystematicFlag*2.*SFb_error[JetPtBin];
      else SFb += SystematicFlag*SFb_error[JetPtBin];

    }

  } else if (TaggerName=="CSVT") {

    //Tagger: CSVT within 20 < pt < 800 GeV, abs(eta) < 2.4, x = pt
    SFb = 0.869965*((1.+(0.0335062*x))/(1.+(0.0304598*x)));
    float SFb_error[nBTagBins] = {
      0.0567059,
      0.0266907,
      0.0263491,
      0.0342831,
      0.0303327,
      0.024608,
      0.0333786,
      0.0317642,
      0.031102,
      0.0295603,
      0.0474663,
      0.0503182,
      0.0580424,
      0.0575776,
      0.0769779,
      0.0898199 };

    if (SystematicFlag!=0) {
      
      if (JetPt>800.) SFb += SystematicFlag*2.*SFb_error[JetPtBin];
      else SFb += SystematicFlag*SFb_error[JetPtBin];

    }

  } else if (TaggerName=="TCHPT") {

    //Tagger: TCHPT within 20 < pt < 800 GeV, abs(eta) < 2.4, x = pt
    SFb = 0.305208*((1.+(0.595166*x))/(1.+(0.186968*x)));
    float SFb_error[nBTagBins] = {
      0.0725549,
      0.0275189,
      0.0279695,
      0.028065,
      0.0270752,
      0.0254934,
      0.0262087,
      0.0230919,
      0.0294829,
      0.0226487,
      0.0272755,
      0.0303747,
      0.051223,
      0.0542895,
      0.0589887,
      0.0584216 };

    if (SystematicFlag!=0) {
      
      if (JetPt>800.) SFb += SystematicFlag*2.*SFb_error[JetPtBin];
      else SFb += SystematicFlag*SFb_error[JetPtBin];

    }
    
  } else if (TaggerName=="JPL") {

    //Tagger: JPL within 20 < pt < 800 GeV, abs(eta) < 2.4, x = pt
    SFb = 0.977721*((1.+(-1.02685e-06*x))/(1.+(-2.56586e-07*x)));
    float SFb_error[nBTagBins] = {
      0.0456879,
      0.0229755,
      0.0229115,
      0.0219184,
      0.0222935,
      0.0189195,
      0.0237255,
      0.0236069,
      0.0159177,
      0.0196792,
      0.0168556,
      0.0168882,
      0.0348084,
      0.0355933,
      0.0476836,
      0.0500367 };

    if (SystematicFlag!=0) {
      
      if (JetPt>800.) SFb += SystematicFlag*2.*SFb_error[JetPtBin];
      else SFb += SystematicFlag*SFb_error[JetPtBin];

    }

  } else if (TaggerName=="JPM") {

    //Tagger: JPM within 20 < pt < 800 GeV, abs(eta) < 2.4, x = pt
    SFb = 0.87887*((1.+(0.0393348*x))/(1.+(0.0354499*x)));
    float SFb_error[nBTagBins] = {
      0.0584144,
      0.0304763,
      0.0311788,
      0.0339226,
      0.0343223,
      0.0303401,
      0.0329372,
      0.0339472,
      0.0368516,
      0.0319189,
      0.0354756,
      0.0347098,
      0.0408868,
      0.0415471,
      0.0567743,
      0.0605397 };

    if (SystematicFlag!=0) {
      
      if (JetPt>800.) SFb += SystematicFlag*2.*SFb_error[JetPtBin];
      else SFb += SystematicFlag*SFb_error[JetPtBin];

    }
    
  } else if (TaggerName=="JPT") {

    //Tagger: JPT within 20 < pt < 800 GeV, abs(eta) < 2.4, x = pt
    SFb = 0.802097*((1.+(0.013219*x))/(1.+(0.0107842*x)));
    float SFb_error[nBTagBins] = {
      0.0673183,
      0.0368276,
      0.037958,
      0.0418136,
      0.0463115,
      0.0409334,
      0.0436405,
      0.0419725,
      0.0451182,
      0.0394386,
      0.0423327,
      0.0393015,
      0.0499883,
      0.0509444,
      0.0780023,
      0.0856582 };

    if (SystematicFlag!=0) {
      
      if (JetPt>800.) SFb += SystematicFlag*2.*SFb_error[JetPtBin];
      else SFb += SystematicFlag*SFb_error[JetPtBin];

    }
  
  }

  return SFb;

}


float BTagSFUtil::ScaleFactorLight(float JetPt, float JetEta, int SystematicFlag, TString DataPeriod) {

  float SFl = -1.;
 
  if (JetPt<20.) { cout << "SFb is not available for jet pt<20 GeV" << endl; return -1.; }

  float x = JetPt;
  if (JetPt>800.) x = 800.;
  if ((TaggerName=="JPL" || TaggerName=="CSVL") && fabs(JetEta)>=1.5 && JetPt>700.) x = 700.;
  if ((TaggerName=="JPM" || TaggerName=="CSVM") && fabs(JetEta)>=1.6 && JetPt>700.) x = 700.;

  if (DataPeriod=="ABCD") {

    // Begin of definition of functions from SF_12ABCD ---------------
    
    if( TaggerName=="CSVL" && fabs(JetEta)>=0.0 && fabs(JetEta)<0.5)
      {
	if( SystematicFlag== 0 ) SFl = ((1.04901+(0.00152181*x))+(-3.43568e-06*(x*x)))+(2.17219e-09*(x*(x*x)));
	if( SystematicFlag==-1 ) SFl = ((0.973773+(0.00103049*x))+(-2.2277e-06*(x*x)))+(1.37208e-09*(x*(x*x)));
	if( SystematicFlag==+1 ) SFl = ((1.12424+(0.00201136*x))+(-4.64021e-06*(x*x)))+(2.97219e-09*(x*(x*x)));
      }
    if( TaggerName=="CSVL" && fabs(JetEta)>=0.5 && fabs(JetEta)<1.0)
      {
	if( SystematicFlag== 0 ) SFl = ((0.991915+(0.00172552*x))+(-3.92652e-06*(x*x)))+(2.56816e-09*(x*(x*x)));
	if( SystematicFlag==-1 ) SFl = ((0.921518+(0.00129098*x))+(-2.86488e-06*(x*x)))+(1.86022e-09*(x*(x*x)));
	if( SystematicFlag==+1 ) SFl = ((1.06231+(0.00215815*x))+(-4.9844e-06*(x*x)))+(3.27623e-09*(x*(x*x)));
      }
    if( TaggerName=="CSVL" && fabs(JetEta)>=1.0 && fabs(JetEta)<1.5)
      {
	if( SystematicFlag== 0 ) SFl = ((0.962127+(0.00192796*x))+(-4.53385e-06*(x*x)))+(3.0605e-09*(x*(x*x)));
	if( SystematicFlag==-1 ) SFl = ((0.895419+(0.00153387*x))+(-3.48409e-06*(x*x)))+(2.30899e-09*(x*(x*x)));
	if( SystematicFlag==+1 ) SFl = ((1.02883+(0.00231985*x))+(-5.57924e-06*(x*x)))+(3.81235e-09*(x*(x*x)));
      }
    if( TaggerName=="CSVL" && fabs(JetEta)>=1.5 && fabs(JetEta)<2.4)
      {
	if( SystematicFlag== 0 ) SFl = ((1.06121+(0.000332747*x))+(-8.81201e-07*(x*x)))+(7.43896e-10*(x*(x*x)));
	if( SystematicFlag==-1 ) SFl = ((0.983607+(0.000196747*x))+(-3.98327e-07*(x*x)))+(2.95764e-10*(x*(x*x)));
	if( SystematicFlag==+1 ) SFl = ((1.1388+(0.000468418*x))+(-1.36341e-06*(x*x)))+(1.19256e-09*(x*(x*x)));
      }
    if( TaggerName=="CSVM" && fabs(JetEta)>=0.0 && fabs(JetEta)<0.8)
      {
	if( SystematicFlag== 0 ) SFl = ((1.06238+(0.00198635*x))+(-4.89082e-06*(x*x)))+(3.29312e-09*(x*(x*x)));
	if( SystematicFlag==-1 ) SFl = ((0.972746+(0.00104424*x))+(-2.36081e-06*(x*x)))+(1.53438e-09*(x*(x*x)));
	if( SystematicFlag==+1 ) SFl = ((1.15201+(0.00292575*x))+(-7.41497e-06*(x*x)))+(5.0512e-09*(x*(x*x)));
      }
    if( TaggerName=="CSVM" && fabs(JetEta)>=0.8 && fabs(JetEta)<1.6)
      {
	if( SystematicFlag== 0 ) SFl = ((1.08048+(0.00110831*x))+(-2.96189e-06*(x*x)))+(2.16266e-09*(x*(x*x)));
	if( SystematicFlag==-1 ) SFl = ((0.9836+(0.000649761*x))+(-1.59773e-06*(x*x)))+(1.14324e-09*(x*(x*x)));
	if( SystematicFlag==+1 ) SFl = ((1.17735+(0.00156533*x))+(-4.32257e-06*(x*x)))+(3.18197e-09*(x*(x*x)));
      }
    if( TaggerName=="CSVM" && fabs(JetEta)>=1.6 && fabs(JetEta)<2.4)
      {
	if( SystematicFlag== 0 ) SFl = ((1.09145+(0.000687171*x))+(-2.45054e-06*(x*x)))+(1.7844e-09*(x*(x*x)));
	if( SystematicFlag==-1 ) SFl = ((1.00616+(0.000358884*x))+(-1.23768e-06*(x*x)))+(6.86678e-10*(x*(x*x)));
	if( SystematicFlag==+1 ) SFl = ((1.17671+(0.0010147*x))+(-3.66269e-06*(x*x)))+(2.88425e-09*(x*(x*x)));
      }
    if( TaggerName=="CSVT" && fabs(JetEta)>=0.0 && fabs(JetEta)<2.4)
      {
	if( SystematicFlag== 0 ) SFl = ((1.01739+(0.00283619*x))+(-7.93013e-06*(x*x)))+(5.97491e-09*(x*(x*x)));
	if( SystematicFlag==-1 ) SFl = ((0.953587+(0.00124872*x))+(-3.97277e-06*(x*x)))+(3.23466e-09*(x*(x*x)));
	if( SystematicFlag==+1 ) SFl = ((1.08119+(0.00441909*x))+(-1.18764e-05*(x*x)))+(8.71372e-09*(x*(x*x)));
      }
    if( TaggerName=="JPL" && fabs(JetEta)>=0.0 && fabs(JetEta)<0.5)
      {
	if( SystematicFlag== 0 ) SFl = ((1.05617+(0.000986016*x))+(-2.05398e-06*(x*x)))+(1.25408e-09*(x*(x*x)));
	if( SystematicFlag==-1 ) SFl = ((0.918762+(0.000749113*x))+(-1.48511e-06*(x*x)))+(8.78559e-10*(x*(x*x)));
	if( SystematicFlag==+1 ) SFl = ((1.19358+(0.00122182*x))+(-2.62078e-06*(x*x)))+(1.62951e-09*(x*(x*x)));
      }
    if( TaggerName=="JPL" && fabs(JetEta)>=0.5 && fabs(JetEta)<1.0)
      {
	if( SystematicFlag== 0 ) SFl = ((1.02884+(0.000471854*x))+(-1.15441e-06*(x*x)))+(7.83716e-10*(x*(x*x)));
	if( SystematicFlag==-1 ) SFl = ((0.893017+(0.000369124*x))+(-8.68577e-07*(x*x)))+(5.79006e-10*(x*(x*x)));
	if( SystematicFlag==+1 ) SFl = ((1.16466+(0.000573985*x))+(-1.43899e-06*(x*x)))+(9.88387e-10*(x*(x*x)));
      }
    if( TaggerName=="JPL" && fabs(JetEta)>=1.0 && fabs(JetEta)<1.5)
      {
	if( SystematicFlag== 0 ) SFl = ((1.02463+(0.000907924*x))+(-2.07133e-06*(x*x)))+(1.37083e-09*(x*(x*x)));
	if( SystematicFlag==-1 ) SFl = ((0.89415+(0.000712877*x))+(-1.57703e-06*(x*x)))+(1.02034e-09*(x*(x*x)));
	if( SystematicFlag==+1 ) SFl = ((1.15511+(0.00110197*x))+(-2.56374e-06*(x*x)))+(1.72152e-09*(x*(x*x)));
      }
    if( TaggerName=="JPL" && fabs(JetEta)>=1.5 && fabs(JetEta)<2.4)
      {
	if( SystematicFlag== 0 ) SFl = ((1.05387+(0.000951237*x))+(-2.35437e-06*(x*x)))+(1.66123e-09*(x*(x*x)));
	if( SystematicFlag==-1 ) SFl = ((0.918611+(0.000781707*x))+(-1.8923e-06*(x*x)))+(1.312e-09*(x*(x*x)));
	if( SystematicFlag==+1 ) SFl = ((1.1891+(0.00112006*x))+(-2.81586e-06*(x*x)))+(2.01249e-09*(x*(x*x)));
      }
    if( TaggerName=="JPM" && fabs(JetEta)>=0.0 && fabs(JetEta)<0.8)
      {
	if( SystematicFlag== 0 ) SFl = ((0.980407+(0.00190765*x))+(-4.49633e-06*(x*x)))+(3.02664e-09*(x*(x*x)));
	if( SystematicFlag==-1 ) SFl = ((0.813164+(0.00127951*x))+(-2.74274e-06*(x*x)))+(1.78799e-09*(x*(x*x)));
	if( SystematicFlag==+1 ) SFl = ((1.14766+(0.00253327*x))+(-6.24447e-06*(x*x)))+(4.26468e-09*(x*(x*x)));
      }
    if( TaggerName=="JPM" && fabs(JetEta)>=0.8 && fabs(JetEta)<1.6)
      {
	if( SystematicFlag== 0 ) SFl = ((1.01783+(0.00183763*x))+(-4.64972e-06*(x*x)))+(3.34342e-09*(x*(x*x)));
	if( SystematicFlag==-1 ) SFl = ((0.860873+(0.00110031*x))+(-2.48023e-06*(x*x)))+(1.73776e-09*(x*(x*x)));
	if( SystematicFlag==+1 ) SFl = ((1.17479+(0.00257252*x))+(-6.81377e-06*(x*x)))+(4.94891e-09*(x*(x*x)));
      }
    if( TaggerName=="JPM" && fabs(JetEta)>=1.6 && fabs(JetEta)<2.4)
      {
	if( SystematicFlag== 0 ) SFl = ((0.866685+(0.00396887*x))+(-1.11342e-05*(x*x)))+(8.84085e-09*(x*(x*x)));
	if( SystematicFlag==-1 ) SFl = ((0.740983+(0.00302736*x))+(-8.12284e-06*(x*x)))+(6.281e-09*(x*(x*x)));
	if( SystematicFlag==+1 ) SFl = ((0.992297+(0.00490671*x))+(-1.41403e-05*(x*x)))+(1.14097e-08*(x*(x*x)));
      }
    if( TaggerName=="JPT" && fabs(JetEta)>=0.0 && fabs(JetEta)<2.4)
      {
	if( SystematicFlag== 0 ) SFl = ((0.89627+(0.00328988*x))+(-8.76392e-06*(x*x)))+(6.4662e-09*(x*(x*x)));
	if( SystematicFlag==-1 ) SFl = ((0.666092+(0.00262465*x))+(-6.5345e-06*(x*x)))+(4.73926e-09*(x*(x*x)));
	if( SystematicFlag==+1 ) SFl = ((1.12648+(0.00394995*x))+(-1.0981e-05*(x*x)))+(8.19134e-09*(x*(x*x)));
      }
    if( TaggerName=="TCHPT" && fabs(JetEta)>=0.0 && fabs(JetEta)<2.4)
      {
	if( SystematicFlag== 0 ) SFl = ((1.1676+(0.00136673*x))+(-3.51053e-06*(x*x)))+(2.4966e-09*(x*(x*x)));
	if( SystematicFlag==-1 ) SFl = ((0.988346+(0.000914722*x))+(-2.37077e-06*(x*x)))+(1.72082e-09*(x*(x*x)));
	if( SystematicFlag==+1 ) SFl = ((1.34691+(0.00181637*x))+(-4.64484e-06*(x*x)))+(3.27122e-09*(x*(x*x)));
      }

    // End of definition of functions from SF_12ABCD ---------------

  } else if (DataPeriod=="AB") { 

    // Begin of definition of functions from SF_12AB ---------------

    if( TaggerName=="CSVL" && fabs(JetEta)>=0.0 && fabs(JetEta)<0.5)
      {
	if( SystematicFlag== 0 ) SFl = ((1.00989+(0.00155686*x))+(-3.72647e-06*(x*x)))+(2.47025e-09*(x*(x*x)));
	if( SystematicFlag==-1 ) SFl = ((0.947488+(0.00105091*x))+(-2.43972e-06*(x*x)))+(1.58902e-09*(x*(x*x)));
	if( SystematicFlag==+1 ) SFl = ((1.07229+(0.00206098*x))+(-5.00971e-06*(x*x)))+(3.35179e-09*(x*(x*x)));
      }
    if( TaggerName=="CSVL" && fabs(JetEta)>=0.5 && fabs(JetEta)<1.0)
      {
	if( SystematicFlag== 0 ) SFl = ((0.958598+(0.00173458*x))+(-4.12744e-06*(x*x)))+(2.83257e-09*(x*(x*x)));
	if( SystematicFlag==-1 ) SFl = ((0.900024+(0.00129392*x))+(-3.01708e-06*(x*x)))+(2.06723e-09*(x*(x*x)));
	if( SystematicFlag==+1 ) SFl = ((1.01716+(0.00217335*x))+(-5.23419e-06*(x*x)))+(3.5986e-09*(x*(x*x)));
      }
    if( TaggerName=="CSVL" && fabs(JetEta)>=1.0 && fabs(JetEta)<1.5)
      {
	if( SystematicFlag== 0 ) SFl = ((0.963113+(0.00163674*x))+(-3.84776e-06*(x*x)))+(2.56918e-09*(x*(x*x)));
	if( SystematicFlag==-1 ) SFl = ((0.90596+(0.00125465*x))+(-2.78863e-06*(x*x)))+(1.78602e-09*(x*(x*x)));
	if( SystematicFlag==+1 ) SFl = ((1.02025+(0.0020171*x))+(-4.90389e-06*(x*x)))+(3.35329e-09*(x*(x*x)));
      }
    if( TaggerName=="CSVL" && fabs(JetEta)>=1.5 && fabs(JetEta)<2.4)
      {
	if( SystematicFlag== 0 ) SFl = ((1.04996+(0.00031979*x))+(-8.43322e-07*(x*x)))+(6.9451e-10*(x*(x*x)));
	if( SystematicFlag==-1 ) SFl = ((0.983472+(0.000169396*x))+(-2.82848e-07*(x*x)))+(1.52744e-10*(x*(x*x)));
	if( SystematicFlag==+1 ) SFl = ((1.11645+(0.000469873*x))+(-1.40321e-06*(x*x)))+(1.23681e-09*(x*(x*x)));
      }
    if( TaggerName=="CSVM" && fabs(JetEta)>=0.0 && fabs(JetEta)<0.8)
      {
	if( SystematicFlag== 0 ) SFl = ((1.02213+(0.00189078*x))+(-4.59419e-06*(x*x)))+(3.0366e-09*(x*(x*x)));
	if( SystematicFlag==-1 ) SFl = ((0.946+(0.000940317*x))+(-1.99048e-06*(x*x)))+(1.18343e-09*(x*(x*x)));
	if( SystematicFlag==+1 ) SFl = ((1.09827+(0.00283897*x))+(-7.19354e-06*(x*x)))+(4.89013e-09*(x*(x*x)));
      }
    if( TaggerName=="CSVM" && fabs(JetEta)>=0.8 && fabs(JetEta)<1.6)
      {
	if( SystematicFlag== 0 ) SFl = ((1.0596+(0.00102926*x))+(-2.70312e-06*(x*x)))+(1.82871e-09*(x*(x*x)));
	if( SystematicFlag==-1 ) SFl = ((0.974966+(0.000545735*x))+(-1.23123e-06*(x*x)))+(7.05661e-10*(x*(x*x)));
	if( SystematicFlag==+1 ) SFl = ((1.14423+(0.00151156*x))+(-4.17277e-06*(x*x)))+(2.95233e-09*(x*(x*x)));
      }
    if( TaggerName=="CSVM" && fabs(JetEta)>=1.6 && fabs(JetEta)<2.4)
      {
	if( SystematicFlag== 0 ) SFl = ((1.04976+(0.000897158*x))+(-3.22829e-06*(x*x)))+(2.71316e-09*(x*(x*x)));
	if( SystematicFlag==-1 ) SFl = ((0.977166+(0.000550586*x))+(-1.91114e-06*(x*x)))+(1.44817e-09*(x*(x*x)));
	if( SystematicFlag==+1 ) SFl = ((1.12232+(0.00124269*x))+(-4.54368e-06*(x*x)))+(3.98079e-09*(x*(x*x)));
      }
    if( TaggerName=="CSVT" && fabs(JetEta)>=0.0 && fabs(JetEta)<2.4)
      {
	if( SystematicFlag== 0 ) SFl = ((0.985589+(0.00302526*x))+(-8.73861e-06*(x*x)))+(6.65051e-09*(x*(x*x)));
	if( SystematicFlag==-1 ) SFl = ((0.93612+(0.00131596*x))+(-4.30052e-06*(x*x)))+(3.45957e-09*(x*(x*x)));
	if( SystematicFlag==+1 ) SFl = ((1.03505+(0.00472994*x))+(-1.31661e-05*(x*x)))+(9.84151e-09*(x*(x*x)));
      }
    if( TaggerName=="JPL" && fabs(JetEta)>=0.0 && fabs(JetEta)<0.5)
      {
	if( SystematicFlag== 0 ) SFl = ((1.01004+(0.000693171*x))+(-1.71673e-06*(x*x)))+(1.13601e-09*(x*(x*x)));
	if( SystematicFlag==-1 ) SFl = ((0.920801+(0.00048556*x))+(-1.14573e-06*(x*x)))+(7.29722e-10*(x*(x*x)));
	if( SystematicFlag==+1 ) SFl = ((1.09929+(0.000899912*x))+(-2.28605e-06*(x*x)))+(1.54241e-09*(x*(x*x)));
      }
    if( TaggerName=="JPL" && fabs(JetEta)>=0.5 && fabs(JetEta)<1.0)
      {
	if( SystematicFlag== 0 ) SFl = ((0.983323+(0.00021632*x))+(-8.21701e-07*(x*x)))+(6.67398e-10*(x*(x*x)));
	if( SystematicFlag==-1 ) SFl = ((0.893817+(0.000139244*x))+(-5.53288e-07*(x*x)))+(4.54312e-10*(x*(x*x)));
	if( SystematicFlag==+1 ) SFl = ((1.07283+(0.000292983*x))+(-1.08908e-06*(x*x)))+(8.80497e-10*(x*(x*x)));
      }
    if( TaggerName=="JPL" && fabs(JetEta)>=1.0 && fabs(JetEta)<1.5)
      {
	if( SystematicFlag== 0 ) SFl = ((1.00787+(0.00067391*x))+(-1.85829e-06*(x*x)))+(1.42239e-09*(x*(x*x)));
	if( SystematicFlag==-1 ) SFl = ((0.922717+(0.000501621*x))+(-1.3493e-06*(x*x)))+(1.02068e-09*(x*(x*x)));
	if( SystematicFlag==+1 ) SFl = ((1.09302+(0.000845356*x))+(-2.36546e-06*(x*x)))+(1.82448e-09*(x*(x*x)));
      }
    if( TaggerName=="JPL" && fabs(JetEta)>=1.5 && fabs(JetEta)<2.4)
      {
	if( SystematicFlag== 0 ) SFl = ((1.03546+(0.000774019*x))+(-2.15928e-06*(x*x)))+(1.6934e-09*(x*(x*x)));
	if( SystematicFlag==-1 ) SFl = ((0.946412+(0.000642931*x))+(-1.74696e-06*(x*x)))+(1.34402e-09*(x*(x*x)));
	if( SystematicFlag==+1 ) SFl = ((1.12449+(0.000904468*x))+(-2.57084e-06*(x*x)))+(2.04473e-09*(x*(x*x)));
      }
    if( TaggerName=="JPM" && fabs(JetEta)>=0.0 && fabs(JetEta)<0.8)
      {
	if( SystematicFlag== 0 ) SFl = ((0.95344+(0.00171952*x))+(-4.71763e-06*(x*x)))+(3.41607e-09*(x*(x*x)));
	if( SystematicFlag==-1 ) SFl = ((0.821103+(0.0011014*x))+(-2.81576e-06*(x*x)))+(2.00088e-09*(x*(x*x)));
	if( SystematicFlag==+1 ) SFl = ((1.08578+(0.00233518*x))+(-6.61409e-06*(x*x)))+(4.83128e-09*(x*(x*x)));
      }
    if( TaggerName=="JPM" && fabs(JetEta)>=0.8 && fabs(JetEta)<1.6)
      {
	if( SystematicFlag== 0 ) SFl = ((0.994557+(0.00176506*x))+(-4.95785e-06*(x*x)))+(3.63594e-09*(x*(x*x)));
	if( SystematicFlag==-1 ) SFl = ((0.874792+(0.000997253*x))+(-2.51511e-06*(x*x)))+(1.75184e-09*(x*(x*x)));
	if( SystematicFlag==+1 ) SFl = ((1.11431+(0.0025305*x))+(-7.39562e-06*(x*x)))+(5.52077e-09*(x*(x*x)));
      }
    if( TaggerName=="JPM" && fabs(JetEta)>=1.6 && fabs(JetEta)<2.4)
      {
	if( SystematicFlag== 0 ) SFl = ((0.850708+(0.00373619*x))+(-1.10196e-05*(x*x)))+(9.0243e-09*(x*(x*x)));
	if( SystematicFlag==-1 ) SFl = ((0.756769+(0.0028841*x))+(-8.02579e-06*(x*x)))+(6.29964e-09*(x*(x*x)));
	if( SystematicFlag==+1 ) SFl = ((0.944553+(0.00458472*x))+(-1.40078e-05*(x*x)))+(1.1758e-08*(x*(x*x)));
      }
    if( TaggerName=="JPT" && fabs(JetEta)>=0.0 && fabs(JetEta)<2.4)
      {
	if( SystematicFlag== 0 ) SFl = ((0.869428+(0.00338068*x))+(-9.51813e-06*(x*x)))+(7.08382e-09*(x*(x*x)));
	if( SystematicFlag==-1 ) SFl = ((0.657456+(0.00279437*x))+(-7.29415e-06*(x*x)))+(5.28578e-09*(x*(x*x)));
	if( SystematicFlag==+1 ) SFl = ((1.08141+(0.00396192*x))+(-1.17308e-05*(x*x)))+(8.88194e-09*(x*(x*x)));
      }
    if( TaggerName=="TCHPT" && fabs(JetEta)>=0.0 && fabs(JetEta)<2.4)
      {
	if( SystematicFlag== 0 ) SFl = ((1.10796+(0.00168368*x))+(-4.50964e-06*(x*x)))+(3.21561e-09*(x*(x*x)));
	if( SystematicFlag==-1 ) SFl = ((0.948977+(0.0011449*x))+(-3.05912e-06*(x*x)))+(2.17813e-09*(x*(x*x)));
	if( SystematicFlag==+1 ) SFl = ((1.26699+(0.00221974*x))+(-5.95423e-06*(x*x)))+(4.25254e-09*(x*(x*x)));
      }
    
    // End of definition of functions from SF_12AB ---------------
    
  } else if (DataPeriod=="C") { 

    // Begin of definition of functions from SF_12C ---------------
    
    if( TaggerName=="CSVL" && fabs(JetEta)>=0.0 && fabs(JetEta)<0.5)
      {
	if( SystematicFlag== 0 ) SFl = ((1.03512+(0.00172098*x))+(-4.10286e-06*(x*x)))+(2.72413e-09*(x*(x*x)));
	if( SystematicFlag==-1 ) SFl = ((0.971321+(0.00117532*x))+(-2.71334e-06*(x*x)))+(1.77294e-09*(x*(x*x)));
	if( SystematicFlag==+1 ) SFl = ((1.0989+(0.0022646*x))+(-5.48834e-06*(x*x)))+(3.67551e-09*(x*(x*x)));
      }
    if( TaggerName=="CSVL" && fabs(JetEta)>=0.5 && fabs(JetEta)<1.0)
      {
	if( SystematicFlag== 0 ) SFl = ((0.977454+(0.00186222*x))+(-4.30874e-06*(x*x)))+(2.82227e-09*(x*(x*x)));
	if( SystematicFlag==-1 ) SFl = ((0.917942+(0.00139264*x))+(-3.13422e-06*(x*x)))+(2.02475e-09*(x*(x*x)));
	if( SystematicFlag==+1 ) SFl = ((1.03695+(0.00232982*x))+(-5.47968e-06*(x*x)))+(3.62048e-09*(x*(x*x)));
      }
    if( TaggerName=="CSVL" && fabs(JetEta)>=1.0 && fabs(JetEta)<1.5)
      {
	if( SystematicFlag== 0 ) SFl = ((0.940154+(0.00214045*x))+(-5.30206e-06*(x*x)))+(3.75872e-09*(x*(x*x)));
	if( SystematicFlag==-1 ) SFl = ((0.885078+(0.00170468*x))+(-4.08896e-06*(x*x)))+(2.85628e-09*(x*(x*x)));
	if( SystematicFlag==+1 ) SFl = ((0.995215+(0.00257376*x))+(-6.5103e-06*(x*x)))+(4.66211e-09*(x*(x*x)));
      }
    if( TaggerName=="CSVL" && fabs(JetEta)>=1.5 && fabs(JetEta)<2.4)
      {
	if( SystematicFlag== 0 ) SFl = ((1.04882+(0.000373418*x))+(-1.00316e-06*(x*x)))+(8.52325e-10*(x*(x*x)));
	if( SystematicFlag==-1 ) SFl = ((0.982642+(0.000211816*x))+(-4.11471e-07*(x*x)))+(2.88443e-10*(x*(x*x)));
	if( SystematicFlag==+1 ) SFl = ((1.11499+(0.000534645*x))+(-1.59409e-06*(x*x)))+(1.41682e-09*(x*(x*x)));
      }
    if( TaggerName=="CSVM" && fabs(JetEta)>=0.0 && fabs(JetEta)<0.8)
      {
	if( SystematicFlag== 0 ) SFl = ((1.0444+(0.00216756*x))+(-5.4224e-06*(x*x)))+(3.69351e-09*(x*(x*x)));
	if( SystematicFlag==-1 ) SFl = ((0.966203+(0.00112979*x))+(-2.56147e-06*(x*x)))+(1.65183e-09*(x*(x*x)));
	if( SystematicFlag==+1 ) SFl = ((1.1226+(0.00320252*x))+(-8.27754e-06*(x*x)))+(5.73519e-09*(x*(x*x)));
      }
    if( TaggerName=="CSVM" && fabs(JetEta)>=0.8 && fabs(JetEta)<1.6)
      {
	if( SystematicFlag== 0 ) SFl = ((1.05203+(0.00138588*x))+(-3.97677e-06*(x*x)))+(3.13655e-09*(x*(x*x)));
	if( SystematicFlag==-1 ) SFl = ((0.966774+(0.000855535*x))+(-2.33883e-06*(x*x)))+(1.86063e-09*(x*(x*x)));
	if( SystematicFlag==+1 ) SFl = ((1.13729+(0.00191432*x))+(-5.61018e-06*(x*x)))+(4.41282e-09*(x*(x*x)));
      }
    if( TaggerName=="CSVM" && fabs(JetEta)>=1.6 && fabs(JetEta)<2.4)
      {
	if( SystematicFlag== 0 ) SFl = ((1.06547+(0.000850114*x))+(-2.76694e-06*(x*x)))+(1.75015e-09*(x*(x*x)));
	if( SystematicFlag==-1 ) SFl = ((0.992673+(0.000455214*x))+(-1.29572e-06*(x*x)))+(3.89704e-10*(x*(x*x)));
	if( SystematicFlag==+1 ) SFl = ((1.13823+(0.00124422*x))+(-4.23813e-06*(x*x)))+(3.11339e-09*(x*(x*x)));
      }
    if( TaggerName=="CSVT" && fabs(JetEta)>=0.0 && fabs(JetEta)<2.4)
      {
	if( SystematicFlag== 0 ) SFl = ((1.00197+(0.00266395*x))+(-6.95018e-06*(x*x)))+(4.91042e-09*(x*(x*x)));
	if( SystematicFlag==-1 ) SFl = ((0.948887+(0.00103466*x))+(-2.88118e-06*(x*x)))+(2.07782e-09*(x*(x*x)));
	if( SystematicFlag==+1 ) SFl = ((1.05505+(0.00428961*x))+(-1.10115e-05*(x*x)))+(7.74319e-09*(x*(x*x)));
      }
    if( TaggerName=="JPL" && fabs(JetEta)>=0.0 && fabs(JetEta)<0.5)
      {
	if( SystematicFlag== 0 ) SFl = ((1.04676+(0.00112324*x))+(-2.52493e-06*(x*x)))+(1.65931e-09*(x*(x*x)));
	if( SystematicFlag==-1 ) SFl = ((0.954301+(0.000859775*x))+(-1.83391e-06*(x*x)))+(1.17383e-09*(x*(x*x)));
	if( SystematicFlag==+1 ) SFl = ((1.13922+(0.00138539*x))+(-3.21336e-06*(x*x)))+(2.14483e-09*(x*(x*x)));
      }
    if( TaggerName=="JPL" && fabs(JetEta)>=0.5 && fabs(JetEta)<1.0)
      {
	if( SystematicFlag== 0 ) SFl = ((1.01822+(0.000554752*x))+(-1.44083e-06*(x*x)))+(9.9442e-10*(x*(x*x)));
	if( SystematicFlag==-1 ) SFl = ((0.925322+(0.00044089*x))+(-1.09668e-06*(x*x)))+(7.37906e-10*(x*(x*x)));
	if( SystematicFlag==+1 ) SFl = ((1.11112+(0.000667913*x))+(-1.78357e-06*(x*x)))+(1.25108e-09*(x*(x*x)));
      }
    if( TaggerName=="JPL" && fabs(JetEta)>=1.0 && fabs(JetEta)<1.5)
      {
	if( SystematicFlag== 0 ) SFl = ((1.00672+(0.000952424*x))+(-2.24525e-06*(x*x)))+(1.49885e-09*(x*(x*x)));
	if( SystematicFlag==-1 ) SFl = ((0.922111+(0.000746635*x))+(-1.68093e-06*(x*x)))+(1.07795e-09*(x*(x*x)));
	if( SystematicFlag==+1 ) SFl = ((1.09131+(0.00115721*x))+(-2.80782e-06*(x*x)))+(1.92031e-09*(x*(x*x)));
      }
    if( TaggerName=="JPL" && fabs(JetEta)>=1.5 && fabs(JetEta)<2.4)
      {
	if( SystematicFlag== 0 ) SFl = ((1.0346+(0.00102282*x))+(-2.61072e-06*(x*x)))+(1.91999e-09*(x*(x*x)));
	if( SystematicFlag==-1 ) SFl = ((0.946009+(0.00085601*x))+(-2.11306e-06*(x*x)))+(1.5133e-09*(x*(x*x)));
	if( SystematicFlag==+1 ) SFl = ((1.12316+(0.00118887*x))+(-3.10773e-06*(x*x)))+(2.32911e-09*(x*(x*x)));
      }
    if( TaggerName=="JPM" && fabs(JetEta)>=0.0 && fabs(JetEta)<0.8)
      {
	if( SystematicFlag== 0 ) SFl = ((0.976097+(0.00191673*x))+(-4.58557e-06*(x*x)))+(3.10331e-09*(x*(x*x)));
	if( SystematicFlag==-1 ) SFl = ((0.838655+(0.00127127*x))+(-2.69532e-06*(x*x)))+(1.73384e-09*(x*(x*x)));
	if( SystematicFlag==+1 ) SFl = ((1.11354+(0.00255978*x))+(-6.47089e-06*(x*x)))+(4.47278e-09*(x*(x*x)));
      }
    if( TaggerName=="JPM" && fabs(JetEta)>=0.8 && fabs(JetEta)<1.6)
      {
	if( SystematicFlag== 0 ) SFl = ((1.01871+(0.00173675*x))+(-4.58934e-06*(x*x)))+(3.38512e-09*(x*(x*x)));
	if( SystematicFlag==-1 ) SFl = ((0.896752+(0.000935708*x))+(-2.09673e-06*(x*x)))+(1.47779e-09*(x*(x*x)));
	if( SystematicFlag==+1 ) SFl = ((1.14066+(0.00253557*x))+(-7.07715e-06*(x*x)))+(5.29295e-09*(x*(x*x)));
      }
    if( TaggerName=="JPM" && fabs(JetEta)>=1.6 && fabs(JetEta)<2.4)
      {
	if( SystematicFlag== 0 ) SFl = ((0.857254+(0.00402389*x))+(-1.15649e-05*(x*x)))+(9.37845e-09*(x*(x*x)));
	if( SystematicFlag==-1 ) SFl = ((0.766329+(0.00306103*x))+(-8.27326e-06*(x*x)))+(6.43552e-09*(x*(x*x)));
	if( SystematicFlag==+1 ) SFl = ((0.948081+(0.00498302*x))+(-1.48512e-05*(x*x)))+(1.23315e-08*(x*(x*x)));
      }
    if( TaggerName=="JPT" && fabs(JetEta)>=0.0 && fabs(JetEta)<2.4)
      {
	if( SystematicFlag== 0 ) SFl = ((0.911203+(0.00300921*x))+(-8.03854e-06*(x*x)))+(5.97264e-09*(x*(x*x)));
	if( SystematicFlag==-1 ) SFl = ((0.698853+(0.00240605*x))+(-5.86773e-06*(x*x)))+(4.23559e-09*(x*(x*x)));
	if( SystematicFlag==+1 ) SFl = ((1.12356+(0.00360804*x))+(-1.01996e-05*(x*x)))+(7.70949e-09*(x*(x*x)));
      }
    if( TaggerName=="TCHPT" && fabs(JetEta)>=0.0 && fabs(JetEta)<2.4)
      {
	if( SystematicFlag== 0 ) SFl = ((1.1393+(0.00148115*x))+(-3.72335e-06*(x*x)))+(2.6087e-09*(x*(x*x)));
	if( SystematicFlag==-1 ) SFl = ((0.977185+(0.000957694*x))+(-2.36635e-06*(x*x)))+(1.65373e-09*(x*(x*x)));
	if( SystematicFlag==+1 ) SFl = ((1.30146+(0.00200233*x))+(-5.07544e-06*(x*x)))+(3.56314e-09*(x*(x*x)));
      }
    
    // End of definition of functions from SF_12C ---------------
    
  } else if (DataPeriod=="D") { 
    
    // Begin of definition of functions from SF_12D ---------------
    
    if( TaggerName=="CSVL" && fabs(JetEta)>=0.0 && fabs(JetEta)<0.5)
      {
	if( SystematicFlag== 0 ) SFl = ((1.1121+(0.00156291*x))+(-3.72267e-06*(x*x)))+(2.54276e-09*(x*(x*x)));
	if( SystematicFlag==-1 ) SFl = ((1.04345+(0.00100049*x))+(-2.27285e-06*(x*x)))+(1.53238e-09*(x*(x*x)));
	if( SystematicFlag==+1 ) SFl = ((1.18074+(0.00212352*x))+(-5.16888e-06*(x*x)))+(3.55347e-09*(x*(x*x)));
      }
    if( TaggerName=="CSVL" && fabs(JetEta)>=0.5 && fabs(JetEta)<1.0)
      {
	if( SystematicFlag== 0 ) SFl = ((1.05107+(0.0018085*x))+(-4.42378e-06*(x*x)))+(3.12722e-09*(x*(x*x)));
	if( SystematicFlag==-1 ) SFl = ((0.986937+(0.00132072*x))+(-3.17261e-06*(x*x)))+(2.25152e-09*(x*(x*x)));
	if( SystematicFlag==+1 ) SFl = ((1.11519+(0.00229425*x))+(-5.67093e-06*(x*x)))+(4.00366e-09*(x*(x*x)));
      }
    if( TaggerName=="CSVL" && fabs(JetEta)>=1.0 && fabs(JetEta)<1.5)
      {
	if( SystematicFlag== 0 ) SFl = ((0.984747+(0.00233796*x))+(-5.84283e-06*(x*x)))+(4.21798e-09*(x*(x*x)));
	if( SystematicFlag==-1 ) SFl = ((0.927306+(0.00186598*x))+(-4.5141e-06*(x*x)))+(3.21483e-09*(x*(x*x)));
	if( SystematicFlag==+1 ) SFl = ((1.04217+(0.0028073*x))+(-7.16639e-06*(x*x)))+(5.2225e-09*(x*(x*x)));
      }
    if( TaggerName=="CSVL" && fabs(JetEta)>=1.5 && fabs(JetEta)<2.4)
      {
	if( SystematicFlag== 0 ) SFl = ((1.0944+(0.000394694*x))+(-1.31095e-06*(x*x)))+(1.29556e-09*(x*(x*x)));
	if( SystematicFlag==-1 ) SFl = ((1.0255+(0.000220197*x))+(-6.45505e-07*(x*x)))+(6.40579e-10*(x*(x*x)));
	if( SystematicFlag==+1 ) SFl = ((1.1633+(0.000568652*x))+(-1.97487e-06*(x*x)))+(1.95111e-09*(x*(x*x)));
      }
    if( TaggerName=="CSVM" && fabs(JetEta)>=0.0 && fabs(JetEta)<0.8)
      {
	if( SystematicFlag== 0 ) SFl = ((1.13626+(0.00209868*x))+(-5.54303e-06*(x*x)))+(3.9911e-09*(x*(x*x)));
	if( SystematicFlag==-1 ) SFl = ((1.05089+(0.00100001*x))+(-2.44384e-06*(x*x)))+(1.72918e-09*(x*(x*x)));
	if( SystematicFlag==+1 ) SFl = ((1.22164+(0.00319447*x))+(-8.63596e-06*(x*x)))+(6.25306e-09*(x*(x*x)));
      }
    if( TaggerName=="CSVM" && fabs(JetEta)>=0.8 && fabs(JetEta)<1.6)
      {
	if( SystematicFlag== 0 ) SFl = ((1.13291+(0.00128329*x))+(-3.79952e-06*(x*x)))+(3.03032e-09*(x*(x*x)));
	if( SystematicFlag==-1 ) SFl = ((1.04112+(0.000728221*x))+(-2.04996e-06*(x*x)))+(1.64537e-09*(x*(x*x)));
	if( SystematicFlag==+1 ) SFl = ((1.22469+(0.0018366*x))+(-5.54498e-06*(x*x)))+(4.4159e-09*(x*(x*x)));
      }
    if( TaggerName=="CSVM" && fabs(JetEta)>=1.6 && fabs(JetEta)<2.4)
      {
	if( SystematicFlag== 0 ) SFl = ((1.18705+(0.000305854*x))+(-1.86925e-06*(x*x)))+(1.79183e-09*(x*(x*x)));
	if( SystematicFlag==-1 ) SFl = ((1.10587+(-8.23503e-05*x))+(-3.06139e-07*(x*x)))+(2.38667e-10*(x*(x*x)));
	if( SystematicFlag==+1 ) SFl = ((1.26821+(0.000693404*x))+(-3.43071e-06*(x*x)))+(3.34622e-09*(x*(x*x)));
      }
    if( TaggerName=="CSVT" && fabs(JetEta)>=0.0 && fabs(JetEta)<2.4)
      {
	if( SystematicFlag== 0 ) SFl = ((1.08603+(0.0027994*x))+(-8.44182e-06*(x*x)))+(6.847e-09*(x*(x*x)));
	if( SystematicFlag==-1 ) SFl = ((1.02837+(0.00104078*x))+(-3.81136e-06*(x*x)))+(3.43028e-09*(x*(x*x)));
	if( SystematicFlag==+1 ) SFl = ((1.14368+(0.00455363*x))+(-1.30615e-05*(x*x)))+(1.0264e-08*(x*(x*x)));
      }
    if( TaggerName=="JPL" && fabs(JetEta)>=0.0 && fabs(JetEta)<0.5)
      {
	if( SystematicFlag== 0 ) SFl = ((1.12238+(0.00152486*x))+(-3.2873e-06*(x*x)))+(2.17918e-09*(x*(x*x)));
	if( SystematicFlag==-1 ) SFl = ((1.02366+(0.0012007*x))+(-2.45347e-06*(x*x)))+(1.58906e-09*(x*(x*x)));
	if( SystematicFlag==+1 ) SFl = ((1.22108+(0.00184736*x))+(-4.11792e-06*(x*x)))+(2.76952e-09*(x*(x*x)));
      }
    if( TaggerName=="JPL" && fabs(JetEta)>=0.5 && fabs(JetEta)<1.0)
      {
	if( SystematicFlag== 0 ) SFl = ((1.10302+(0.000874045*x))+(-1.99863e-06*(x*x)))+(1.39584e-09*(x*(x*x)));
	if( SystematicFlag==-1 ) SFl = ((1.00257+(0.000719596*x))+(-1.5641e-06*(x*x)))+(1.07029e-09*(x*(x*x)));
	if( SystematicFlag==+1 ) SFl = ((1.20346+(0.00102756*x))+(-2.43131e-06*(x*x)))+(1.72172e-09*(x*(x*x)));
      }
    if( TaggerName=="JPL" && fabs(JetEta)>=1.0 && fabs(JetEta)<1.5)
      {
	if( SystematicFlag== 0 ) SFl = ((1.06244+(0.00149626*x))+(-3.55121e-06*(x*x)))+(2.51004e-09*(x*(x*x)));
	if( SystematicFlag==-1 ) SFl = ((0.973582+(0.00122407*x))+(-2.81096e-06*(x*x)))+(1.94803e-09*(x*(x*x)));
	if( SystematicFlag==+1 ) SFl = ((1.15128+(0.0017669*x))+(-4.28856e-06*(x*x)))+(3.07303e-09*(x*(x*x)));
      }
    if( TaggerName=="JPL" && fabs(JetEta)>=1.5 && fabs(JetEta)<2.4)
      {
	if( SystematicFlag== 0 ) SFl = ((1.10432+(0.00129851*x))+(-3.26353e-06*(x*x)))+(2.32516e-09*(x*(x*x)));
	if( SystematicFlag==-1 ) SFl = ((1.00996+(0.00109673*x))+(-2.66726e-06*(x*x)))+(1.84117e-09*(x*(x*x)));
	if( SystematicFlag==+1 ) SFl = ((1.19864+(0.00149937*x))+(-3.85934e-06*(x*x)))+(2.81243e-09*(x*(x*x)));
      }
    if( TaggerName=="JPM" && fabs(JetEta)>=0.0 && fabs(JetEta)<0.8)
      {
	if( SystematicFlag== 0 ) SFl = ((1.0075+(0.00257791*x))+(-5.91599e-06*(x*x)))+(4.142e-09*(x*(x*x)));
	if( SystematicFlag==-1 ) SFl = ((0.867433+(0.00177435*x))+(-3.62606e-06*(x*x)))+(2.46206e-09*(x*(x*x)));
	if( SystematicFlag==+1 ) SFl = ((1.14757+(0.00337847*x))+(-8.19968e-06*(x*x)))+(5.82223e-09*(x*(x*x)));
      }
    if( TaggerName=="JPM" && fabs(JetEta)>=0.8 && fabs(JetEta)<1.6)
      {
	if( SystematicFlag== 0 ) SFl = ((1.02335+(0.00274819*x))+(-7.08829e-06*(x*x)))+(5.44469e-09*(x*(x*x)));
	if( SystematicFlag==-1 ) SFl = ((0.903546+(0.0017276*x))+(-3.96021e-06*(x*x)))+(3.00594e-09*(x*(x*x)));
	if( SystematicFlag==+1 ) SFl = ((1.14313+(0.00376554*x))+(-1.02094e-05*(x*x)))+(7.88496e-09*(x*(x*x)));
      }
    if( TaggerName=="JPM" && fabs(JetEta)>=1.6 && fabs(JetEta)<2.4)
      {
	if( SystematicFlag== 0 ) SFl = ((0.898963+(0.00428943*x))+(-1.12357e-05*(x*x)))+(8.35894e-09*(x*(x*x)));
	if( SystematicFlag==-1 ) SFl = ((0.804719+(0.00322419*x))+(-7.66523e-06*(x*x)))+(5.18187e-09*(x*(x*x)));
	if( SystematicFlag==+1 ) SFl = ((0.993109+(0.00535115*x))+(-1.48031e-05*(x*x)))+(1.15468e-08*(x*(x*x)));
      }
    if( TaggerName=="JPT" && fabs(JetEta)>=0.0 && fabs(JetEta)<2.4)
      {
	if( SystematicFlag== 0 ) SFl = ((0.887609+(0.00411151*x))+(-1.10861e-05*(x*x)))+(8.50678e-09*(x*(x*x)));
	if( SystematicFlag==-1 ) SFl = ((0.681054+(0.00325624*x))+(-8.17194e-06*(x*x)))+(6.14789e-09*(x*(x*x)));
	if( SystematicFlag==+1 ) SFl = ((1.09416+(0.00496121*x))+(-1.39877e-05*(x*x)))+(1.08667e-08*(x*(x*x)));
      }
    if( TaggerName=="TCHPT" && fabs(JetEta)>=0.0 && fabs(JetEta)<2.4)
      {
	if( SystematicFlag== 0 ) SFl = ((1.25209+(0.00136201*x))+(-3.90275e-06*(x*x)))+(3.1283e-09*(x*(x*x)));
	if( SystematicFlag==-1 ) SFl = ((1.07448+(0.000816207*x))+(-2.38176e-06*(x*x)))+(1.97453e-09*(x*(x*x)));
	if( SystematicFlag==+1 ) SFl = ((1.42974+(0.00190546*x))+(-5.41793e-06*(x*x)))+(4.28148e-09*(x*(x*x)));
      }

    // End of definition of functions from SF_12D ---------------
    
  } 
  
  return SFl;

}

float BTagSFUtil::ScaleFactorJet(int JetFlavor, float JetPt, float JetEta, int SystematicFlag, TString DataPeriod) {
  
  float SF = -1.;
 
  if (JetPt<20.) { cout << "SFb is not available for jet pt<20 GeV" << endl; return -1.; }

  if (abs(JetFlavor)==5) SF = ScaleFactorB(JetPt, SystematicFlag);
  else if (abs(JetFlavor)==4) SF = ScaleFactorB(JetPt, 2*SystematicFlag);
  else SF = ScaleFactorLight(JetPt, JetEta, SystematicFlag, DataPeriod);

  //if (SF==-1.) cout << "Jet parameter out of BTV prescriptions" << endl;

  return SF;

}

float BTagSFUtil::TagEfficiencyB(float JetPt, float JetEta) {

  if (TaggerName=="CSVL") {
    if (30<=JetPt && JetPt<50) {
      if (JetEta>=2.2) return 0.809008; // 3312.87 0.682939
      else if (JetEta>=2) return 0.852787; // 3060.58 0.64046
      else if (JetEta>=1.8) return 0.858175; // 3707.54 0.572956
      else if (JetEta>=1.6) return 0.848624; // 4706.23 0.522456
      else if (JetEta>=1.4) return 0.846757; // 5615.71 0.480692
      else if (JetEta>=1.2) return 0.860785; // 6459.96 0.430701
      else if (JetEta>=1) return 0.838642; // 7239.39 0.432347
      else if (JetEta>=0.8) return 0.856607; // 8095.13 0.389532
      else if (JetEta>=0.6) return 0.858616; // 8676.53 0.374048
      else if (JetEta>=0.4) return 0.856543; // 9048.09 0.368517
      else if (JetEta>=0.2) return 0.850265; // 9378.28 0.368449
      else if (JetEta>=0) return 0.837786; // 9619.12 0.375875
    } else if (50<=JetPt && JetPt<80) {
      if (JetEta>=2.2) return 0.787633; // 2567.98 0.807067
      else if (JetEta>=2) return 0.83771; // 2598.2 0.723363
      else if (JetEta>=1.8) return 0.832037; // 3464.23 0.635147
      else if (JetEta>=1.6) return 0.833712; // 4794.31 0.537744
      else if (JetEta>=1.4) return 0.833563; // 6163.37 0.474444
      else if (JetEta>=1.2) return 0.841055; // 7730.32 0.41585
      else if (JetEta>=1) return 0.842563; // 9031.9 0.383235
      else if (JetEta>=0.8) return 0.846668; // 10456.4 0.352356
      else if (JetEta>=0.6) return 0.855472; // 11673.4 0.325448
      else if (JetEta>=0.4) return 0.85725; // 12383.7 0.314352
      else if (JetEta>=0.2) return 0.85215; // 13135.8 0.309699
      else if (JetEta>=0) return 0.837702; // 13479.1 0.317593
    } else if (80<=JetPt && JetPt<120) {
      if (JetEta>=2) return 0.769223; // 2105.09 0.918305
      else if (JetEta>=1.8) return 0.816982; // 1671.31 0.945855
      else if (JetEta>=1.6) return 0.828224; // 2412.84 0.767874
      else if (JetEta>=1.4) return 0.831575; // 3437.13 0.638345
      else if (JetEta>=1.2) return 0.837462; // 4820.41 0.531396
      else if (JetEta>=1) return 0.842355; // 5947.54 0.472519
      else if (JetEta>=0.8) return 0.852549; // 7018.78 0.423207
      else if (JetEta>=0.6) return 0.865364; // 8107.63 0.379081
      else if (JetEta>=0.4) return 0.863202; // 8724.63 0.367893
      else if (JetEta>=0.2) return 0.857301; // 9481.23 0.359206
      else if (JetEta>=0) return 0.833508; // 9576.16 0.380677
    } else if (120<=JetPt && JetPt<160) {
      if (JetEta>=1.4) return 0.796303; // 2478.95 0.808905
      else if (JetEta>=1.2) return 0.815799; // 1581.29 0.974838
      else if (JetEta>=1) return 0.829592; // 2047.94 0.830843
      else if (JetEta>=0.8) return 0.853015; // 2451.91 0.715093
      else if (JetEta>=0.6) return 0.849798; // 2855.79 0.668548
      else if (JetEta>=0.4) return 0.870126; // 3139.45 0.599964
      else if (JetEta>=0.2) return 0.850251; // 3397.08 0.612213
      else if (JetEta>=0) return 0.835425; // 3555.59 0.621842
    } else if (160<=JetPt && JetPt<210) {
      if (JetEta>=1) return 0.793887; // 1792.65 0.9554
      else if (JetEta>=0.6) return 0.835242; // 1816.48 0.87039
      else if (JetEta>=0.2) return 0.852439; // 2407.34 0.722851
      else if (JetEta>=-0.2) return 0.838569; // 1332.46 1.00794
    } else if (210<=JetPt && JetPt<100000) {
      if (JetEta>=0.4) return 0.778819; // 1749.14 0.992385
      else if (JetEta>=-0.2) return 0.806058; // 1133.46 1.1744
    }
  }

  if (TaggerName=="CSVM") {
    if (30<=JetPt && JetPt<50) {
      if (JetEta>=2.2) return 0.457318; // 3312.87 0.865525
      else if (JetEta>=2) return 0.581291; // 3060.58 0.891767
      else if (JetEta>=1.8) return 0.621905; // 3707.54 0.796378
      else if (JetEta>=1.6) return 0.608839; // 4706.23 0.711365
      else if (JetEta>=1.4) return 0.609522; // 5615.71 0.651015
      else if (JetEta>=1.2) return 0.615207; // 6459.96 0.605354
      else if (JetEta>=1) return 0.636686; // 7239.39 0.565266
      else if (JetEta>=0.8) return 0.676667; // 8095.13 0.519877
      else if (JetEta>=0.6) return 0.678029; // 8676.53 0.501602
      else if (JetEta>=0.4) return 0.68552; // 9048.09 0.488122
      else if (JetEta>=0.2) return 0.686303; // 9378.28 0.479128
      else if (JetEta>=0) return 0.665883; // 9619.12 0.480929
    } else if (50<=JetPt && JetPt<80) {
      if (JetEta>=2.2) return 0.468871; // 2567.98 0.984761
      else if (JetEta>=2) return 0.610662; // 2598.2 0.956594
      else if (JetEta>=1.8) return 0.648496; // 3464.23 0.811176
      else if (JetEta>=1.6) return 0.651667; // 4794.31 0.688093
      else if (JetEta>=1.4) return 0.645824; // 6163.37 0.609197
      else if (JetEta>=1.2) return 0.666929; // 7730.32 0.536055
      else if (JetEta>=1) return 0.692456; // 9031.9 0.485579
      else if (JetEta>=0.8) return 0.720003; // 10456.4 0.439089
      else if (JetEta>=0.6) return 0.737833; // 11673.4 0.407071
      else if (JetEta>=0.4) return 0.741772; // 12383.7 0.393289
      else if (JetEta>=0.2) return 0.731069; // 13135.8 0.386875
      else if (JetEta>=0) return 0.723047; // 13479.1 0.385439
    } else if (80<=JetPt && JetPt<120) {
      if (JetEta>=1.8) return 0.59996; // 3778.31 0.797011
      else if (JetEta>=1.6) return 0.681846; // 2412.84 0.948194
      else if (JetEta>=1.4) return 0.680889; // 3437.13 0.79508
      else if (JetEta>=1.2) return 0.696981; // 4820.41 0.661916
      else if (JetEta>=1) return 0.730068; // 5947.54 0.575626
      else if (JetEta>=0.8) return 0.7469; // 7018.78 0.518975
      else if (JetEta>=0.6) return 0.77531; // 8107.63 0.463535
      else if (JetEta>=0.4) return 0.765666; // 8724.63 0.453486
      else if (JetEta>=0.2) return 0.766997; // 9481.23 0.434155
      else if (JetEta>=0) return 0.741348; // 9576.16 0.447479
    } else if (120<=JetPt && JetPt<160) {
      if (JetEta>=1.4) return 0.639172; // 2478.95 0.964551
      else if (JetEta>=1) return 0.701581; // 3629.46 0.759506
      else if (JetEta>=0.8) return 0.758373; // 2451.91 0.864494
      else if (JetEta>=0.6) return 0.74889; // 2855.79 0.811479
      else if (JetEta>=0.4) return 0.77955; // 3139.45 0.739861
      else if (JetEta>=0.2) return 0.753272; // 3397.08 0.73966
      else if (JetEta>=0) return 0.736507; // 3555.59 0.738783
    } else if (160<=JetPt && JetPt<210) {
      if (JetEta>=0.8) return 0.66771; // 2630.44 0.918414
      else if (JetEta>=0.4) return 0.740332; // 2140.21 0.947752
      else if (JetEta>=0) return 0.739091; // 2579.1 0.864687
    } else if (210<=JetPt && JetPt<100000) {
      if (JetEta>=0.2) return 0.64097; // 2312.74 0.997519
      else if (JetEta>=-0.2) return 0.704848; // 569.568 1.91116
    }
  }

  cout << "BTagSFUtil: Jet pt (" << JetPt << ") or jet eta (" << JetEta << ") out of range!" << endl;

  if (TaggerName=="CSVL") return 0.80;
  if (TaggerName=="CSVM") return 0.65;

  cout << "BTagSFUtil: Tagger (" << TaggerName << ") not valid!" << endl;

  return 0.0;

}

float BTagSFUtil::TagEfficiencyC(float JetPt, float JetEta) {

  if (TaggerName=="CSVL") {
    if (30<=JetPt && JetPt<50) {
      if (JetEta>=0) return 0.473677; // 2671.38 0.966049
    } else if (50<=JetPt && JetPt<80) {
      if (JetEta>=-0.2) return 0.392501; // 1795.58 1.15237
    } else if (80<=JetPt && JetPt<120) {
      if (JetEta>=-0.2) return 0.359924; // 1061.43 1.47324
    } else if (120<=JetPt && JetPt<160) {
      if (JetEta>=-0.2) return 0.37042; // 484.451 2.19406
    } else if (160<=JetPt && JetPt<210) {
      if (JetEta>=-0.2) return 0.359324; // 302.053 2.76071
    } else if (210<=JetPt && JetPt<100000) {
      if (JetEta>=-0.2) return 0.345628; // 276.093 2.86213
    }
  }

  if (TaggerName=="CSVM") {
    if (30<=JetPt && JetPt<50) {
      if (JetEta>=1) return 0.137101; // 1367.57 0.93009
      else if (JetEta>=-0.2) return 0.172955; // 1302.01 1.04815
    } else if (50<=JetPt && JetPt<80) {
      if (JetEta>=0.4) return 0.163397; // 1413.84 0.983288
      else if (JetEta>=-0.2) return 0.173224; // 381.629 1.93721
    } else if (80<=JetPt && JetPt<120) {
      if (JetEta>=-0.2) return 0.18961; // 1061.43 1.20318
    } else if (120<=JetPt && JetPt<160) {
      if (JetEta>=-0.2) return 0.188993; // 484.451 1.77873
    } else if (160<=JetPt && JetPt<210) {
      if (JetEta>=-0.2) return 0.197597; // 302.053 2.2911
    } else if (210<=JetPt && JetPt<100000) {
      if (JetEta>=-0.2) return 0.204649; // 276.093 2.42804
    }
  }
  
  cout << "BTagSFUtil: Jet pt (" << JetPt << ") or jet eta (" << JetEta << ") out of range!" << endl;

  if (TaggerName=="CSVL") return 0.35;
  if (TaggerName=="CSVM") return 0.16;

  cout << "BTagSFUtil: Tagger (" << TaggerName << ") not valid!" << endl;

  return 0.0;

}

float BTagSFUtil::TagEfficiencyLight(float JetPt, float JetEta) {

  if (TaggerName=="CSVL") {
    if (30<=JetPt && JetPt<50) {
      if (JetEta>=2.2) return 0.420442; // 3926.52 0.787767
      else if (JetEta>=2) return 0.338665; // 2983.58 0.866418
      else if (JetEta>=1.8) return 0.247388; // 3171.73 0.766173
      else if (JetEta>=1.6) return 0.238397; // 3493.16 0.720951
      else if (JetEta>=1.4) return 0.238432; // 3989.68 0.674632
      else if (JetEta>=1.2) return 0.240449; // 4269.08 0.654067
      else if (JetEta>=1) return 0.148722; // 4592.34 0.525057
      else if (JetEta>=0.8) return 0.165212; // 4739.99 0.539412
      else if (JetEta>=0.6) return 0.144527; // 5046.85 0.494956
      else if (JetEta>=0.4) return 0.148472; // 5155.2 0.495221
      else if (JetEta>=0.2) return 0.127898; // 4976.4 0.473432
      else if (JetEta>=0) return 0.127603; // 5197.03 0.462818
    } else if (50<=JetPt && JetPt<80) {
      if (JetEta>=2.2) return 0.321751; // 2480.17 0.938024
      else if (JetEta>=2) return 0.238807; // 1867.7 0.986545
      else if (JetEta>=1.8) return 0.137103; // 1896.63 0.78979
      else if (JetEta>=1.6) return 0.123762; // 2209.89 0.700518
      else if (JetEta>=1.4) return 0.14158; // 2504.66 0.696588
      else if (JetEta>=1.2) return 0.113777; // 2764.43 0.603942
      else if (JetEta>=1) return 0.0912294; // 2834.32 0.540841
      else if (JetEta>=0.8) return 0.0903294; // 2960.95 0.526794
      else if (JetEta>=0.6) return 0.101207; // 2992.92 0.5513
      else if (JetEta>=0.4) return 0.0799342; // 3249.78 0.475717
      else if (JetEta>=0.2) return 0.0784796; // 3148.24 0.479288
      else if (JetEta>=0) return 0.0772237; // 3205.01 0.471529
    } else if (80<=JetPt && JetPt<120) {
      if (JetEta>=2) return 0.252252; // 2326.08 0.900497
      else if (JetEta>=1.6) return 0.114915; // 2269.35 0.669469
      else if (JetEta>=1.4) return 0.110387; // 1359.65 0.849856
      else if (JetEta>=1.2) return 0.0830715; // 1533.19 0.704848
      else if (JetEta>=1) return 0.0778623; // 1585.36 0.672974
      else if (JetEta>=0.8) return 0.0768658; // 1608.92 0.664097
      else if (JetEta>=0.6) return 0.0676013; // 1712.71 0.606647
      else if (JetEta>=0.4) return 0.0830071; // 1676.9 0.673732
      else if (JetEta>=0.2) return 0.073603; // 1714.09 0.63071
      else if (JetEta>=0) return 0.0620921; // 1822.07 0.565348
    } else if (120<=JetPt && JetPt<160) {
      if (JetEta>=1.6) return 0.181298; // 2011.89 0.858929
      else if (JetEta>=1.2) return 0.0864793; // 1317.64 0.774315
      else if (JetEta>=1) return 0.0666236; // 781.072 0.892271
      else if (JetEta>=0.8) return 0.0800225; // 756.503 0.986482
      else if (JetEta>=0.6) return 0.0664807; // 814.199 0.87306
      else if (JetEta>=0.4) return 0.0665918; // 817.127 0.872171
      else if (JetEta>=0.2) return 0.0609753; // 725.88 0.888143
      else if (JetEta>=0) return 0.0751063; // 808.285 0.927047
    } else if (160<=JetPt && JetPt<210) {
      if (JetEta>=1.4) return 0.164064; // 1524.31 0.948543
      else if (JetEta>=1) return 0.0921369; // 838.988 0.998502
      else if (JetEta>=0.6) return 0.0978261; // 908.33 0.985714
      else if (JetEta>=0.2) return 0.0780631; // 955.181 0.868022
      else if (JetEta>=-0.2) return 0.0832897; // 501.078 1.23441
    } else if (210<=JetPt && JetPt<100000) {
      if (JetEta>=1.2) return 0.129322; // 1701.89 0.81339
      else if (JetEta>=0.8) return 0.0753734; // 962.893 0.850753
      else if (JetEta>=0.4) return 0.0895906; // 923.429 0.939827
      else if (JetEta>=0) return 0.0755069; // 1017.49 0.828285
    }
  }
  
  if (TaggerName=="CSVM") {
    if (30<=JetPt && JetPt<50) {
      if (JetEta>=2.2) return 0.0218702; // 3926.52 0.233411
      else if (JetEta>=2) return 0.0252778; // 2983.58 0.28737
      else if (JetEta>=1.8) return 0.0187535; // 3171.73 0.24087
      else if (JetEta>=1.6) return 0.0207562; // 3493.16 0.241219
      else if (JetEta>=1.4) return 0.0172815; // 3989.68 0.206318
      else if (JetEta>=1.2) return 0.0127035; // 4269.08 0.171403
      else if (JetEta>=1) return 0.0157108; // 4592.34 0.183503
      else if (JetEta>=0.8) return 0.0164909; // 4739.99 0.184979
      else if (JetEta>=0.6) return 0.0132922; // 5046.85 0.161207
      else if (JetEta>=0.4) return 0.0145052; // 5155.2 0.16652
      else if (JetEta>=0.2) return 0.0105695; // 4976.4 0.144964
      else if (JetEta>=0) return 0.0133768; // 5197.03 0.159358
    } else if (50<=JetPt && JetPt<80) {
      if (JetEta>=2.2) return 0.0121253; // 2480.17 0.219764
      else if (JetEta>=2) return 0.0195525; // 1867.7 0.320375
      else if (JetEta>=1.8) return 0.0159969; // 1896.63 0.288088
      else if (JetEta>=1.6) return 0.0163013; // 2209.89 0.269375
      else if (JetEta>=1.4) return 0.014961; // 2504.66 0.242567
      else if (JetEta>=1.2) return 0.0082208; // 2764.43 0.171736
      else if (JetEta>=1) return 0.0112669; // 2834.32 0.198252
      else if (JetEta>=0.8) return 0.0136138; // 2960.95 0.21296
      else if (JetEta>=0.6) return 0.0154717; // 2992.92 0.225598
      else if (JetEta>=0.4) return 0.00991515; // 3249.78 0.173804
      else if (JetEta>=0.2) return 0.0119935; // 3148.24 0.194008
      else if (JetEta>=0) return 0.0114284; // 3205.01 0.187751
    } else if (80<=JetPt && JetPt<120) {
      if (JetEta>=2.2) return 0.0171949; // 1340.69 0.355034
      else if (JetEta>=2) return 0.0208581; // 985.459 0.45524
      else if (JetEta>=1.8) return 0.0210796; // 1078.96 0.437323
      else if (JetEta>=1.6) return 0.00915257; // 1190.28 0.276026
      else if (JetEta>=1.4) return 0.0140617; // 1359.65 0.319322
      else if (JetEta>=1.2) return 0.00958679; // 1533.19 0.248855
      else if (JetEta>=1) return 0.0140639; // 1585.36 0.295742
      else if (JetEta>=0.8) return 0.0137454; // 1608.92 0.290273
      else if (JetEta>=0.6) return 0.0116382; // 1712.71 0.259154
      else if (JetEta>=0.4) return 0.0115792; // 1676.9 0.261251
      else if (JetEta>=0.2) return 0.01354; // 1714.09 0.279147
      else if (JetEta>=0) return 0.00933232; // 1822.07 0.225256
    } else if (120<=JetPt && JetPt<160) {
      if (JetEta>=2.2) return 0.0113618; // 548.511 0.452532
      else if (JetEta>=2) return 0.0362423; // 448.5 0.882492
      else if (JetEta>=1.8) return 0.0201319; // 471.407 0.646886
      else if (JetEta>=1.6) return 0.0141313; // 543.408 0.506334
      else if (JetEta>=1.4) return 0.018753; // 651.027 0.53165
      else if (JetEta>=1.2) return 0.0215288; // 666.74 0.56209
      else if (JetEta>=1) return 0.00903814; // 781.072 0.338628
      else if (JetEta>=0.8) return 0.0245823; // 756.503 0.562991
      else if (JetEta>=0.6) return 0.0149224; // 814.199 0.424902
      else if (JetEta>=0.4) return 0.00705541; // 817.127 0.292805
      else if (JetEta>=0.2) return 0.0195207; // 725.88 0.513493
      else if (JetEta>=0) return 0.0126196; // 808.285 0.392629
    } else if (160<=JetPt && JetPt<210) {
      if (JetEta>=2.2) return 0.00722232; // 305.489 0.484469
      else if (JetEta>=2) return 0.0057698; // 254.765 0.474519
      else if (JetEta>=1.8) return 0.0139492; // 284.348 0.695503
      else if (JetEta>=1.6) return 0.0255133; // 315.103 0.88827
      else if (JetEta>=1.4) return 0.00946272; // 364.607 0.507027
      else if (JetEta>=1.2) return 0.00377937; // 400.467 0.306622
      else if (JetEta>=1) return 0.0125736; // 438.595 0.532047
      else if (JetEta>=0.8) return 0.017102; // 468.029 0.599296
      else if (JetEta>=0.6) return 0.0240959; // 440.421 0.730703
      else if (JetEta>=0.4) return 0.0137613; // 472.588 0.535894
      else if (JetEta>=0.2) return 0.0218971; // 482.462 0.666275
      else if (JetEta>=0) return 0.0172466; // 501.078 0.581596
    } else if (210<=JetPt && JetPt<100000) {
      if (JetEta>=2) return 0.0246663; // 358.693 0.818968
      else if (JetEta>=1.6) return 0.0373008; // 519.881 0.831098
      else if (JetEta>=1.4) return 0.0338009; // 383.362 0.92298
      else if (JetEta>=1.2) return 0.0128993; // 441.087 0.53728
      else if (JetEta>=1) return 0.0189693; // 454.311 0.640015
      else if (JetEta>=0.8) return 0.0177886; // 508.606 0.586115
      else if (JetEta>=0.6) return 0.0139127; // 474.895 0.537483
      else if (JetEta>=0.4) return 0.0145154; // 448.404 0.564813
      else if (JetEta>=0.2) return 0.0210017; // 501.926 0.640028
      else if (JetEta>=0) return 0.0180327; // 515.528 0.586074
    }
  }
  
  cout << "BTagSFUtil: Jet pt (" << JetPt << ") or jet eta (" << JetEta << ") out of range!" << endl;

  if (TaggerName=="CSVL") return 0.10;
  if (TaggerName=="CSVM") return 0.02;

  cout << "BTagSFUtil: Tagger (" << TaggerName << ") not valid!" << endl;

  return 0.0;

}

float BTagSFUtil::JetTagEfficiency(int JetFlavor, float JetPt, float JetEta) {

  if (abs(JetFlavor)==5) return TagEfficiencyB(JetPt, JetEta);
  else if (abs(JetFlavor)==4) return TagEfficiencyC(JetPt, JetEta);
  else return TagEfficiencyLight(JetPt, JetEta);

}

bool BTagSFUtil::IsTagged(float JetDiscriminant, int JetFlavor, float JetPt, float JetEta, int SystematicFlag, TString DataPeriod) {
  
  bool isBTagged = JetDiscriminant>TaggerCut;

  if (JetFlavor==-999999) return isBTagged; // Data: no correction needed

  bool newBTag = isBTagged;

  float Btag_SF = ScaleFactorJet(JetFlavor, JetPt, JetEta, SystematicFlag, DataPeriod);

  if (Btag_SF == 1) return newBTag; //no correction needed 

  //throw die
  float coin = rand_->Uniform(1.);    

  //cout << "##### " << coin << endl;
 
  if(Btag_SF > 1){  // use this if SF>1

    if( !isBTagged ) {

      float Btag_eff = JetTagEfficiency(JetFlavor, JetPt, fabs(JetEta));

      //fraction of jets that need to be upgraded
      float mistagPercent = (1.0 - Btag_SF) / (1.0 - (Btag_SF/Btag_eff) );

      //upgrade to tagged
      if( coin < mistagPercent ) {newBTag = true;}
    }

  }else{  // use this if SF<1
      
    //downgrade tagged to untagged
    if( isBTagged && coin > Btag_SF ) {newBTag = false;}

  }

  return newBTag;

}

