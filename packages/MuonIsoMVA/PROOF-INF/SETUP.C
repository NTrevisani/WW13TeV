Int_t SETUP() {
  // Load library
  //gROOT->ProcessLine(".L MuonIsoMVA.C+");
  //cout << "Loading MuonIsoMVA_C.so..." << endl;

  if (gSystem->Load("MuonIsoMVA_C") == -1)
     return -1;

  /*
  // Enlarge include path (for MuonIsoMVA.h)
  TString pwdpath = "-I";
  pwdpath+=gSystem->pwd(); //path to local dir
  TString curpath = gSystem->GetIncludePath(); //current path

  //Some info
  if (gProofDebugLevel > 0) {
    TString infomess="Current PWD path is \"";
    infomess+=pwdpath;
    infomess+="\"";
    Info("SETUP.C",infomess.Data());
    infomess = "Current include path is \"";
    infomess+=curpath;
    infomess+="\"";
    Info("SETUP.C",infomess.Data());
    //
  }

  if (!curpath.Contains(pwdpath)) {
    gSystem->AddIncludePath(pwdpath);
    
    if (gProofDebugLevel > 0)
      cout << ">> Include path set to: \"" << gSystem->GetIncludePath() 
	   << "\"" << endl;
  }
  */

  return 0;
}
