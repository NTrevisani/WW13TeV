//Clone the latino tree into a new one without the PAF-hated branches
//run typing:  root -l 'macroNoElement.C("original folder","your folder", copy all files: 1 / copy only missing files: 0)'
//e.g. root -l 'macroNoElement.C("/gpfs/csic_projects/cms/piedra/latino/","/gpfs/csic_projects/cms/trevisanin/newLatino/",0)'

#include "TFile.h"
#include "TH1.h"
#include "TCanvas.h"
#include "TROOT.h"
#include "TTree.h"
#include "TChain.h"
#include "TString.h"
#include "TLegend.h"
#include "TObjArray.h"
#include <algorithm>
#include "TCut.h"
#include <iostream>
#include <fstream>
#include <string>

TString latinoVar = "";
TString darkVar = "";
Float_t range = 1000.;

using namespace std;

int NoElement(TString name, Int_t Long){

  TFile* DY = new TFile(name,"read");

  if (DY -> GetListOfKeys()->Contains("latino") == 0){
    cout<<"cannot find a tree named 'latino'"<<endl;
    return 0;
  }

  TTree *tDark = (TTree*) DY -> Get("latino");
 
  TObjArray *tl = tDark -> GetListOfBranches();

  cout<<"ciao"<<endl;
  TString nBranch = tl->First()->GetName();

  Int_t cont = 0;
  tDark->SetBranchStatus("*",0);

  while (tl -> After(tl->FindObject(nBranch)) != 0){                                                                                           
    nBranch = tl -> At(cont) -> GetName();
    if (nBranch.Contains("v_") || nBranch.Contains( "fCoordinate")) //skip the structured branches
      tDark->SetBranchStatus(nBranch,0);
    else 
      tDark->SetBranchStatus(nBranch,1);
      ++cont;    
  }
  
  TString newName = name;
  newName.Remove(0,Long);//38);

  cout<<newName<<endl;

  TFile *newfile = new TFile(newName,"recreate");
  
  cout<<"trabajando..."<<endl;
  TTree *newtree = tDark->CloneTree();
  newtree -> SetDirectory(0);

  newtree->Print();
  newtree->Write();
  newfile->Close();

  delete tDark;
  delete newtree;
  delete newfile;

  return 1;
}

void macroNoElement(TString startFolder = "/gpfs/csic_projects/cms/piedra/latino/", 
		    TString arrivalFolder = "/gpfs/csic_projects/cms/trevisanin/newLatino/",
		    Int_t copyAll = 0){

  gSystem -> Exec("mkdir " + arrivalFolder);

  TFile *disposable = new TFile(arrivalFolder + "disp.root","recreate");
  disposable->cd();
  disposable->Close();

  TString command = "ls ";///gpfs/csic_projects/cms/piedra/latino/";
  command = command + startFolder;
  command = command + "*.root> inputNoElement.tmp"; 
  gSystem -> Exec(command); 

  command = "ls ";///gpfs/csic_projects/cms/trevisanin/newLatino/";
  command = command + arrivalFolder;
  command = command + "> inputHere.tmp";
  gSystem -> Exec(command);
  
  ifstream inFile("inputNoElement.tmp");

  std::string line;
  std::string lineHere;
  
  //run the code only if you don't already have the file here in the folder (if copyAll = 1 run the code for all the files)
  while(getline(inFile,line)){
    cout<<line<<endl;
    Int_t pair = 1;
    ifstream inFileHere("inputHere.tmp");
    while(getline(inFileHere,lineHere)){
      TString newLine = line;
      newLine.Remove(0.,startFolder.Length());
      TString lineHere_ = lineHere;
      cout<<newLine<<" "<<lineHere_<<endl;
      if (copyAll == 0 && newLine == lineHere_) pair = 0;
    }
    if (pair == 1){
      TString move = "mv ";
      move = move + newLine + " " + arrivalFolder;
      cout<<move<<endl;
      if (NoElement(line,startFolder.Length()) == 1)
      gSystem -> Exec(move);
    }
    else cout<<newLine<<" ya estaba!"<<endl;
  }
  gSystem -> Exec("rm " + arrivalFolder + "disp.root");
}

