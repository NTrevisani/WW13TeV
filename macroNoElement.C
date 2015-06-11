//Clone the latino tree into a new one without the PAF-hated branches
//run typing:  root -l macroNoElement.C                                                                                                                  

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

void NoElement(TString name){

  TFile* DY = new TFile(name,"read");

  TTree *tDark = (TTree*) DY -> Get("latino");
  cout<<"ciao"<<endl;

  
  TObjArray *tl = tDark->GetListOfBranches();
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
  newName.Remove(0,38);

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
  //delete tl;
}

void macroNoElement(){

  TString command = "ls /gpfs/csic_projects/cms/piedra/latino/";
  command = command + "*.root> inputNoElement.tmp"; 
  gSystem -> Exec(command); 

  TString command = "ls /gpfs/csic_projects/cms/trevisanin/newLatino/";
  command = command + "> inputHere.tmp";
  gSystem -> Exec(command);
  
  ifstream inFile("inputNoElement.tmp");

  std::string line;
  std::string lineHere;
  
  //run the code only if you don't have already the tree here in the folder
  while(getline(inFile,line)){
    Float_t pair = 0;
    //cout<<line<<endl;
    ifstream inFileHere("inputHere.tmp");
    while(getline(inFileHere,lineHere)){
      TString newLine = line;
      newLine.Remove(0,31);
      TString lineHere_ = lineHere;
      //cout<<newLine<<" "<<lineHere_<<endl;
      if(newLine == lineHere_) pair = 1;
    }
    newLine.Remove(0,7);
    if (pair == 0){
      TString move = "mv ";
      move = move + newLine + " latino/";
      cout<<move<<endl;
      NoElement(line);
      gSystem -> Exec(move);
    }
    else cout<<newLine<<" ya estaba!"<<endl;
  }
}

