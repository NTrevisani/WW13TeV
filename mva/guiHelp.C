//root -l

void guiHelp(){
gROOT->LoadMacro("$ROOTSYS/tmva/test/TMVAGui.C");
TMVAGui("myAnalysisFile.root");
return;
}