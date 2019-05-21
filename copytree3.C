
//Copytree-file

#include "TFile.h"
#include "TTree.h"
#include "TChain.h"
//#include "TBrowser.h"
//#include "TH2.h"
//#include "TRandom.h"
//#include "TClassTable.h"
#include "TSystem.h"
#include "TROOT.h"

void copytree3() {
   //Get old file, old tree and set top branch address
   TFile * fMC=TFile::Open("/data/lhcb/marin/lb2pkgamma/MC/2012/15102203/2hG-S21/radiative2hG_MC2012-Lb2L1520gamma_HighPt-15102203-Py8Sim09dReco14c_S21.root");
   TTree * tMC = (TTree*)fMC->Get("pkGTupleMC/DecayTree");
   //gRoot->cd();

   //Create a new file + a clone of old tree in new file
   TFile *newfile = new TFile("/exp/LHCb/volle/MCsignal.root","recreate");
   TTree * treecopy = tMC->CopyTree("B_BKGCAT==0");
   newfile->Write();
   newfile->Close();
   //return 0;
}
