
// copytree3.C

#include "TFile.h"
#include "TTree.h"
#include "TChain.h"
//#include "TBrowser.h"
//#include "TH2.h"
//#include "TRandom.h"
//#include "TClassTable.h"
#include "TSystem.h"
#include "TROOT.h"

//________________________________MAIN SCRIPT______________________________

void copytree3() {
   //Get old file, old tree and set top branch address
   TFile * fD=TFile::Open("/data/lhcb/marin/lb2pkgamma/Data/2012/2hG-S21/radiative2hG_R14S21_MagUp_tmp.root");
   TTree * tD = (TTree*)fD->Get("pkGTuple/DecayTree");
   //gRoot->cd();

   //Create a new file + a clone of old tree in new file
   TFile *newfile = new TFile("/exp/LHCb/volle/Data_RS.root","recreate");
   TTree * treecopy = tD->CopyTree( "pplus_ProbNNp>0.2 && Kminus_ProbNNk>0.2 && B_PT>4000 && Lambda_1520_0_PT>1500 && gamma_PT>3000 && pplus_PT>1000 && B_FDCHI2_OWNPV>100 && pplus_IPCHI2_OWNPV>50 && Kminus_IPCHI2_OWNPV>40 && B_M>6120");
   newfile->Write();
   newfile->Close();
   //return 0;
}
