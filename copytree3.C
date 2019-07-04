
// copytree3.C

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


//________________________________MAIN SCRIPT______________________________

void copytree3() {
   //Get old file, old tree and set top branch address
   TFile * fD=TFile::Open("/sps/lhcb/marin/lb2pkgamma/radiative2hG_R14S21_MagDown_9.root");
   TTree * tD = (TTree*)fD->Get("pkGTuple/DecayTree");
   //TFile * fM=TFile::Open("/data/lhcb/marin/lb2pkgamma/MC/2012/15102203/2hG-S21/radiative2hG_MC2012-Lb2L1520gamma_HighPt-15102203-Py8Sim09dReco14c_S21.root");
   //TTree * tM = (TTree*)fM->Get("pkGTupleMC/DecayTree");
   
   //gRoot->cd();

   //Create a new file + a clone of old tree in new file
   TFile *newfile = new TFile("/sps/lhcb/volle/MagUp/MagDown_9.root","recreate");
   TTree * treecopy = tD->CopyTree( "pplus_ProbNNp>0.2 && Kminus_ProbNNk>0.2 && B_PT>4000 && Lambda_1520_0_PT>1500 && gamma_PT>3000 && pplus_PT>1000 && B_FDCHI2_OWNPV>100 && pplus_IPCHI2_OWNPV>50 && Kminus_IPCHI2_OWNPV>40 && (B_L0PhotonDecision_TOS || B_L0ElectronDecision_TOS) && (B_Hlt1TrackAllL0Decision_TOS || B_Hlt1TrackPhotonDecision_TOS) && (B_Hlt2TopoRad2BodyBBDTDecision_TOS || B_Hlt2TopoRad2plus1BodyBBDTDecision_TOS || B_Hlt2RadiativeTopoTrackDecision_TOS || B_Hlt2RadiativeTopoPhotonDecision_TOS)"); //Data

   //TFile *newfile = new TFile("/sps/lhcb/volle/MCsignal_TrigSel.root","recreate");
   //TTree * treecopy = tM->CopyTree( "pplus_ProbNNp>0.2 && Kminus_ProbNNk>0.2 && B_PT>4000 && Lambda_1520_0_PT>1500 && gamma_PT>3000 && pplus_PT>1000 && B_FDCHI2_OWNPV>100 && pplus_IPCHI2_OWNPV>50 && Kminus_IPCHI2_OWNPV>40 && (B_L0PhotonDecision_TOS || B_L0ElectronDecision_TOS) && (B_Hlt1TrackAllL0Decision_TOS || B_Hlt1TrackPhotonDecision_TOS) && (B_Hlt2TopoRad2BodyBBDTDecision_TOS || B_Hlt2TopoRad2plus1BodyBBDTDecision_TOS || B_Hlt2RadiativeTopoTrackDecision_TOS || B_Hlt2RadiativeTopoPhotonDecision_TOS) && B_BKGCAT==0"); //MC
   newfile->Write();
   newfile->Close();
   //return 0;
}
