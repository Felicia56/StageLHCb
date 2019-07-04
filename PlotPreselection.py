#!/usr/bin/env python

'''
Project        : Plot Histograms for left and right sideband
Script name    : PlotHistograms.py 
Author         : Felicia Volle
'''

# Python and ROOT import
import sys    # exit
import time   # time accounting
import getopt # command line parser
import os, pickle
from ROOT import gStyle,gROOT,TLegend,TPad,TCanvas,TH1F,TH2F,TFile,TDirectory,TCut,THStack,TLatex,TTree


#________________________________SET VARIABLES________________________________


TFileMCName = "/data/lhcb/marin/lb2pkgamma/MC/2012/15102203/2hG-S21/radiative2hG_MC2012-Lb2L1520gamma_HighPt-15102203-Py8Sim09dReco14c_S21.root"##"/sps/lhcb/volle/MCsignal_TrigSel.root"
TFileDataName = "/sps/lhcb/marin/lb2pkgamma/radiative2hG_R14S21_MagUp.root"#"/sps/lhcb/volle/Data_all_WithoutCut.root"#"/sps/lhcb/volle/Data_all.root"#"/data/lhcb/marin/lb2pkgamma/Data/2012/2hG-S21/radiative2hG_R14S21_MagUp_tmp.root"#"/sps/lhcb/volle/Data_full_TrigSel.root"
MCTreeName = "pkGTupleMC/DecayTree"#
DataTreeName = "pkGTuple/DecayTree"#"DecayTree"#


CutMC = TCut("B_BKGCAT==0 && (B_L0PhotonDecision_TOS || B_L0ElectronDecision_TOS) && (B_Hlt1TrackAllL0Decision_TOS || B_Hlt1TrackPhotonDecision_TOS) && (B_Hlt2TopoRad2BodyBBDTDecision_TOS || B_Hlt2TopoRad2plus1BodyBBDTDecision_TOS || B_Hlt2RadiativeTopoTrackDecision_TOS || B_Hlt2RadiativeTopoPhotonDecision_TOS) && pplus_ProbNNp>0.2")#"pplus_ProbNNp>0.2 && Kminus_ProbNNk>0.2 && B_PT>4000 && Lambda_1520_0_PT>1500 && gamma_PT>3000 && pplus_PT>1000 && B_FDCHI2_OWNPV>100 && pplus_IPCHI2_OWNPV>50 && Kminus_IPCHI2_OWNPV>40 && B_BKGCAT==0")

CutData_D = TCut("(B_L0PhotonDecision_TOS || B_L0ElectronDecision_TOS) && (B_Hlt1TrackAllL0Decision_TOS || B_Hlt1TrackPhotonDecision_TOS) && (B_Hlt2TopoRad2BodyBBDTDecision_TOS || B_Hlt2TopoRad2plus1BodyBBDTDecision_TOS || B_Hlt2RadiativeTopoTrackDecision_TOS || B_Hlt2RadiativeTopoPhotonDecision_TOS) && pplus_ProbNNp>0.2")#"pplus_ProbNNp>0.2 && Kminus_ProbNNk>0.2 && B_PT>4000 && Lambda_1520_0_PT>1500 && gamma_PT>3000 && pplus_PT>1000 && B_FDCHI2_OWNPV>100 && pplus_IPCHI2_OWNPV>50 && Kminus_IPCHI2_OWNPV>40 && B_M<5120")


#_____________________________DEFINITION OF FUNCTIONS________________________ 


def setupHistograms():
    histograms = [
        ["pplus_ProbNNp","Prob_{NN}(p)","Entries","pplus_ProbNNp>>h(100,0,1)"],

        ]
    return histograms

def GetAndDrawHistograms(DataFileName,DataTreeName,MCFileName,MCTreeName,noStack=True):
    fMC = TFile(MCFileName)
    tMC = fMC.Get(MCTreeName)
    fData = TFile(DataFileName)
    tData = fData.Get(DataTreeName)

    histoListe=setupHistograms()
    HistDic = {}

    for i,h in enumerate(histoListe):
        legende=TLegend(0.7, 0.7, 0.89, 0.89)
        CanvasName = "c"+str(i)
        CanvasName = TCanvas(CanvasName,'Defining histogram size') 
        g=THStack("","")
        
        tMC.Draw(histoListe[i][3],CutMC)
        hist_MC = tMC.GetHistogram().Clone("hist_MC")
        hist_MC.SetLineColor(2)
        hist_MC.SetLineWidth(1)
        #hist_MC.SetMarkerStyle(21)
        legende.AddEntry(hist_MC, "MC", "lp")
        hist_MC.GetXaxis().SetTitle(h[1])
        hist_MC.GetYaxis().SetTitle(h[2])
        hist_MC.Scale(1/hist_MC.Integral())
        MaxMC = hist_MC.GetMaximum()
        g.Add(hist_MC)
        
        tData.Draw(histoListe[i][3],CutData_D)
        hist_D = tData.GetHistogram().Clone("hist_D")
        hist_D.SetLineColor(9)
        hist_D.SetLineWidth(1)
        #hist_D.SetMarkerStyle(21)
        legende.AddEntry(hist_D, "Data", "lp")
        #hist_LS.GetXaxis().SetTitle(h[1])
        #hist_LS.GetYaxis().SetTitle(h[2])
        hist_D.Scale(1/hist_D.Integral())
        MaxD = hist_D.GetMaximum()
        g.Add(hist_D)
        

        if noStack:
            if (MaxMC > MaxD):
                g.SetMaximum(MaxMC*1.1)
            else:
                g.SetMaximum(MaxD*1.1)
            g.Draw("nostackhist")
        else:
            g.SetMaximum(g.GetMaximum()*1.1)
            g.Draw("hist")

        legende.Draw()
        g.GetXaxis().SetTitle(h[1])
        g.GetYaxis().SetTitle(h[2])
        g.GetYaxis().SetTitleOffset(1.4)
        g.GetXaxis().SetTitleOffset(1.1)
        my_Latex = TLatex()
        my_Latex.SetTextSize(0.04)
        my_Latex.DrawLatexNDC(0.13,0.85, "LHCb MC & Data 2012")
        CanvasName.Update()

        outDir = '/sps/lhcb/volle/Plot_Presel'
        if not os.path.isdir(outDir):
            os.makedirs(outDir)
        CanvasName.SaveAs("{}/{}_after.pdf".format(outDir,histoListe[i][0]))
        CanvasName.SaveAs("{}/{}_after.png".format(outDir,histoListe[i][0]))
        #CanvasName.SaveAs("{}/{}.root".format(outDir,histoListe[i][0]))


#__________________________________MAIN PART_________________________________


if __name__ == '__main__':
    GetAndDrawHistograms(TFileDataName,DataTreeName,TFileMCName,MCTreeName,noStack=True)
