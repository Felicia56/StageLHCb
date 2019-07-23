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


TFileMCName = "/sps/lhcb/volle/MCsignal_TrigSel.root"#"/data/lhcb/marin/lb2pkgamma/MC/2012/15102203/2hG-S21/radiative2hG_MC2012-Lb2L1520gamma_HighPt-15102203-Py8Sim09dReco14c_S21.root"#
TFileDataName = "/sps/lhcb/volle/Data_all.root"#"/data/lhcb/marin/lb2pkgamma/Data/2012/2hG-S21/radiative2hG_R14S21_MagUp_tmp.root"#"/sps/lhcb/volle/Data_full_TrigSel.root"
MCTreeName = "DecayTree"#"pkGTupleMC/DecayTree"#
DataTreeName = "DecayTree"#"pkGTuple/DecayTree"#"DecayTree"#


CutMC = TCut("pplus_ProbNNp>0.2 && Kminus_ProbNNk>0.2 && B_PT>4000 && Lambda_1520_0_PT>1500 && gamma_PT>3000 && pplus_PT>1000 && B_FDCHI2_OWNPV>100 && pplus_IPCHI2_OWNPV>50 && Kminus_IPCHI2_OWNPV>40 && B_BKGCAT==0")

CutData_LS = TCut("pplus_ProbNNp>0.2 && Kminus_ProbNNk>0.2 && B_PT>4000 && Lambda_1520_0_PT>1500 && gamma_PT>3000 && pplus_PT>1000 && B_FDCHI2_OWNPV>100 && pplus_IPCHI2_OWNPV>50 && Kminus_IPCHI2_OWNPV>40 && B_M<5120")

CutData_RS = TCut("pplus_ProbNNp>0.2 && Kminus_ProbNNk>0.2 && B_PT>4000 && Lambda_1520_0_PT>1500 && gamma_PT>3000 && pplus_PT>1000 && B_FDCHI2_OWNPV>100 && pplus_IPCHI2_OWNPV>50 && Kminus_IPCHI2_OWNPV>40 && B_M>6120")


#_____________________________DEFINITION OF FUNCTIONS________________________ 


def setupHistograms():
    histograms = [
        ["pplus_PT","P_{T}(p)","Entries / {MeV}","pplus_PT>>h(100,900,5600)"],
        ["Kminus_PT","P_{T}(K^{-})","Entries / {MeV}","Kminus_PT>>h(100,0,5000)"],
        ["gamma_PT","P_{T}(#gamma)","Entries / {MeV}","gamma_PT>>h(100,2800,11300)"],
        ["Lambda_1520_0_PT","P_{T}(#Lambda*)","Entries / {MeV}","Lambda_1520_0_PT>>h(100,1000,10000)"],
        ["B_PT","P_{T}(#Lambda_{b})","Entries / {MeV}","B_PT>>h(100,3500,20000)"],
        ["B_M","M(#Lambda_{b})","Entries / {MeV}","B_M>>h(100,3600,7100)"],
        ["beta","#beta","Entries / {}","(-(gamma_P-Kminus_P-pplus_P)/(gamma_P+Kminus_P+pplus_P))>>h(100,-0.95,0.95)"],
        ["MomCons1","P_{tot,1}","Entries / {MeV}","-(B_P-gamma_P-Lambda_1520_0_P)>>h(100,-10,400)"],
        #["MomCons1X","P_{x,1}","Entries / {MeV}","-(B_PX-gamma_PX-Lambda_1520_0_PX)>>h(100,-100,100)"],
        #["MomCons1Y","P_{y,1}","Entries / {MeV}","-(B_PY-gamma_PY-Lambda_1520_0_PY)>>h(100,-100,100)"],
        #["MomCons1Z","P_{z,1}","Entries / {MeV}","-(B_PZ-gamma_PZ-Lambda_1520_0_PZ)>>h(100,-15,15)"],
        #["MomCons1E","E_{tot,1}","Entries / {MeV}","-(B_PE-gamma_PE-Lambda_1520_0_PE)>>h(100,-1,1)"],
        ["MomCons2","P_{tot,2}","Entries / {MeV}","-(Lambda_1520_0_P-Kminus_P-pplus_P)>>h(100,-10,40)"],
        ["Sum_Kminus_p_eta","#eta(K^{-})+#eta(p^{+})","Entries / {}","(atanh(pplus_PZ/pplus_P)+atanh(Kminus_PZ/Kminus_P))>>h(100,3.5,10)"],
        ["Diff_Kminus_p_eta","#eta(K^{-})-#eta(p^{+})","Entries / {}","(atanh(Kminus_PZ/Kminus_P)-atanh(pplus_PZ/pplus_P))>>h(100,-1,1)"],
        ["Kminus_eta","#eta(K^{-})","Entries / {}","(atanh(Kminus_PZ/Kminus_P))>>h(100,1.5,5)"],
        ["pplus_eta","#eta(p^{+})","Entries / {}","(atanh(pplus_PZ/pplus_P))>>h(100,1.5,5)"],
        ["gamma_eta","#eta(#gamma)","Entries / {}","(atanh(gamma_PZ/gamma_P))>>h(100,1.5,5.5)"],
        ["Lambda_1520_0_eta","#eta(#Lambda*)","Entries / {}","(atanh(Lambda_1520_0_PZ/Lambda_1520_0_P))>>h(100,1.5,5)"],
        ["B_FDCHI2_OWNPV","#chi^{2}_{FD}(#Lambda_{b})","Entries / {}","B_FDCHI2_OWNPV>>h(100,0,2000)"],
        ["Lambda_1520_0_FDCHI2_OWNPV","#chi^{2}_{FD}(#Lambda*)","Entries / {}","Lambda_1520_0_FDCHI2_OWNPV>>h(100,0,2000)"],
        ["pplus_IPCHI2_OWNPV","#chi^{2}_{IP}(p^{+})","Entries / {}","pplus_IPCHI2_OWNPV>>h(100,0,1000)"],
        ["Kminus_IPCHI2_OWNPV","#chi^{2}_{IP}(K^{-})","Entries / {}","Kminus_IPCHI2_OWNPV>>h(100,0,700)"],
        ["Lambda_1520_0_IPCHI2_OWNPV","#chi^{2}_{IP}(#Lambda*)","Entries / {}","Lambda_1520_0_IPCHI2_OWNPV>>h(100,0,1000)"],
        ["B_IPCHI2_OWNPV","#chi^{2}_{IP}(#Lambda_{b})","Entries / {}","B_IPCHI2_OWNPV>>h(100,0,10)"],
        ["Lambda_1520_0_OWNPV_CHI2","#chi^{2}_{vertex}(#Lambda*)","Entries / {}","Lambda_1520_0_OWNPV_CHI2>>h(100,0,90)"],
        ["B_OWNPV_CHI2","#chi^{2}_{vertex}(#Lambda_{b})","Entries / {}","B_OWNPV_CHI2>>h(100,0,90)"],
        ["B_BMassFit_chi2_per_nDOF","#chi^{2}_{DTF}/n_{dof}","Entries / {}","B_BMassFit_chi2/B_BMassFit_nDOF>>h(100,0,500)"],
        ["B_BMassFit_chi2_per_nDOF_zoom","#chi^{2}_{DTF,BM}/n_{dof}","Entries / {}","B_BMassFit_chi2/B_BMassFit_nDOF>>h(100,0,20)"],
        ["B_PVFit_chi2_per_nDOF","#chi^{2}_{DTF,PV}/n_{dof}","Entries / {}","B_PVFit_chi2/B_PVFit_nDOF>>h(100,0,5)"],
        ["Lambda_1520_0_DIRA_OWNPV","DIRA(#Lambda*)","Entries / {}","Lambda_1520_0_DIRA_OWNPV>>h(100,0.99,1.004)"],
        ["B_DIRA_OWNPV","DIRA(#Lambda_{b})","Entries / {}","B_DIRA_OWNPV>>h(100,0.9998,1.00007)"],
        ["B_PVFit_nPV","nPV(#Lambda_{b},PVFit)","Entries / {}","B_PVFit_nPV>>h(100,0,9)"], 
        ["Angle_p_K","#Theta(p,K^{-})-#pi","Entries / {rad}","acos(pplus_CosTheta)+acos(Kminus_CosTheta)-TMath::Pi()>>h(100,-0.005,0.005)"] 
        

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
        hist_MC.SetLineColor(6)
        hist_MC.SetLineWidth(1)
        #hist_MC.SetMarkerStyle(21)
        legende.AddEntry(hist_MC, "MC", "lp")
        hist_MC.GetXaxis().SetTitle(h[1])
        hist_MC.GetYaxis().SetTitle(h[2])
        hist_MC.Scale(1/hist_MC.Integral())
        MaxMC = hist_MC.GetMaximum()
        g.Add(hist_MC)
        
        tData.Draw(histoListe[i][3],CutData_LS)
        hist_LS = tData.GetHistogram().Clone("hist_LS")
        hist_LS.SetLineColor(3)
        hist_LS.SetLineWidth(1)
        #hist_LS.SetMarkerStyle(21)
        legende.AddEntry(hist_LS, "Data LS", "lp")
        #hist_LS.GetXaxis().SetTitle(h[1])
        #hist_LS.GetYaxis().SetTitle(h[2])
        hist_LS.Scale(1/hist_LS.Integral())
        MaxLS = hist_LS.GetMaximum()
        g.Add(hist_LS)
        
        tData.Draw(histoListe[i][3],CutData_RS)
        hist_RS = tData.GetHistogram().Clone("hist_RS")
        hist_RS.SetLineColor(4)
        hist_RS.SetLineWidth(1)
        #hist_RS.SetMarkerStyle(21)
        legende.AddEntry(hist_RS, "Data RS", "lp")
        hist_RS.Scale(1/hist_RS.Integral())
        MaxRS = hist_RS.GetMaximum()
        g.Add(hist_RS)

        if noStack:
            if (MaxMC > MaxLS and MaxMC > MaxRS):
                g.SetMaximum(MaxMC*1.1)
            elif (MaxLS > MaxMC and MaxLS > MaxRS):
                g.SetMaximum(MaxLS*1.1)
            else:
                g.SetMaximum(MaxRS*1.1)
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

        outDir = '/sps/lhcb/volle/Plots_RS_LS'
        if not os.path.isdir(outDir):
            os.makedirs(outDir)
        CanvasName.SaveAs("{}/{}.pdf".format(outDir,histoListe[i][0]))
        CanvasName.SaveAs("{}/{}.png".format(outDir,histoListe[i][0]))
        #CanvasName.SaveAs("{}/{}.root".format(outDir,histoListe[i][0]))


#__________________________________MAIN PART_________________________________


if __name__ == '__main__':
    GetAndDrawHistograms(TFileDataName,DataTreeName,TFileMCName,MCTreeName,noStack=True)
    
    
    
