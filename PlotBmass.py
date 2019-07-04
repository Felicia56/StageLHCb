#!/usr/bin/env python

'''
Project        : Plot Bmass
Script name    : PlotBmass.py 
Author         : Felicia Volle
'''

# Python and ROOT import
import sys    # exit
import time   # time accounting
import getopt # command line parser
import os, pickle
from ROOT import gStyle,gROOT,TLegend,TPad,TCanvas,TH1F,TH2F,TFile,TDirectory,TCut,THStack,TLatex,TTree


#________________________________SET VARIABLES________________________________


TFileMCName = "/sps/lhcb/volle/MCsignal_TrigSel.root"
TFileDataName = "/sps/lhcb/volle/Data_all.root"#TrigSel.root"
MCTreeName = "DecayTree"
DataTreeName = "DecayTree"

comment = "all"

CutMC = TCut("")#"pplus_ProbNNp>0.2 && Kminus_ProbNNk>0.2 && B_PT>4000 && Lambda_1520_0_PT>1500 && gamma_PT>3000 && pplus_PT>1000 && B_FDCHI2_OWNPV>100 && pplus_IPCHI2_OWNPV>50 && Kminus_IPCHI2_OWNPV>40 && B_BKGCAT==0")

CutData = TCut("")#"pplus_ProbNNp>0.2 && Kminus_ProbNNk>0.2 && B_PT>4000 && Lambda_1520_0_PT>1500 && gamma_PT>3000 && pplus_PT>1000 && B_FDCHI2_OWNPV>100 && pplus_IPCHI2_OWNPV>50 && Kminus_IPCHI2_OWNPV>40 && B_M<5120")


#_____________________________DEFINITION OF FUNCTIONS________________________ 


def setupHistograms():
    histograms = [
        ["B_M","M(#Lambda_{b})","Entries / {MeV}","B_M>>h(100,3600,7100)"]
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
        
        #tMC.Draw(histoListe[i][3],CutMC)
        #hist_MC = tMC.GetHistogram().Clone("hist_MC")
        #hist_MC.SetLineColor(6)
        #hist_MC.SetLineWidth(1)
        ##hist_MC.SetMarkerStyle(21)
        #legende.AddEntry(hist_MC, "MC", "lp")
        #hist_MC.GetXaxis().SetTitle(h[1])
        #hist_MC.GetYaxis().SetTitle(h[2])
        #hist_MC.Scale(1/hist_MC.Integral())
        #MaxMC = hist_MC.GetMaximum()
        #g.Add(hist_MC)
        
        tData.Draw(histoListe[i][3],CutData)
        hist_data = tData.GetHistogram().Clone("hist_data")
        hist_data.SetLineColor(9)
        hist_data.SetLineWidth(1)
        #hist_data.SetMarkerStyle(21)
        legende.AddEntry(hist_data, "Data", "lp")
        hist_data.GetXaxis().SetTitle(h[1])
        hist_data.GetYaxis().SetTitle(h[2])
        #hist_data.Scale(1/hist_data.Integral())
        MaxData = hist_data.GetMaximum()
        g.Add(hist_data)
        
        #tData.Draw(histoListe[i][3],CutData_RS)
        #hist_RS = tData.GetHistogram().Clone("hist_RS")
        #hist_RS.SetLineColor(4)
        #hist_RS.SetLineWidth(1)
        ##hist_RS.SetMarkerStyle(21)
        #legende.AddEntry(hist_RS, "Data RS", "lp")
        #hist_RS.Scale(1/hist_RS.Integral())
        #MaxRS = hist_RS.GetMaximum()
        #g.Add(hist_RS)

        #if noStack:
        #    if (MaxMC > MaxData):
        #        g.SetMaximum(MaxMC*1.1)
        #    else:
        #        g.SetMaximum(MaxData*1.1)
        #    g.Draw("nostackhist")
        #else:
        #    g.SetMaximum(g.GetMaximum()*1.1)
        #    g.Draw("hist")
        g.SetMaximum(MaxData*1.1)
        g.Draw("hist")
        #legende.Draw()
        g.GetXaxis().SetTitle(h[1])
        g.GetYaxis().SetTitle(h[2])
        g.GetYaxis().SetTitleOffset(1.4)
        g.GetXaxis().SetTitleOffset(1.1)
        my_Latex = TLatex()
        my_Latex.SetTextSize(0.04)
        my_Latex.DrawLatexNDC(0.13,0.85, "LHCb Data 2012")
        CanvasName.Update()

        outDir = '/sps/lhcb/volle/BMass'
        if not os.path.isdir(outDir):
            os.makedirs(outDir)
        CanvasName.SaveAs("{}/{}_{}.pdf".format(outDir,histoListe[i][0],comment))
        CanvasName.SaveAs("{}/{}_{}.png".format(outDir,histoListe[i][0],comment))
        #CanvasName.SaveAs("{}/{}_{}.root".format(outDir,histoListe[i][0],comment))


#__________________________________MAIN PART_________________________________


if __name__ == '__main__':
    GetAndDrawHistograms(TFileDataName,DataTreeName,TFileMCName,MCTreeName,noStack=True)
    
    
    
