#!/usr/bin/env python

'''
Project         : perform a splot
Script name     : SPlot.py
Author          : Felicia Volle
'''

# Python and ROOT import
import sys    # exit
import time   # time accounting
import getopt # command line parser
import os, pickle
from ROOT import gStyle,gROOT,TLegend,TPad,TCanvas,TH1F,TH2F,TFile,TDirectory,TCut,THStack,TLatex,TTree


#________________________________SET VARIABLES________________________________


TFileSweightName = "/users/LHCb/volle/Stage/sData.root"
TFileDataName = "/sps/lhcb/volle/FitBMass/v1_Try2/BDTG_BKGfile_0.92.root"#v2/SPlot_output.root"#TrigSel.root"
SweightTreeName = "sTree"
DataTreeName = "DecayTree"


#_____________________________DEFINITION OF FUNCTIONS________________________ 


def setupHistograms():
    histograms = [
        ["B_M","M(#Lambda_{b})","Entries / {MeV}","B_M>>h(100,3600,7100)"],
        ["Lambda_1520_0_M","M(#Lambda^{*})","Entries / {MeV}","Lambda_1520_0_M>>h(100,1400,2600)"],
        ["Lambda_1520_0_cosThetaH","cos(#theta_{H,#Lambda^{*}})","Entries","Lambda_1520_0_cosThetaH>>h(100,-1,1)"]
        ]
    return histograms

def DrawSplotHistograms(DataFileName,DataTreeName,SFileName,STreeName,noStack=True):
    #fS = TFile(SFileName)
    #tS = fS.Get(STreeName)
    fData = TFile(DataFileName)
    tData = fData.Get(DataTreeName)

    tData.AddFriend(STreeName,SFileName)
    histoListe=setupHistograms()
    HistDic = {}

    for i,h in enumerate(histoListe):
        legende=TLegend(0.75, 0.83, 0.97, 0.97)
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
        
        tData.Draw(histoListe[i][3])
        hist_data = tData.GetHistogram().Clone("hist_data")
        hist_data.SetLineColor(5)
        hist_data.SetLineWidth(1)
        #hist_data.SetMarkerStyle(21)
        legende.AddEntry(hist_data, "Data", "lp")
        hist_data.GetXaxis().SetTitle(h[1])
        hist_data.GetYaxis().SetTitle(h[2])
        MaxData = hist_data.GetMaximum()
        g.Add(hist_data)
        
        tData.Draw(histoListe[i][3],"nsig_sw")
        hist_sw = tData.GetHistogram().Clone("hist_sw")
        hist_sw.SetLineColor(4)
        hist_sw.SetLineWidth(1)
        ##hist_RS.SetMarkerStyle(21)
        legende.AddEntry(hist_sw, "sweighted Data", "lp")
        MaxSW = hist_sw.GetMaximum()
        g.Add(hist_sw)

        if noStack:
            if (MaxSW > MaxData):
                g.SetMaximum(MaxSW*1.2)
            else:
                g.SetMaximum(MaxData*1.2)
            g.Draw("nostackhist")
        else:
            g.SetMaximum(g.GetMaximum()*1.2)
            g.Draw("hist")

        legende.Draw()
        g.GetXaxis().SetTitle(h[1])
        g.GetYaxis().SetTitle(h[2])
        g.GetYaxis().SetTitleOffset(1.4)
        g.GetXaxis().SetTitleOffset(1.1)
        my_Latex = TLatex()
        my_Latex.SetTextSize(0.04)
        my_Latex.DrawLatexNDC(0.13,0.85, "LHCb Data 2012")
        CanvasName.Update()

        outDir = '/sps/lhcb/volle/FitBMass/v1_Try2'
        if not os.path.isdir(outDir):
            os.makedirs(outDir)
        CanvasName.SaveAs("{}/{}_both.pdf".format(outDir,histoListe[i][0]))#_both
        CanvasName.SaveAs("{}/{}_both.png".format(outDir,histoListe[i][0]))#_both
        #CanvasName.SaveAs("{}/{}.root".format(outDir,histoListe[i][0]))


#__________________________________MAIN PART_________________________________


if __name__ == '__main__':
    DrawSplotHistograms(TFileDataName,DataTreeName,TFileSweightName,SweightTreeName,noStack=True)
    
