#!/usr/bin/env python

'''
Project        : Plot Correlation of <BDTG> and Mass Lb
Script name    : PlotCorrBDTGBMass.py 
Author         : Felicia Volle
'''

# Python and ROOT import
import sys    # exit
import time   # time accounting
import getopt # command line parser
import os, pickle
import ROOT
from array import array
from ROOT import gStyle,gROOT,TLegend,TPad,TCanvas,TH1F,TH2F,TFile,TDirectory,TCut,THStack,TLatex,TTree


#________________________________SET VARIABLES________________________________

comment = "Secure"

#_____________________________DEFINITION OF FUNCTIONS________________________ 


def setupHistograms():
    histograms = [
        ["B_M","M(#Lambda_{b})","Entries / {MeV}","B_M>>h(100,3600,7100)"]
        ]
    return histograms
B_M_list = [5120,5170,5220,5270,5320,5370,5420,5470,5520,5570,5620,5670,5720,5770,5820,5870,5920,5970,6020,6070,6120]#[3720,3920,4120,4320,4520,4720,4920,5120,5320,5520,5720,5920,6120,6320,6520,6720,6920,7120,7320]

def DrawCorrelationWithBDTG(noStack=True):
    
    f = TFile("/sps/lhcb/volle/Data_BDTG_Secure.root")
    t = f.Get("DecayTree")
    fM = TFile("/sps/lhcb/volle/MCsignal_BDTG_Secure.root")#"all_BDTG.root")
    tM = fM.Get("DecayTree")

    x1 = []
    y1 = []
    ex1 = []
    ey1 = []

    x2 = []
    y2 = []
    ex2 = []
    ey2 = []

    outDir = '/sps/lhcb/volle/CheckBDTG'
    if not os.path.isdir(outDir):
        os.makedirs(outDir)


    for i in range(len(B_M_list)-1):
        #c1 = TCanvas("c1",'Defining histogram size') 
        condition = "B_M>{} && B_M<{}".format(B_M_list[i],B_M_list[i+1])
        Mass_val = (B_M_list[i]+B_M_list[i+1])/2
        #print Mass_val
        x1.append(Mass_val)
        x2.append(Mass_val)
        #fnew = ROOT.TFile("{}/Root_{}.root".format(outDir,i),"recreate")
        #tnew = t.CopyTree(condition)
        #fnew.Write()
        if B_M_list[i]>=6120 or B_M_list[i+1]<=5120:
            t.Draw("BDTG",condition)
            bdtg_mean = t.GetHistogram().GetMean()
            bdtg_mean_err = t.GetHistogram().GetMeanError()
            y1.append(bdtg_mean)
            #y2.append(1.2)
            #ey2.append(0)
            #ex2.append(0)
            ey1.append(bdtg_mean_err)
            ex1.append(0)
        else:
            tM.Draw("BDTG",condition)
            bdtg_mean = tM.GetHistogram().GetMean()
            bdtg_mean_err = tM.GetHistogram().GetMeanError()
            y2.append(bdtg_mean)
            ey2.append(bdtg_mean_err)
            ex2.append(0)
            #y1.append(1.2)
            #ey1.append(0)
            #ex1.append(0)
        

    legende=TLegend(0.7, 0.7, 0.89, 0.89)
    c2 = TCanvas("c2",'c2') 

    n = len(x1)
    xS = array('d',x1)
    yS = array('d',y1)
    exS = array('d',ex1)
    eyS = array('d',ey1)

    xM = array('d',x2)
    yM = array('d',y2)
    eyM = array('d',ey2)
    exM = array('d',ex2)
       
    #gS = ROOT.TGraphErrors(n,xS,yS,exS,eyS)
    #legende.AddEntry(gS, "Data", "lp")

    #gS.SetMarkerStyle(20) 
    #gS.SetMarkerColor(4)
    #gS.SetLineColor(4)

    gM = ROOT.TGraphErrors(n,xM,yM,exM,eyM)
    gM.SetMinimum(0.1)
    gM.SetMaximum(1.1)

    legende.AddEntry(gM, "MC", "lp")
    gM.SetMarkerStyle(20) 
    gM.SetMarkerColor(6)
    gM.SetLineColor(6)

    #gS.Draw()
    gM.Draw()

    gM.SetTitle("")
    gM.GetXaxis().SetTitle("M(#Lambda_{b})")
    gM.GetYaxis().SetTitle("<BDTG>")
    gM.GetYaxis().SetTitleOffset(1.4)
    gM.GetXaxis().SetTitleOffset(1.1)

    c2.Update()

    #legende.Draw()
    #my_Latex = TLatex()
    #my_Latex.SetTextSize(0.04)
    #my_Latex.DrawLatexNDC(0.13,0.85, "LHCb Data 2012")
    c2.Update()

    c2.SaveAs("{}/BDTG_{}.pdf".format(outDir,comment))
    c2.SaveAs("{}/BDTG_{}.png".format(outDir,comment))
    #c2.SaveAs("{}/BDTG_{}.root".format(outDir,comment))


#__________________________________MAIN PART_________________________________


if __name__ == '__main__':
    DrawCorrelationWithBDTG(noStack=True)
    
