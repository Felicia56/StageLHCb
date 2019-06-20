#!/usr/bin/env python

''' 
_____________________________________________________________________________

Project        : Plot the figure of Merit
Script name    : FigureOfMerit.py 
Author         : Felicia Volle
_____________________________________________________________________________

''' 

# Python and ROOT import
import sys    # exit
import time,datetime   # time accounting
import argparse
import os
import ROOT
from array import array
from simpleFit import simpleFit


#____________________________VARIABLE DEFINITION______________________________

f = ROOT.TFile("/sps/lhcb/volle/Data_RS_BDTG.root")
t = f.Get("DecayTree")
fSig = ROOT.TFile("/sps/lhcb/volle/MCsignal_2_BDTG.root")
tSig = fSig.Get("DecayTree")

BDTG_val_list = [x * 0.1 for x in range(-10,10)]#[-1,-0.8,-0.6,-0.4,-0.2,0,0.2,0.4,0.6,0.8]

n = len(BDTG_val_list)
x1 = []
y1 = []
ex1 = [x * 0 for x in range(-10,10)]
ey1 = []
x2 = []
y2 = []
ex2 = [x * 0 for x in range(-10,10)]
ey2 = []


#____________________________FUNCTION DEFINITION______________________________

def GetBackgroundYieldDict(t,version,xmin,xmax):
    Nbkg = {}
    #for i in BDTG_val_list:
        
    FileName = "{}/BDTG_txtfile.txt".format(newDir)
    if os.path.exists(FileName):
        ftxt = open(FileName,'r')
        Text = ftxt.read()
        print Textmmmm
        TextList = Text.split(";")
        print TextList
        for s in range(len(TextList)-1):
            txtlist =[]
            txtlist = TextList[s].split(" : ")
            Nbkg.update({txtlist[0] : txtlist[1]})
    else:
        ftxt = open(FileName,'w')
        for i in BDTG_val_list:
            CutVal = "BDTG>="+str(i)
            OutputList = []
            ds, model, chi2, Nbkg_SR, Nbkg_SR_error = simpleFit(t, version, CutVal, xmin, xmax)
            OutputList.append(Nbkg_SR)
            OutputList.append(Nbkg_SR_error)
            OutputList.append(chi2)
            ftxt.write('{} : {};'.format(i,OutputList))
            Nbkg.update({i : OutputList})
    ftxt.close()
    return Nbkg

def GetSignalEfficiencyAndYieldDict(fSig,tSig):
    NsigE = {} #for Efficiency calculation
    NsigY = {} #for the calculation of the expected events (yield)
    for j in BDTG_val_list:
        SignalYieldList = []
        CutVal = "BDTG>="+str(j) 
        CUT = CutVal + " && B_M>5120 && B_M<6120"

       # FileName = "{}/BDTG_SIGfile_{}.root".format(newDir,j)
       # if os.path.exists(FileName):
       #     fnew = ROOT.TFile(FileName)
       #     tnew = fnew.Get("DecayTree")
       # else:
       #     fnew = ROOT.TFile(FileName,"recreate")
       #     tnew = tSig.CopyTree(CUT)
       #     fnew.Write()
        
        NsigMC = tSig.GetEntries(CUT)
        NsigE.update({j : NsigMC})
        
        NsigYield = NsigMC * 3*10**(15) * 568*10**(-6) * 2 * 0.17 * 3.39*10**(-5)
        NsigYield_error = NsigYield * ROOT.TMath.Sqrt( ((0.48*10**(-5))/(3.39*10**(-5)))**2 )
        SignalYieldList.append(NsigYield)
        SignalYieldList.append(NsigYield_error)
        NsigY.update({j : SignalYieldList})
    return NsigE, NsigY

def DrawPunzi(Nbkg,NsigE):
    b = 0
    for b in BDTG_val_list:
        x1.append(b)
        eff_sig = NsigE[b]/NsigE[-1]
        print 'eff_sig = {}'.format(eff_sig)
        y_val = eff_sig/( 2.5 + ROOT.TMath.Sqrt(Nbkg[b][0]) )
        print 'y_val = {}'.format(y_val)
        dp_sur_db = eff_sig/( ((2.5 + ROOT.TMath.Sqrt(Nbkg[b][0]))**2) * (2*ROOT.TMath.Sqrt(Nbkg[b][0])) )
        y_val_err = dp_sur_db * Nbkg[b][1]
        print 'y_val_err = {}'.format(y_val_err)
        y1.append(y_val)
        ey1.append(y_val_err)
        b = b+1
       
    xP = array('f',x1)
    yP = array('f',y1)
    exP = array('f',ex1)
    eyP = array('f',ey1)

    cP = ROOT.TCanvas('cP','cP')
    ROOT.gROOT.SetStyle('Plain')
    #ROOT.gROOT.SetOptTitle(0)
    gP = ROOT.TGraphErrors(n,xP,yP,exP,eyP)
    gP.SetMarkerStyle(20)
    gP.SetMarkerColor(4)
    gP.SetLineColor(4)
    gP.Draw()
    gP.GetXaxis().SetTitle("BDTG cut value")
    gP.GetYaxis().SetTitle("f_{PUNZI}")
    cP.Update()
    cP.SaveAs('{}/FigureOfMerit_punzi.eps'.format(newDir))
    #return xP, yP, exP, eyP

def DrawSignificance(Nbkg,NsigY):
    b = 0
    for b in BDTG_val_list:
        x2.append(b)
        print 'b = {}'.format(b)
        y_val2 = NsigY[b][0]/(NsigY[b][0] + Nbkg[b][0])
        print 'y_val2 = {}'.format(y_val2)
        df_sur_ds = Nbkg[b][0]/( (NsigY[b][0] + Nbkg[b][0])**2 )
        df_sur_db = NsigY[b][0]/((NsigY[b][0] + Nbkg[b][0])**2)
        y_val_err2 = ROOT.TMath.Sqrt((df_sur_ds*NsigY[b][1])**2 + (df_sur_db*Nbkg[b][1])**2)
        print 'y_val_err2 = {}'.format(y_val_err2)
        y2.append(y_val2)
        ey2.append(y_val_err2)
        b = b+1
       
    xS = array('f',x2)
    yS = array('f',y2)
    exS = array('f',ex2)
    eyS = array('f',ey2)

    cS = ROOT.TCanvas('cS','cS')
    ROOT.gROOT.SetStyle('Plain')
    #ROOT.gROOT.SetOptTitle(0)
    gS = ROOT.TGraphErrors(n,xS,yS,exS,eyS)
    gS.SetMarkerStyle(20)
    gS.SetMarkerColor(4)
    gS.SetLineColor(4)
    gS.Draw()
    gS.GetXaxis().SetTitle("BDTG cut value")
    gS.GetYaxis().SetTitle("s/(s+b)")
    cS.Update()
    cS.SaveAs('{}/FigureOfMerit_signif.eps'.format(newDir))
    #return xS, yS, exS, eyS




#____________________________________MAIN_____________________________________


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    #parser.add_argument("file", action="store", type=str)
    #parser.add_argument("-t", "--tree", default="DecayTree", action="store", type=str)
    #parser.add_argument("-m", "--mean", default=5620., action="store", type=float)
    parser.add_argument("-n", "--xmin", default=6120., action="store", type=float)
    parser.add_argument("-x", "--xmax", default=9120., action="store", type=float)
    #parser.add_argument("-c", "--cuts", default="", action="store", type=str)
    parser.add_argument("-v", "--version", default="v1",action="store", type=str)
    args = parser.parse_args()

    newDir = '/sps/lhcb/volle/FitResults/{}'.format(args.version)
    if not os.path.isdir(newDir):
        os.makedirs(newDir)

    #Save for all BDTG_cut_values the Number of bkg events, error & chi2 of fit
    Nbkg = GetBackgroundYieldDict(t,args.version,args.xmin,args.xmax) 
    print "Nbkg = "
    print Nbkg

    #Save for all BDTG_cut_values the Number of Signal events & their error
    NsigE, NsigY =  GetSignalEfficiencyAndYieldDict(fSig,tSig)
    print "NsigE = "
    print NsigE
    print "NsigY = "
    print NsigY
    
    #Draw Punzi Figure of merit
    DrawPunzi(Nbkg,NsigE)
    
    #Draw Significance Figure of merit
    DrawSignificance(Nbkg,NsigY)



  
        
        
        


    


