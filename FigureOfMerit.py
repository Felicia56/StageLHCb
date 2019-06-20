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

BDTG_val_list = [-1,-0.9,-0.8,-0.7,-0.6,-0.5,-0.4,-0.3,-0.2,-0.1,0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9]#[x * 0.1 for x in range(-10,10)]

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
        print Text
        TextList = Text.split("];")
        print TextList
        for s in range(len(TextList)-1):
            txtlist = []
            txtlist2 = []
            txtlist = TextList[s].split(" : [")
            txtlist2 = txtlist[1].split(", ")
            txtlist3 = [float(txtlist2[0]),float(txtlist2[1]),float(txtlist2[2])]
            Nbkg.update({float(txtlist[0]) : txtlist3})
    else:
        ftxt = open(FileName,'w')
        for i in BDTG_val_list:
            CutVal = "BDTG>="+str(i)
            OutputList = []
            ds, model, chi2, Nbkg_SR, Nbkg_SR_error = simpleFit(t, version, CutVal, xmin, xmax)
            OutputList.append(float(Nbkg_SR))
            OutputList.append(float(Nbkg_SR_error))
            OutputList.append(float(chi2))
            ftxt.write('{} : {};'.format(i,OutputList))
            Nbkg.update({i : OutputList})
    ftxt.close()
    return Nbkg

def GetSignalEfficiencyAndYieldDict(fSig,tSig):
    NsigE = {} #for Efficiency calculation
    NsigY = {} #for the calculation of the expected events (yield)
    Sum_Nsig = 0
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
        Sum_Nsig = Sum_Nsig + NsigMC

    Var_Nsig = Sum_Nsig/len(BDTG_val_list) #is the error of all MC event numbers
        
    for j in BDTG_val_list:
        SignalYieldList = []
        CutVal = "BDTG>="+str(j) 
        CUT = CutVal + " && B_M>5120 && B_M<6120"

        NsigMC = tSig.GetEntries(CUT)
        NsigE.update({j : NsigMC})
        
        #Lumi=3*10**(15);Crosssection-bb=568*10**(-6)*2;f_Lb=0.17;BR=3.39*10**(-5)
        NsigYield = float(NsigMC) * 3*10**(15) * 568*10**(-6) * 2 * 0.17 * 3.39*10**(-5)
        NsigYield_error = NsigYield * ROOT.TMath.Sqrt( (0.48/3.39)**2 + (float(Var_Nsig)/NsigMC)**2)
        SignalYieldList.append(NsigYield)
        SignalYieldList.append(NsigYield_error)
        NsigY.update({j : SignalYieldList})

    return NsigE, NsigY, Var_Nsig

def DrawPunzi(Nbkg, NsigE, Var_Nsig):
    b = 0
    for b in BDTG_val_list:
        x1.append(b)
        eff_sig = float(NsigE[b])/NsigE[-1]
        print 'eff_sig = {}'.format(eff_sig)
        y_val = eff_sig/( 2.5 + ROOT.TMath.Sqrt(Nbkg[b][0]) )
        print 'y_val = {}'.format(y_val)
        dp_sur_db = eff_sig/( ((2.5 + ROOT.TMath.Sqrt(Nbkg[b][0]))**2) * (2*ROOT.TMath.Sqrt(Nbkg[b][0])) )
        dp_sur_deff = 1/(2.5 + ROOT.TMath.Sqrt(Nbkg[b][0]))
        y_val_err = y_val*ROOT.TMath.Sqrt( (dp_sur_db * Nbkg[b][1])**2 + (dp_sur_deff * Var_Nsig/NsigE[-1])**2 )
        print 'y_val_err = {}'.format(y_val_err)
        y1.append(y_val)
        ey1.append(y_val_err)
        b = b+1
       
    xP = array('d',x1)
    yP = array('d',y1)
    exP = array('d',ex1)
    eyP = array('d',ey1)

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

def DrawSignificance(Nbkg, NsigY, Var_Nsig):
    b = 0
    for b in BDTG_val_list:
        x2.append(b)
        print 'b = {}'.format(b)
        y_val2 = float(NsigY[b][0])/ ROOT.TMath.Sqrt( float(NsigY[b][0]) + Nbkg[b][0] )
        print 'y_val2 = {}'.format(y_val2)
        df_sur_ds = (float(NsigY[b][0]) + 2*Nbkg[b][0])/(2* (ROOT.TMath.Sqrt(float(NsigY[b][0]) + Nbkg[b][0]))**3 )
        df_sur_db = (float(NsigY[b][0]))/(2* (ROOT.TMath.Sqrt(float(NsigY[b][0]) + Nbkg[b][0]))**3 )
        y_val_err2 = ROOT.TMath.Sqrt((df_sur_ds*float(NsigY[b][1]))**2 + (df_sur_db*Nbkg[b][1])**2)
        print 'y_val_err2 = {}'.format(y_val_err2)
        print 'Error on signal events = {}'.format(df_sur_ds*float(NsigY[b][1]))
        print 'Error on bkg events = {}'.format(df_sur_db*Nbkg[b][1])
        y2.append(y_val2)
        ey2.append(y_val_err2)
        b = b+1
       
    xS = array('d',x2)
    yS = array('d',y2)
    exS = array('d',ex2)
    eyS = array('d',ey2)

    cS = ROOT.TCanvas('cS','cS')
    ROOT.gROOT.SetStyle('Plain')
    #ROOT.gROOT.SetOptTitle(0)
    gS = ROOT.TGraphErrors(n,xS,yS,exS,eyS)
    gS.SetMarkerStyle(20)
    gS.SetMarkerColor(4)
    gS.SetLineColor(4)
    gS.Draw()
    gS.GetXaxis().SetTitle("BDTG cut value")
    gS.GetYaxis().SetTitle("s/srt(s+b)")
    cS.Update()
    cS.SaveAs('{}/FigureOfMerit_signif.eps'.format(newDir))


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
    #print "Nbkg = "
    #print Nbkg

    #Save for all BDTG_cut_values the Number of Signal events & their error
    NsigE, NsigY, Var_Nsig =  GetSignalEfficiencyAndYieldDict(fSig,tSig)
    #print "NsigE = "
    #print NsigE
    #print "NsigY = "
    #print NsigY
    
    #Draw Punzi Figure of merit
    DrawPunzi(Nbkg, NsigE, Var_Nsig)
    
    #Draw Significance Figure of merit
    DrawSignificance(Nbkg, NsigY, Var_Nsig)



  
        
        
        


    


