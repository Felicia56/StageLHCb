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
from ROOT import TLatex
from array import array
from simpleFit import simpleFit


#____________________________VARIABLE DEFINITION______________________________

f = ROOT.TFile("/sps/lhcb/volle/Data_all_BDTG.root")#"full_BDTG.root")
t = f.Get("DecayTree")
fSig = ROOT.TFile("/sps/lhcb/volle/MCsignal_all_BDTG.root")#_2_BDTG.root")
tSig = fSig.Get("DecayTree")

#ftot = ROOT.TFile("/data/lhcb/marin/lb2pkgamma/MC/2012/15102203/2hG-S21/radiative2hG_MC2012-Lb2L1520gamma_HighPt-15102203-Py8Sim09dReco14c_S21.root")
#ttot = ftot.Get("pkGTupleMC/DecayTree")
#Get the total Number of MC events before the preselection
#Ntot = ttot.GetEntries("B_BKGCAT==0")#before selection
Ntot = (508316+505855)/0.2329
print "Ntot = {}".format(Ntot)

BDTG_val_list2 = [ 0.5, 0.52, 0.54, 0.56, 0.58, 0.6, 0.62, 0.64, 0.66, 0.68, 0.7, 0.72, 0.74, 0.76, 0.78, 0.8, 0.82, 0.84, 0.86, 0.88, 0.90, 0.92, 0.94, 0.96,0.98]#vZOOM[-1,-0.9,-0.8,-0.7,-0.6,-0.5,-0.4,-0.3,-0.2,-0.1,0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9]#[0, 0.02, 0.04, 0.06, 0.08, 0.1, 0.12, 0.14, 0.16, 0.18, 0.2, 0.22, 0.24, 0.26, 0.28, 0.3, 0.32, 0.34, 0.36, 0.38, 0.4, 0.42, 0.44, 0.46, 0.48, 0.5, 0.52, 0.54, 0.56, 0.58, 0.6, 0.62, 0.64, 0.66, 0.68, 0.7, 0.72, 0.74, 0.76, 0.78, 0.8, 0.82, 0.84, 0.86, 0.88, 0.90, 0.92, 0.94, 0.96,0.98]#v4
#[0, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.90, 0.95]#v4
BDTG_val_list = [-1,-0.9,-0.8,-0.7,-0.6,-0.5,-0.4,-0.3,-0.2,-0.1,0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9]
#[x * 0.1 for x in range(-10,10)]

n = len(BDTG_val_list2)
x1 = []
y1 = []
ex1 = [x * 0 for x in range(0,len(BDTG_val_list2))]
ey1 = []
x2 = []
y2 = []
ex2 = [x * 0 for x in range(0,len(BDTG_val_list2))]
ey2 = []
xb = []
yb = []
exb = [x * 0 for x in range(0,len(BDTG_val_list2))]
eyb = []
xs = []
ys = []
exs = [x * 0 for x in range(0,len(BDTG_val_list2))]
eys = []


#____________________________FUNCTION DEFINITION______________________________


def GetBackgroundYieldDict(t,version,xmini,xmin,xmid,xmax):
    Nbkg = {}
    #for i in BDTG_val_list:
        
    FileName = "{}/BDTG_txtfile.txt".format(newDir)
    if os.path.exists(FileName):
        ftxt = open(FileName,'r')
        Text = ftxt.read()
        #print Text
        TextList = Text.split("];")
        #print TextList
        for s in range(len(TextList)-1):
            txtlist = []
            txtlist2 = []
            txtlist = TextList[s].split(" : [")
            txtlist2 = txtlist[1].split(", ")
            txtlist3 = [float(txtlist2[0]),float(txtlist2[1]),float(txtlist2[2])]
            Nbkg.update({float(txtlist[0]) : txtlist3})
    else:
        ftxt = open(FileName,'w')
        for i in BDTG_val_list2:
            CutVal = "BDTG>="+str(i)
            OutputList = []
            ds, model, chi2, Nbkg_SR, Nbkg_SR_error = simpleFit(t, version, CutVal, xmini, xmin, xmid, xmax)
            OutputList.append(float(Nbkg_SR))
            OutputList.append(float(Nbkg_SR_error))
            OutputList.append(float(chi2))
            ftxt.write('{} : {};'.format(i,OutputList))
            Nbkg.update({i : OutputList})
    ftxt.close()
    return Nbkg

def GetSignalEfficiencyAndYieldDict(fSig,tSig,xmin,xmid):
    NsigE = {} #for Efficiency calculation
    NsigY = {} #for the calculation of the expected events (yield)
    #Sum_Nsig = 0
    #Exp_Nsig = 0

    CUT0 = "B_M>{} && B_M<{}".format(xmin,xmid)
    #NsigMC_CUT0 = tSig.GetEntries(CUT0)
    #print "NsigMC_CUT0 = {}".format(NsigMC_CUT0)
    #NsigE.update({-1 : NsigMC_CUT0})        

    #for j in BDTG_val_list:
    #    CutVal = "BDTG>="+str(j) 
    #    CUT = CutVal + " && "+ CUT0

       # FileName = "{}/BDTG_SIGfile_{}.root".format(newDir,j)
       # if os.path.exists(FileName):
       #     fnew = ROOT.TFile(FileName)
       #     tnew = fnew.Get("DecayTree")
       # else:
       #     fnew = ROOT.TFile(FileName,"recreate")
       #     tnew = tSig.CopyTree(CUT)
       #     fnew.Write()
        
    #    NsigMC = tSig.GetEntries(CUT)
    #    Sum_Nsig = Sum_Nsig + NsigMC
    #    Exp_Nsig = Exp_Nsig + NsigMC**2
    
    #Var_Nsig = Sum_Nsig/len(BDTG_val_list)#Exp_Nsig/Sum_Nsig #is the error of all MC event numbers
    #print "Var_Nsig = {}".format(Var_Nsig)

    for j in BDTG_val_list2:
        SignalYieldList = []
        EffList = []
        CutVal = "BDTG>="+str(j) 
        CUT = CutVal +" && "+ CUT0

        NsigMC = tSig.GetEntries(CUT)
        NsigEff = NsigMC/float(Ntot)
        EffList.append(NsigEff)
        Deff = ROOT.TMath.Sqrt( float(NsigMC)*(1-( float(NsigMC)/Ntot )) /float(Ntot)**2 + (NsigEff*0.0008/0.2329)**2) #binomial error & efficiency-generated-error
        #ROOT.TMath.Sqrt( float(NsigMC) * (NsigMC+Ntot) / Ntot**3 ) #poisson error
        
        EffList.append(Deff)
        NsigE.update({j : EffList})

        print "NsigEff = {}".format(NsigEff)
        print "NsigEff_error = {}".format(Deff)
        
        #Lumi=2*10**(15);Crosssection-bb=(298+-2+-36)*10**(-6)*2;f_Lb=0.17;BR=3.39*10**(-5)
        NsigYield = NsigEff * 2*10**(15) * 298*10**(-6) * 2 * 0.17 * 3.39*10**(-5)
        NsigYield_error = NsigYield * ROOT.TMath.Sqrt( (0.48/3.39)**2 + (Deff/NsigEff)**2 + (36/298)**2)
        print "NsigYield = {}".format(NsigYield)
        print "NsigYield_error = {}".format(NsigYield_error)
        print "MC error = {}".format(NsigYield*Deff/NsigEff)
        print "BR error = {}".format(NsigYield*0.48/3.39)
        SignalYieldList.append(NsigYield)
        SignalYieldList.append(NsigYield_error)
        NsigY.update({j : SignalYieldList})

    return NsigE, NsigY

def DrawPunzi(Nbkg, NsigE):
    b = 0
    for b in BDTG_val_list2:
        x1.append(b)
        y_val = NsigE[b][0]/( 2.5 + ROOT.TMath.Sqrt(Nbkg[b][0]) )
        #print 'y_val = {}'.format(y_val)
        dp_sur_db = NsigE[b][0]/( (( 2.5 + ROOT.TMath.Sqrt(Nbkg[b][0]) )**2) * (2*ROOT.TMath.Sqrt(Nbkg[b][0])) )
        dp_sur_deff = 1/(2.5 + ROOT.TMath.Sqrt(Nbkg[b][0]))
        y_val_err = ROOT.TMath.Sqrt( (dp_sur_db * Nbkg[b][1])**2 + (dp_sur_deff * NsigE[b][1])**2 )
        #print 'y_val_err = {}'.format(y_val_err)
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
    gP.GetYaxis().SetTitleOffset(1.2)
    cP.Update()
    cP.SaveAs('{}/FigureOfMerit_punzi.pdf'.format(newDir))
    cP.SaveAs('{}/FigureOfMerit_punzi.png'.format(newDir))

def DrawSignificance(Nbkg, NsigY):
    b = 0
    for b in BDTG_val_list2:
        x2.append(b)
        print 'b = {}'.format(b)
        y_val2 = float(NsigY[b][0])/ ROOT.TMath.Sqrt( float(NsigY[b][0]) + Nbkg[b][0] )
        print 'y_val2 = {}'.format(y_val2)
        df_sur_ds = (float(NsigY[b][0]) + 2*Nbkg[b][0])/(2* (ROOT.TMath.Sqrt( float(NsigY[b][0]) + Nbkg[b][0] ))**3 )
        df_sur_db = (float(NsigY[b][0]))/(2* (ROOT.TMath.Sqrt(float(NsigY[b][0]) + Nbkg[b][0]))**3 )
        y_val_err2 = ROOT.TMath.Sqrt((df_sur_ds*float(NsigY[b][1]))**2 + (df_sur_db*Nbkg[b][1])**2)
        #print 'y_val_err2 = {}'.format(y_val_err2)
        #print 'Error on signal events = {}'.format(df_sur_ds*float(NsigY[b][1]))
        #print 'Error on bkg events = {}'.format(df_sur_db*Nbkg[b][1])
        y2.append(y_val2)
        ey2.append(y_val_err2)
        b = b+1
    print max(y2)
    imax = 10
    for i in range(len(BDTG_val_list2)):
        if y2[i] == max(y2):
            print "Maximum of Significance is: {} +- {} at {}".format(y2[i],ey2[i],x2[i])
            print y2
            imax = i
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
    gS.GetYaxis().SetTitle("s/#sqrt{s+b}")
    gS.GetYaxis().SetTitleOffset(1.2)
    cS.Update()
    cS.SaveAs('{}/FigureOfMerit_signif.pdf'.format(newDir))
    cS.SaveAs('{}/FigureOfMerit_signif.png'.format(newDir))
    return imax

def DrawBackgroundYield(Yield, imax):
    b = 0
    for b in BDTG_val_list2:
        xb.append(b)
        yb.append(Yield[b][0])
        eyb.append(Yield[b][1])
        b = b+1
    print "Background yield at imax = {} is {}+-{}".format(imax,yb[imax],eyb[imax])   
    xB = array('d',xb)
    yB = array('d',yb)
    exB = array('d',exb)
    eyB = array('d',eyb)

    cP = ROOT.TCanvas('cP','cP')
    ROOT.gROOT.SetStyle('Plain')
    gP = ROOT.TGraphErrors(n,xB,yB,exB,eyB)
    gP.SetMarkerStyle(20)
    gP.SetMarkerColor(4)
    gP.SetLineColor(4)
    gP.Draw()
    gP.GetXaxis().SetTitle("BDTG cut value")
    gP.GetYaxis().SetTitle("N_{bkg}")
    gP.GetYaxis().SetTitleOffset(1.2)
    cP.Update()
    cP.SaveAs('{}/Bkg_Yield.pdf'.format(newDir))
    cP.SaveAs('{}/Bkg_Yield.png'.format(newDir))

def DrawSignalYield(Yield,imax):
    i = 0
    for i in BDTG_val_list2:
        xs.append(i)
        ys.append(Yield[i][0])
        eys.append(Yield[i][1])
        i = i+1
    
    print "Signal yield at imax = {} is {}+-{}".format(imax,ys[imax],eys[imax])    
    xS = array('d',xs)
    yS = array('d',ys)
    exS = array('d',exs)
    eyS = array('d',eys)

    c = ROOT.TCanvas('c','c')
    ROOT.gROOT.SetStyle('Plain')
    g = ROOT.TGraphErrors(n,xS,yS,exS,eyS)
    g.SetMarkerStyle(20)
    g.SetMarkerColor(4)
    g.SetLineColor(4)
    g.Draw()
    g.GetXaxis().SetTitle("BDTG cut value")
    g.GetYaxis().SetTitle("N_{signal}")
    g.GetYaxis().SetTitleOffset(1.2)
    c.Update()
    c.SaveAs('{}/Signal_Yield.pdf'.format(newDir))
    c.SaveAs('{}/Signal_Yield.png'.format(newDir))


#____________________________________MAIN_____________________________________


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    #parser.add_argument("file", action="store", type=str)
    #parser.add_argument("-t", "--tree", default="DecayTree", action="store", type=str)
    #(xmini)______FIT______(xmin)______SIGNAL______(xmid)_______FIT_______(xmax)
    parser.add_argument("-m", "--xmini", default=4420., action="store", type=float)
    parser.add_argument("-n", "--xmin", default=5120., action="store", type=float)
    parser.add_argument("-d", "--xmid", default=6120., action="store", type=float)
    parser.add_argument("-x", "--xmax", default=6620., action="store", type=float)
    #parser.add_argument("-c", "--cuts", default="", action="store", type=str)
    parser.add_argument("-v", "--version", default="v5_all",action="store", type=str)
    args = parser.parse_args()

    newDir = '/sps/lhcb/volle/FitResults/{}'.format(args.version)
    if not os.path.isdir(newDir):
        os.makedirs(newDir)

    #Save for all BDTG_cut_values the Number of bkg events, error & chi2 of fit
    Nbkg = GetBackgroundYieldDict(t,args.version,args.xmini,args.xmin,args.xmid,args.xmax) 

    #Save for all BDTG_cut_values the Number of Signal events & their error
    NsigE, NsigY =  GetSignalEfficiencyAndYieldDict(fSig,tSig,args.xmin,args.xmid)
    
    #Draw Punzi Figure of merit
    DrawPunzi(Nbkg, NsigE)
    
    #Draw Significance Figure of merit
    imax = DrawSignificance(Nbkg, NsigY)
    
    #Draw background yield distribution
    DrawBackgroundYield(Nbkg, imax)
    DrawSignalYield(NsigY,imax)



  
        
        
        


    


