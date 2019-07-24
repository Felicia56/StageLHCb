#!/usr/bin/env python
# =============================================================================
# @file   FitBMass.py
# @author Felicia Volle
# @date   23.07.2019
# =============================================================================
"""This script fits a simple pdf to a given dataset using RooFit."""

# imports
import os,sys
#from uncertainties import ufloat
import ROOT
from ROOT import kTRUE

RooFit         = ROOT.RooFit
RooRealVar     = ROOT.RooRealVar
RooArgList     = ROOT.RooArgList
RooArgSet      = ROOT.RooArgSet
RooDataSet     = ROOT.RooDataSet
RooGaussian    = ROOT.RooGaussian
RooExponential = ROOT.RooExponential
RooAddPdf      = ROOT.RooAddPdf
RooExtendPdf   = ROOT.RooExtendPdf

name = "BifurcatedCB.cxx"
ROOT.gSystem.CompileMacro(name, 'k')
ROOT.gSystem.Load(name.replace('.c', '_c'))

# alle prints auskommentiert -> bei FigureOfMerit kein ewig langer Text

# definition of functions for this script
def FitBMassSigAndBkg(tree, treeSig, version, cuts,meanV, xmini, xmin, xmid, xmax): #mean, xmin = 6120, xmax = 9000):
    """
    This function fits the "B_M" variable of a given TTree
    with a model formed by a Gaussian and an exponential pdf.
    All shape parameters are allowed to float in the fit. The 
    initial and range values are hardcoded in the code, except
    for the initial value of the Gaussian mean and the range
    of the B_M variable to be used.
    Returns the composed model (RooAbsPdf), the residual plot object
    and its chi2 value (float)
    Definition of the arguments:
    :tree: type TTree
    the root TTree that contains the variable to be fitted
    :cuts: type str
    optional cuts to apply to the TTree before fitting
    :mean: type float
    initial value for the Gaussian mean that will be floated
    during the fit
    :xmin: type float, optional
    minimum value of the B_M range to be fitted. Default: 4000
    :xmax: type float, optional
    maximum value of the B_M range to be fitted. Default: 7000
    """
    
    #___define variables and pdfs_______________________________________________


    if xmini > 0:
        B_M = RooRealVar("B_M","B_M", xmini, xmax)
        B_M.setRange("SignalRegion",xmin,xmid)
        B_M.setRange("FullRange", xmini, xmax)
    else:
        B_M = RooRealVar("B_M","B_M", xmin, xmax)
        B_M.setRange("SignalRegion",xmin,xmid)
        B_M.setRange("FullRange", xmin, xmax)
 
    #mean  = RooRealVar("mean", "mean",  meanV, 5020, 6220)
    #sigma = RooRealVar("sigma", "sigma", 100, 0, 200)
    #gauss = RooGaussian("gauss", "gauss", B_M, mean, sigma)

    # Gaussian core
    mean  = RooRealVar("mean", "mean",  meanV, 5020, 6220)#, PlotLabel="#mu")
    sigma = RooRealVar("sigma", "sigma", 100, 1, 200)#, PlotLabel="#sigma")
    # Low tail
    alphaL= RooRealVar("alphaL", "alphaL", 2, 0, 10)#, PlotLabel="#alpha_{L}")
    nL = RooRealVar("nL", "nL", 5, 0, 10)#, PlotLabel="n_{L}")
    # High tail
    alphaH = RooRealVar("alphaH", "alphaH", 2, 0, 10)#, PlotLabel="#alpha_{R}")
    nH = RooRealVar("nH", "nH", 5, 0, 15)#, PlotLabel="n_{R}")

    SignalModel = ROOT.BifurcatedCB("SignalModel", "SignalModel", B_M, mean, sigma, alphaL, nL, alphaH, nH)


    # Bkg - exponential
    tau = RooRealVar("tau", "tau", -0.0033, -0.1, 0.)
    exp = RooExponential("exp", "exp", B_M, tau)

    
    #___define coefficiencts___________________________________________________


    nsig = RooRealVar("nsig", "nsig", 10000, 10, 500000)#value, MinValue, MaxValue
    nsigMC = RooRealVar("nsig", "nsig", 10000, 10, 500000)#value, MinValue, MaxValue
    nbkg = RooRealVar("nbkg", "nbkg", 4241, 0, 500000)
    
    #___define dataset________________________________________________________


    outDir = '/sps/lhcb/volle/FitBMass'
    newDir = '/sps/lhcb/volle/FitBMass/{}'.format(version)
    a,b = cuts.split("=")

    if not os.path.isdir(newDir):
        os.makedirs(newDir)
        
    CUTS = cuts + " && B_M>{} && B_M<{}".format(args.xmini,args.xmax) 

    fnew = ROOT.TFile("{}/Dataset_{}.root".format(newDir,b),"recreate")
    tnew = tree.CopyTree(CUTS)
    fnew.Write()
    
    ds = RooDataSet("data", "dataset with x", tnew, RooArgSet(B_M))

    #fds = ROOT.TFile("{}/Dataset.root".format(newDir,b),"recreate")
    #tds = ds.store().tree()
    #fds.Write()
    
    fnewMC = ROOT.TFile("{}/BDTG_MCfile_{}.root".format(newDir,b),"recreate")
    tnewMC = treeSig.CopyTree(cuts)
    fnewMC.Write()
    
    dsMC = RooDataSet("data", "dataset with x", tnewMC, RooArgSet(B_M))

    
    #___MC fit and fixing tails_________________________________________________

    SigModel = RooAddPdf("SigModel","SigModel",RooArgList(SignalModel),RooArgList(nsigMC))

    c0 = ROOT.TCanvas("c0","c0")
    
    massFrameMC = B_M.frame(xmini,xmax,100)
    massFrameMC.SetTitle("MC fit for BDTG variable larger or equal to {}".format(b))
    dsMC.plotOn(massFrameMC)
    
    fitResultsMC = SigModel.fitTo(dsMC, RooFit.Range("FullRange"), RooFit.Save())#fits model to B_M by using LS and RS

    SigModel.plotOn(massFrameMC, RooFit.FillColor(4), RooFit.VisualizeError(fitResultsMC, 1), RooFit.Name("curve_model"))

    SigModel.paramOn(massFrameMC, RooFit.Layout(0.59,0.97,0.92), RooFit.AutoPrecision(1) ) #RooFit.Parameters( RooArgSet(tau) ) ) #RooFit.Label("Fit results")

    massFrameMC.Draw()
    massFrameMC.GetXaxis().SetTitle("M(#Lambda_{b})")
    
    c0.Update()

    c0.SaveAs("{}/MCfit_BDTG_cut_{}.pdf".format(newDir,b))
    c0.SaveAs("{}/MCfit_BDTG_cut_{}.png".format(newDir,b))
    #c0.SaveAs("{}/MCfit_BDTG_cut_{}.root".format(newDir,b)

    #Set tails constant
    alphaL.setConstant(kTRUE)
    nL.setConstant(kTRUE)

    alphaH.setConstant(kTRUE)
    nH.setConstant(kTRUE)


    #___build model____________________________________________________________


    suma = RooArgList()
    coeff = RooArgList()
    
    #suma.add(gauss)
    suma.add(SignalModel)
    suma.add(exp)
    
    coeff.add(nsig)
    coeff.add(nbkg)
    
    model = RooAddPdf("model", "model", suma, coeff)

    
    #___plot dataset and fit__________________________________________________

    
    c1 = ROOT.TCanvas("c1","c1")
    
    massFrame = B_M.frame(xmini,xmax,100)
    massFrame.SetTitle("BDTG variable larger or equal to {}".format(b))
    ds.plotOn(massFrame)
    
    fitResults = model.fitTo(ds, RooFit.Range("FullRange"), RooFit.Save())#fits model to B_M by using LS and RS

    model.plotOn(massFrame, RooFit.Components("exp"), RooFit.Range("FullRange"), RooFit.NormRange("FullRange"), RooFit.FillColor(3), RooFit.VisualizeError(fitResults, 1))
    model.plotOn(massFrame, RooFit.FillColor(6), RooFit.VisualizeError(fitResults, 1), RooFit.Name("curve_model"))
    ###model.plotOn(massFrame, RooFit.Components("gauss"), RooFit.LineColor(3),RooFit.VisualizeError(fitResults, 1))
    
    #model.plotOn(massFrame, RooFit.Components("gauss"), RooFit.Range("FullRange"), RooFit.NormRange("FullRange"), RooFit.FillColor(6), RooFit.VisualizeError(fitResults, 1))
    
    model.paramOn(massFrame, RooFit.Layout(0.59,0.97,0.92), RooFit.AutoPrecision(1) ) #RooFit.Parameters( RooArgSet(tau) ) ) #RooFit.Label("Fit results")
    #paramList = RooArgSet(tau, nbkg)
    #paramList.printLatex(RooFit.Format("NEU"), RooFit.AutoPrecision(1), RooFit.VerbatimName())
    #ds.statOn(massFrame)
    #model.paramOn(massFrame, RooFit.Layout=(0.54,0.92,0.92), Parameters = RooArgSet(nbkg_sig, tau))#nsig, nbkg, mean, sigma, tau))
    #chi2 = massFrame.chiSquare()

    massFrame.Draw()

    massFrame.GetXaxis().SetTitle("M(#Lambda_{b})")
    
    c1.Update()

    c1.SaveAs("{}/BDTG_cut_{}.pdf".format(newDir,b))
    c1.SaveAs("{}/BDTG_cut_{}.png".format(newDir,b))
    #c1.SaveAs("{}/BDTG_cut_{}.root".format(newDir,b)
    
    #Chi2 calculation of Fit
    c2 = ROOT.TCanvas("c2","c2")
    
    if xmini>0:
        Frame = B_M.frame(xmini,xmax,100)#60 for xmid = 6120 & xmax = 9120 same bin size than model 
        ds.plotOn(Frame)
        model.plotOn(Frame, RooFit.Range("FullRange"), RooFit.NormRange("FullRange"))
    else:
        Frame = B_M.frame(xmid,xmax,200)#60 for xmid = 6120 & xmax = 9120 same bin size than model 
        ds.plotOn(Frame)
        model.plotOn(Frame, RooFit.Range("FullRange"), RooFit.NormRange("FullRange"))
    Frame.Draw()
    chi2 = Frame.chiSquare()

    #c2.SaveAs("{}/BDTG_cut_{}.pdf".format(newDir,b))
    #c2.SaveAs("{}/BDTG_cut_{}.root".format(newDir,b)
        

    #___print results__________________________________________________________


    print "{} has been fit to {} with a chi2 = {}".format(model.GetName(),tnew.GetName(), chi2)
    print "Total number of entries in the data is: {}".format(ds.numEntries())
    print "Number of sig entries is: {:.0f} +- {:.0f}".format(nsig.getValV(),nsig.getError())
    print "Fit result bkg is: {:.0f} +- {:.0f} in the full range by fitting the interval [{},{}]".format(nbkg.getValV(),nbkg.getError(),xmin,xmax)
    

    #___compute S/B with error propagation from uncertainties module___________


    #B_M.setRange("SignalRegion",5120,6120)
    l = RooArgSet(B_M)


    Nbkg_sig_fraction = exp.createIntegral(l,  RooFit.NormSet(l), RooFit.Range("SignalRegion") )
    Nbkg_sig_fraction_error = Nbkg_sig_fraction.getPropagatedError(fitResults)
#    print "Fraction of bkg events in the signal region is: {} +/- {}".format(Nbkg_sig_fraction.getVal(),Nbkg_sig_fraction_error)
    nbkg_sig = Nbkg_sig_fraction.getVal() * nbkg.getVal() 
    nbkg_sig_error = nbkg_sig * ROOT.TMath.Sqrt( (nbkg.getError()/nbkg.getValV())**2 + (Nbkg_sig_fraction_error/Nbkg_sig_fraction.getVal())**2 ) 
    print "Number of bkg events in the signal region is: {:.0f} +/- {:.0f}".format(nbkg_sig,nbkg_sig_error)

    signif = nsig.getValV()/nbkg_sig 
    signif_err = signif * ROOT.TMath.Sqrt( (nbkg_sig_error/nbkg_sig)**2 + (nsig.getError()/nsig.getValV())**2 )
    print "S/B = {:.2f} +/- {:.2f}".format(signif, signif_err)
  
    #_____splot_____________________________________________________________

    print "Starting the splot analysis"
    
    mean.setConstant(kTRUE)
    sigma.setConstant(kTRUE)
    tau.setConstant(kTRUE)

    sData = ROOT.RooStats.SPlot("sData","sData", ds, model, coeff)
    print "Finished to create sData"
    
    
    #Enable to export RooDataset as a TTree
    ROOT.RooDataSet.setDefaultStorageType(ROOT.RooAbsData.Tree)

    sFile = ROOT.TFile("sData.root","RECREATE")
    sFile.cd()
    
    sTree = ROOT.RooStats.GetAsTTree('sTree','sTree',sData.GetSDataSet())
    sTree.Print()

    #Check sWeights
    print "Check SWeights:"
    print "Nsig = {} and from sWeights = {}".format(nsig.getVal(),sData.GetYieldFromSWeight("nsig"))
    print "Nbkg = {} and from sWeights = {}".format(nbkg.getVal(),sData.GetYieldFromSWeight("nbkg"))

    #Go back to default storage type
    ROOT.RooDataSet.setDefaultStorageType(ROOT.RooAbsData.Vector)
    sTree.Write()
    sFile.Close()

    return ds, model #nsig, nbkg_sig, nbkg_sig_error, nbkg_sig2, nbkg_sig_error2


#______MAIN FUNCTION___________________________________________________________


if __name__=="__main__":
    import argparse
    parser = argparse.ArgumentParser()
    #parser.add_argument("file", action="store", type=str)
    #parser.add_argument("-t", "--tree", default="DecayTree", action="store", type=str)
    #(xmini)______FIT______(xmin)______SIGNAL______(xmid)_______FIT_______(xmax)
    parser.add_argument("-a", "--mean", default=5620, action="store", type=float)
    parser.add_argument("-m", "--xmini", default=4820., action="store", type=float)#4420
    parser.add_argument("-n", "--xmin", default=5120., action="store", type=float)
    parser.add_argument("-d", "--xmid", default=6120., action="store", type=float)
    parser.add_argument("-x", "--xmax", default=6420., action="store", type=float)#6620
    parser.add_argument("-c", "--cuts", default="", action="store", type=str)
    parser.add_argument("-v", "--version", default="v1_Try2", action="store", type=str)
    args = parser.parse_args()

    # sanity check
    #if not os.path.exists(args.file):
    #    print "File doesn't exist! Exiting..."
    #    exit()

    # read data
    f = ROOT.TFile("/sps/lhcb/volle/Data_BDTG_Secure.root")
    t = f.Get("DecayTree")
    fSig = ROOT.TFile("/sps/lhcb/volle/MCsignal_BDTG_Secure.root")
    tSig = fSig.Get("DecayTree")

    #f = ROOT.TFile("/sps/lhcb/volle/Data_all_BDTG.root")
    #t = f.Get("DecayTree")

    #meanVal = 5620
    BDTCut = "BDTG>=0.92"#0.72 #0.92

    ds, model = FitBMassSigAndBkg(t, tSig, args.version,  BDTCut,args.mean, args.xmini, args.xmin, args.xmid, args.xmax) #args.mean, args.xmin, args.xmax)#ds,model,Plot,chi2
    

#EOF

