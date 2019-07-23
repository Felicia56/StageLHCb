#!/usr/bin/env python
# =============================================================================
# @file   simpleFit.py
# @author Felicia Volle
# @date   23.07.2019
# =============================================================================
"""This script fits a simple pdf to a given dataset using RooFit."""

# imports
import os
#from uncertainties import ufloat
import ROOT

RooFit         = ROOT.RooFit
RooRealVar     = ROOT.RooRealVar
RooArgList     = ROOT.RooArgList
RooArgSet      = ROOT.RooArgSet
RooDataSet     = ROOT.RooDataSet
###RooGaussian    = ROOT.RooGaussian
RooExponential = ROOT.RooExponential
RooAddPdf      = ROOT.RooAddPdf
RooExtendPdf   = ROOT.RooExtendPdf

# alle prints auskommentiert -> bei FigureOfMerit kein ewig langer Text

# definition of functions for this script
def simpleFit(tree, version, cuts, xmini, xmin, xmid, xmax): #mean, xmin = 6120, xmax = 9000):
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
    
    #___define variables and pdfs______________________________________________
 
    if xmini > 0:
        B_M = RooRealVar("B_M","B_M", xmini, xmax)#, 3620, 7120)
        B_M.setRange("SignalRegion",xmin,xmid)
        B_M.setRange("FullRange", xmini, xmax)
        B_M.setRange("LeftSideBand", xmini, xmin)
        B_M.setRange("RightSideBand", xmid, xmax)
    else:
        B_M = RooRealVar("B_M","B_M", xmin, xmax)
        B_M.setRange("SignalRegion",xmin,xmid)
        B_M.setRange("FullRange", xmin, xmax)
        B_M.setRange("OuterRange",xmid,xmax)

    ###mean  = RooRealVar("mean", "mean",  mean, mean-50, mean+50)
    ###sigma = RooRealVar("sigma", "sigma", 80, 10, 150)
    ###gauss = RooGaussian("gauss", "gauss", B_M, mean, sigma)
    
    tau = RooRealVar("tau", "tau", -0.01, -0.01, 0.)
    exp = RooExponential("exp", "exp", B_M, tau)

    
    #___define coefficiencts___________________________________________________


    #nsig = RooRealVar("nsig", "nsig", 1000, 0, 20000)#value, MinValue, MaxValue
    nbkg = RooRealVar("nbkg", "nbkg", 20000, 0, 1000000)
    

    #___build model____________________________________________________________


    suma = RooArgList()
    coeff = RooArgList()
    
    ###suma.add(gauss)
    suma.add(exp)
    
    ###coeff.add(nsig)
    coeff.add(nbkg)
    
    model = RooAddPdf("model", "model", suma, coeff)

    
    #___define dataset________________________________________________________


    outDir = '/sps/lhcb/volle/FitResults'
    newDir = '/sps/lhcb/volle/FitResults/{}'.format(version)
    a,b = cuts.split("=")

    if not os.path.isdir(newDir):
        os.makedirs(newDir)
        

    fnew = ROOT.TFile("{}/BDTG_BKGfile_{}.root".format(newDir,b),"recreate")
    tnew = tree.CopyTree(cuts)
    fnew.Write()
    
    ds = RooDataSet("data", "dataset with x", tnew, RooArgSet(B_M))

    
    #___plot dataset and fit__________________________________________________

    
    c1 = ROOT.TCanvas("c1","c1")
    
    massFrame = B_M.frame()
    massFrame.SetTitle("BDTG variable larger or equal to {}".format(b))
    ds.plotOn(massFrame)
    
    if xmini > 0:
        fitResults = model.fitTo(ds, RooFit.Range("LeftSideBand,RightSideBand"), RooFit.Save())#fits model to B_M by using LS and RS
    else:
        fitResults = model.fitTo(ds, RooFit.Range("OuterRange"), RooFit.Save())#fits model to B_M
    #model.plotOn(massFrame, RooFit.LineColor(2), RooFit.VisualizeError(fitResults, 1), RooFit.Name("curve_model"))
    ###model.plotOn(massFrame, RooFit.Components("gauss"), RooFit.LineColor(3),RooFit.VisualizeError(fitResults, 1))
    
    if xmini>0:
        model.plotOn(massFrame, RooFit.Components("exp"), RooFit.Range("FullRange"), RooFit.NormRange("LeftSideBand,RightSideBand"), RooFit.FillColor(9), RooFit.VisualizeError(fitResults, 1))
        model.plotOn(massFrame, RooFit.Components("exp"), RooFit.Range("LeftSideBand,RightSideBand"), RooFit.NormRange("LeftSideBand,RightSideBand"), RooFit.FillColor(2), RooFit.VisualizeError(fitResults, 1))
    else:
        model.plotOn(massFrame, RooFit.Components("exp"), RooFit.Range("FullRange"), RooFit.NormRange("OuterRange"), RooFit.FillColor(9), RooFit.VisualizeError(fitResults, 1)) 
        model.plotOn(massFrame, RooFit.Components("exp"), RooFit.Range("OuterRange"), RooFit.NormRange("OuterRange"), RooFit.FillColor(2), RooFit.VisualizeError(fitResults, 1))
    
    model.paramOn(massFrame, RooFit.Layout(0.59,0.97,0.92), RooFit.AutoPrecision(1) ) #RooFit.Parameters( RooArgSet(tau) ) ) #RooFit.Label("Fit results")
    #paramList = RooArgSet(tau, nbkg)
    #paramList.printLatex(RooFit.Format("NEU"), RooFit.AutoPrecision(1), RooFit.VerbatimName())
    #ds.statOn(massFrame)
    #model.paramOn(massFrame, RooFit.Layout=(0.54,0.92,0.92), Parameters = RooArgSet(nbkg_sig, tau))#nsig, nbkg, mean, sigma, tau))
    #chi2 = massFrame.chiSquare()

    massFrame.Draw()

    #Chi2 calculation of Fit
    if xmini>0:
        Frame = B_M.frame(xmini,xmax,200)#60 for xmid = 6120 & xmax = 9120 same bin size than model 
        ds.plotOn(Frame)
        model.plotOn(Frame, RooFit.Components("exp"), RooFit.Range("LeftSideBand,RightSideBand"), RooFit.NormRange("LeftSideBand,RightSideBand"))
    else:
        Frame = B_M.frame(xmid,xmax,200)#60 for xmid = 6120 & xmax = 9120 same bin size than model 
        ds.plotOn(Frame)
        model.plotOn(Frame, RooFit.Components("exp"), RooFit.Range("OuterRange"), RooFit.NormRange("OuterRange"))
    chi2 = Frame.chiSquare()


    #___print results__________________________________________________________


#    print "{} has been fit to {} with a chi2 = {}".format(model.GetName(),tnew.GetName(), chi2)
#    print "Total number of entries in the data is: {}".format(ds.numEntries())
    ###print "Number of sig entries is: {:.0f} +- {:.0f}".format(nsig.getValV(),nsig.getError())
#    print "Fit result bkg is: {:.0f} +- {:.0f} in the full range by fitting the interval [{},{}]".format(nbkg.getValV(),nbkg.getError(),xmin,xmax)
    

    #___compute S/B with error propagation from uncertainties module___________

    ###nsig = ufloat(nsig.getValV(), nsig.getError())
    ###nbkg = ufloat(nbkg.getValV(), nbkg.getError())
    ###signif = nsig/nbkg
    ###print "S/B = {:.2f} +/- {:.2f}".format(signif.nominal_value, signif.std_dev)
    
    #___calculate integral in signal region____________________________________


    #B_M.setRange("SignalRegion",5120,6120)
    l = RooArgSet(B_M)


    Nbkg_sig_fraction = model.createIntegral(l, RooFit.NormSet(l), RooFit.Range("SignalRegion") )
    Nbkg_sig_fraction_error = Nbkg_sig_fraction.getPropagatedError(fitResults)
#    print "Fraction of bkg events in the signal region is: {} +/- {}".format(Nbkg_sig_fraction.getVal(),Nbkg_sig_fraction_error)
    nbkg_sig = Nbkg_sig_fraction.getVal() * nbkg.getVal() 
    nbkg_sig_error = nbkg_sig * ROOT.TMath.Sqrt( (nbkg.getError()/nbkg.getValV())**2 + (Nbkg_sig_fraction_error/Nbkg_sig_fraction.getVal())**2 ) 
#    print "Number of bkg events in the signal region is: {:.0f} +/- {:.0f}".format(nbkg_sig,nbkg_sig_error)
        
    c1.SaveAs("{}/BDTG_cut_{}.pdf".format(newDir,b))
    c1.SaveAs("{}/BDTG_cut_{}.png".format(newDir,b))
    #c1.SaveAs("{}/BDTG_cut_{}.root".format(newDir,b)
    
    return ds, model, chi2, nbkg_sig, nbkg_sig_error


#______MAIN FUNCTION___________________________________________________________


if __name__=="__main__":
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("file", action="store", type=str)
    #parser.add_argument("-t", "--tree", default="DecayTree", action="store", type=str)
    #(xmini)______FIT______(xmin)______SIGNAL______(xmid)_______FIT_______(xmax)
    parser.add_argument("-m", "--xmini", default=-1., action="store", type=float)
    parser.add_argument("-n", "--xmin", default=5120., action="store", type=float)
    parser.add_argument("-d", "--xmid", default=6120., action="store", type=float)
    parser.add_argument("-x", "--xmax", default=9120., action="store", type=float)
    parser.add_argument("-c", "--cuts", default="", action="store", type=str)
    parser.add_argument("-v", "--version", default="v1", action="store", type=str)
    args = parser.parse_args()

    # sanity check
    if not os.path.exists(args.file):
        print "File doesn't exist! Exiting..."
        exit()

    # read data
    f = ROOT.TFile(args.file)
    t = f.Get("DecayTree")

    ds, model, chi2, Nbkg_SR, Nbkg_SR_error = simpleFit(t, args.version,  args.cuts, args.xmini, args.xmin, args.xmid, args.xmax) #args.mean, args.xmin, args.xmax)#ds,model,Plot,chi2
    

#EOF

