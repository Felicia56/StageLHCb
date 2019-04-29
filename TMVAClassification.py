#!/usr/bin/env python
# @(#)root/tmva $Id: TMVAClassification.py 38475 2011-03-17 10:46:00Z evt $
# ------------------------------------------------------------------------------ #
# Project      : TMVA - a Root-integrated toolkit for multivariate data analysis #
# Package      : TMVA                                                            #
# Python script: TMVAClassification.py                                           #
#                                                                                #
# This python script provides examples for the training and testing of all the   #
# TMVA classifiers through PyROOT.                                               #
#                                                                                #
# The Application works similarly, please see:                                   #
#    TMVA/macros/TMVAClassificationApplication.C                                 #
# For regression, see:                                                           #
#    TMVA/macros/TMVARegression.C                                                #
#    TMVA/macros/TMVARegressionpplication.C                                      #
# and translate to python as done here.                                          #
#                                                                                #
# As input data is used a toy-MC sample consisting of four Gaussian-distributed  #
# and linearly correlated input variables.                                       #
#                                                                                #
# The methods to be used can be switched on and off via the prompt command, for  #
# example:                                                                       #
#                                                                                #
#    python TMVAClassification.py --methods Fisher,Likelihood                    #
#                                                                                #
# The output file "TMVA.root" can be analysed with the use of dedicated          #
# macros (simply say: root -l <../macros/macro.C>), which can be conveniently    #
# invoked through a GUI that will appear at the end of the run of this macro.    #
#                                                                                #
# for help type "python TMVAClassification.py --help"                            #
# ------------------------------------------------------------------------------ #

# --------------------------------------------
# Standard python import
import sys    # exit
import time   # time accounting
import getopt # command line parser
import os

# --------------------------------------------

# Print usage help
def usage():
    # Default settings for command line arguments
    #DEFAULT_OUTFNAME = "TMVA.root"
    DEFAULT_OPTNAME  = ""
    DEFAULT_SIGFNAME = "Lb2L0Gamma_presel_sig_mcmatch_LL_presel.root"
    DEFAULT_BKGFNAME = "Lb2L0Gamma_presel_bkg_LL_presel.root"
    DEFAULT_CHANNEL  = "LL"
    DEFAULT_TREESIG  = "Lb2L0GammaTuple/DecayTree"
    DEFAULT_TREEBKG  = "Lb2L0GammaTuple/DecayTree"
    DEFAULT_METHODS  = "BDT"
    print " "
    print "Usage: python %s [options]" % sys.argv[0]
    print "  -m | --methods    : gives methods to be run (default: %s)" % DEFAULT_METHODS
    print "  -s | --sigfile    : name of signal ROOT file (default: '%s')" % DEFAULT_SIGFNAME
    print "  -b | --bkgfile    : name of bkg ROOT file (default: '%s')" % DEFAULT_BKGFNAME
    #print "  -o | --outputfile : name of output ROOT file containing results (default: '%s')" % DEFAULT_OUTFNAME
    print "  -o | --optname    : name of the options for training to be used in the output files (default: '%s')" % DEFAULT_OPTNAME
    print "  -c | --channel    : name of the channel to be trained (default: '%s')" % DEFAULT_CHANNEL
    print "  -t | --inputtrees : input ROOT Trees for signal and background (default: '%s %s')" \
          % (DEFAULT_TREESIG, DEFAULT_TREEBKG)
    print "  -v | --verbose"
    print "  -? | --usage      : print this help message"
    print "  -h | --help       : print this help message"
    print " "

# Main routine
def TMVAClassification(methods, sigfname, bkgfname, optname, channel, trees="DecayTree,DecayTree", verbose=False):
    # Print methods
    mlist = methods.replace(' ',',').split(',')
    print "=== TMVAClassification: use method(s)..."
    for m in mlist:
        if m.strip() != '':
            print "=== - <%s>" % m.strip()
    # Define trees
    trees = trees.split(",")
    if len(trees)-trees.count('') != 2:
        print "ERROR: need to give two trees (each one for signal and background)"
        print trees
        sys.exit(1)
    treeNameSig = trees[0]
    treeNameBkg = trees[1]
    # Print output file and directory
    outfname = "TMVA_%s_%s.root" %(channel, optname)
    myWeightDirectory = "weights_%s_%s" %(channel, optname)
    print "=== TMVAClassification: output will be written to:"
    print "=== %s" %outfname
    print "=== %s" %myWeightDirectory

    # Import ROOT classes
    from ROOT import gSystem, gROOT, gApplication, TFile, TTree, TCut
    
    # check ROOT version, give alarm if 5.18 
    if gROOT.GetVersionCode() >= 332288 and gROOT.GetVersionCode() < 332544:
        print "*** You are running ROOT version 5.18, which has problems in PyROOT such that TMVA"
        print "*** does not run properly (function calls with enums in the argument are ignored)."
        print "*** Solution: either use CINT or a C++ compiled version (see TMVA/macros or TMVA/examples),"
        print "*** or use another ROOT version (e.g., ROOT 5.19)."
        sys.exit(1)
    
    # Logon not automatically loaded through PyROOT (logon loads TMVA library) load also GUI
    #gROOT.SetMacroPath( "./" )
    #gROOT.Macro       ( "./tmva/test/TMVAlogon.C" )
    #gROOT.LoadMacro   ( "./tmva/test/TMVAGui.C" ) ###Is this really necessary??
    
    # Import TMVA classes from ROOT
    from ROOT import TMVA

    # Output file
    outputFile = TFile( outfname, 'RECREATE' )
    
    # Create instance of TMVA factory (see TMVA/macros/TMVAClassification.C for more factory options)
    # All TMVA output can be suppressed by removing the "!" (not) in 
    # front of the "Silent" argument in the option string
    factory = TMVA.Factory( "TMVAClassification", outputFile, 
                            "!V:!Silent:Color:DrawProgressBar:Transformations=I:AnalysisType=Classification" )

    # Set verbosity
    factory.SetVerbose( verbose )
    
    # Load data
    dataloader = TMVA.DataLoader("dataset_%s" % channel)
    
    # If you wish to modify default settings 
    # (please check "src/Config.h" to see all available global options)
    #    gConfig().GetVariablePlotting()).fTimesRMS = 8.0
    (TMVA.gConfig().GetIONames()).fWeightFileDir = myWeightDirectory
    # Define the input variables that shall be used for the classifier training
    # note that you may also use variable expressions, such as: "3*var1/var2*abs(var3)"
    # [all types of expressions that can also be parsed by TTree::Draw( "expression" )]

    print "*** Training on channel:"
    print "*** %s" % channel
    print "***"

    if channel == "LL":
        dataloader.AddVariable( "piminus_PT",                         "P_{T}(#pi^{-})",                                "MeV",  'F' );
        dataloader.AddVariable( "pplus_PT",                           "P_{T}(p^{+})",                                  "MeV",  'F' );
        dataloader.AddVariable( "piminus_IPCHI2_OWNPV",               "IP #chi^{2}(#pi^{-})",                         ""    ,  'F' );
        dataloader.AddVariable( "pplus_IPCHI2_OWNPV",                 "IP #chi^{2}(p^{+})",                             ""  ,  'F' );
        dataloader.AddVariable( "Lambda_b0_DOCA_12",                  "DOCA(p-#pi)",                              "mm",  'F' );
        dataloader.AddVariable( "gamma_PT",                           "PT(#gamma)",                              "MeV",  'F' );
        dataloader.AddVariable( "Lambda0_PT",                         "P_{T}(#Lambda)",                                  "MeV",  'F' );
        dataloader.AddVariable( "Lambda0_IP_OWNPV",                   "IP(#Lambda)",                             "mm",     'F' );
        dataloader.AddVariable( "Lambda0_IPCHI2_OWNPV",               "IP#chi^{2}(#Lambda)",                             "",     'F' );
        dataloader.AddVariable( "Lambda0_FDCHI2_OWNPV",                   "FD #chi^{2}(#Lambda)",                             "",     'F' );
        dataloader.AddVariable( "Lambda_b0_PT",                               "P_{T}(#Lambda_{b})",                                      "MeV",  'F' );
        dataloader.AddVariable( "Lambda_b0_MTDOCA",                   "MTDOCA(#Lambda_{b})",                                      "mm",  'F' );

    if channel == "DD":
        dataloader.AddVariable( "piminus_PT",                         "P_{T}(#pi^{-})",                                "MeV",  'F' );
        dataloader.AddVariable( "pplus_PT",                           "P_{T}(p^{+})",                                  "MeV",  'F' );
        dataloader.AddVariable( "piminus_IPCHI2_OWNPV",               "IP #chi^{2}(#pi^{-})",                         ""    ,  'F' );
        dataloader.AddVariable( "pplus_IPCHI2_OWNPV",                 "IP #chi^{2}(p^{+})",                             ""  ,  'F' );
        dataloader.AddVariable( "Lambda_b0_DOCA_12",                  "DOCA(p-#pi)",                              "mm",  'F' );
        dataloader.AddVariable( "gamma_PT",                           "PT(#gamma)",                              "MeV",  'F' );
        dataloader.AddVariable( "Lambda0_PT",                         "P_{T}(#Lambda)",                                  "MeV",  'F' );
        dataloader.AddVariable( "Lambda0_IP_OWNPV",                   "IP(#Lambda)",                             "mm",     'F' );
        #dataloader.AddVariable( "Lambda0_IPCHI2_OWNPV",               "IP#chi^{2}(#Lambda)",                             "",     'F' );
        dataloader.AddVariable( "Lambda0_FD_OWNPV",                   "FD (#Lambda)",                             "mm",     'F' );
        dataloader.AddVariable( "Lambda_b0_PT",                               "P_{T}(#Lambda_{b})",                                      "MeV",  'F' );
        dataloader.AddVariable( "Lambda_b0_MTDOCA",                   "MTDOCA(#Lambda_{b})",                                      "mm",  'F' );

    # Read input data
    if gSystem.AccessPathName( sigfname ) != 0: print "Can not find %s" % sigfname
    if gSystem.AccessPathName( bkgfname ) != 0: print "Can not find %s" % bkgfname

    inputSig = TFile.Open( sigfname )
    inputBkg = TFile.Open( bkgfname )

    # Get the signal and background trees for training
    signal      = inputSig.Get( treeNameSig )
    background  = inputBkg.Get( treeNameBkg )
    
    # Global event weights (see below for setting event-wise weights)
    signalWeight     = 1.0
    backgroundWeight = 1.0

    # ====== register trees ====================================================
    #
    # the following method is the prefered one:
    # you can add an arbitrary number of signal or background trees
    dataloader.AddSignalTree    ( signal,     signalWeight     )
    dataloader.AddBackgroundTree( background, backgroundWeight )

    # To give different trees for training and testing, do as follows:
    #    dataloader.AddSignalTree( signalTrainingTree, signalTrainWeight, "Training" )
    #    dataloader.AddSignalTree( signalTestTree,     signalTestWeight,  "Test" )
    
    # Use the following code instead of the above two or four lines to add signal and background 
    # training and test events "by hand"
    # NOTE that in this case one should not give expressions (such as "var1+var2") in the input 
    #      variable definition, but simply compute the expression before adding the event
    #
    #    # --- begin ----------------------------------------------------------
    #    
    # ... *** please lookup code in TMVA/macros/TMVAClassification.C ***
    #    
    #    # --- end ------------------------------------------------------------
    #
    # ====== end of register trees ==============================================    
            
    # Set individual event weights (the variables must exist in the original TTree)
    #    for signal    : dataloader.SetSignalWeightExpression    ("weight1*weight2");
    #    for background: dataloader.SetBackgroundWeightExpression("weight1*weight2");
#dataloader.SetBackgroundWeightExpression( "weight" )

    # Apply additional cuts on the signal and background sample. 
    # example for cut: mycut = TCut( "abs(var1)<0.5 && abs(var2-0.5)<1" )
    mycutSig = TCut( "" ) 
    mycutBkg = TCut( "" ) 
    
    # Here, the relevant variables are copied over in new, slim trees that are
    # used for TMVA training and testing
    # "SplitMode=Random" means that the input events are randomly shuffled before
    # splitting them into training and test samples
    dataloader.PrepareTrainingAndTestTree( mycutSig, mycutBkg,
                                        "nTrain_Signal=0:nTrain_Background=0:SplitMode=Random:NormMode=NumEvents:!V" )

    # --------------------------------------------------------------------------------------------------

    # ---- Book MVA methods
    #
    # please lookup the various method configuration options in the corresponding cxx files, eg:
    # src/MethoCuts.cxx, etc, or here: http://tmva.sourceforge.net/optionRef.html
    # it is possible to preset ranges in the option string in which the cut optimisation should be done:
    # "...:CutRangeMin[2]=-1:CutRangeMax[2]=1"...", where [2] is the third input variable

    # Cut optimisation
    if "Cuts" in mlist:
        factory.BookMethod( dataloader, TMVA.Types.kCuts, "Cuts",
                            "!H:!V:FitMethod=MC:EffSel:SampleSize=200000:VarProp=FSmart" )

    if "CutsD" in mlist:
        factory.BookMethod( dataloader, TMVA.Types.kCuts, "CutsD",
                            "!H:!V:FitMethod=MC:EffSel:SampleSize=200000:VarProp=FSmart:VarTransform=Decorrelate" )

    if "CutsPCA" in mlist:
        factory.BookMethod( dataloader, TMVA.Types.kCuts, "CutsPCA",
                            "!H:!V:FitMethod=MC:EffSel:SampleSize=200000:VarProp=FSmart:VarTransform=PCA" )

    if "CutsGA" in mlist:
        factory.BookMethod( dataloader, TMVA.Types.kCuts, "CutsGA",
                            "H:!V:FitMethod=GA:CutRangeMin[0]=-10:CutRangeMax[0]=10:VarProp[1]=FMax:EffSel:Steps=30:Cycles=3:PopSize=400:SC_steps=10:SC_rate=5:SC_factor=0.95" )

    if "CutsSA" in mlist:
        factory.BookMethod( dataloader, TMVA.Types.kCuts, "CutsSA",
                            "!H:!V:FitMethod=SA:EffSel:MaxCalls=150000:KernelTemp=IncAdaptive:InitialTemp=1e+6:MinTemp=1e-6:Eps=1e-10:UseDefaultScale" )

    # Likelihood ("naive Bayes estimator")
    if "Likelihood" in mlist:
        factory.BookMethod( dataloader, TMVA.Types.kLikelihood, "Likelihood",
                            "H:!V:!TransformOutput:PDFInterpol=Spline2:NSmoothSig[0]=20:NSmoothBkg[0]=20:NSmoothBkg[1]=10:NSmooth=1:NAvEvtPerBin=50" )

    # Decorrelated likelihood
    if "LikelihoodD" in mlist:
        factory.BookMethod( dataloader, TMVA.Types.kLikelihood, "LikelihoodD",
                            "!H:!V:TransformOutput:PDFInterpol=Spline2:NSmoothSig[0]=20:NSmoothBkg[0]=20:NSmooth=5:NAvEvtPerBin=50:VarTransform=Decorrelate" )

    # PCA-transformed likelihood
    if "LikelihoodPCA" in mlist:
        factory.BookMethod( dataloader, TMVA.Types.kLikelihood, "LikelihoodPCA",
                            "!H:!V:!TransformOutput:PDFInterpol=Spline2:NSmoothSig[0]=20:NSmoothBkg[0]=20:NSmooth=5:NAvEvtPerBin=50:VarTransform=PCA" ) 

    # Use a kernel density estimator to approximate the PDFs
    if "LikelihoodKDE" in mlist:
        factory.BookMethod( dataloader, TMVA.Types.kLikelihood, "LikelihoodKDE",
                            "!H:!V:!TransformOutput:PDFInterpol=KDE:KDEtype=Gauss:KDEiter=Adaptive:KDEFineFactor=0.3:KDEborder=None:NAvEvtPerBin=50" ) 

    # Use a variable-dependent mix of splines and kernel density estimator
    if "LikelihoodMIX" in mlist:
        factory.BookMethod( dataloader, TMVA.Types.kLikelihood, "LikelihoodMIX",
                            "!H:!V:!TransformOutput:PDFInterpolSig[0]=KDE:PDFInterpolBkg[0]=KDE:PDFInterpolSig[1]=KDE:PDFInterpolBkg[1]=KDE:PDFInterpolSig[2]=Spline2:PDFInterpolBkg[2]=Spline2:PDFInterpolSig[3]=Spline2:PDFInterpolBkg[3]=Spline2:KDEtype=Gauss:KDEiter=Nonadaptive:KDEborder=None:NAvEvtPerBin=50" ) 

    # Test the multi-dimensional probability density estimator
    # here are the options strings for the MinMax and RMS methods, respectively:
    #      "!H:!V:VolumeRangeMode=MinMax:DeltaFrac=0.2:KernelEstimator=Gauss:GaussSigma=0.3" );
    #      "!H:!V:VolumeRangeMode=RMS:DeltaFrac=3:KernelEstimator=Gauss:GaussSigma=0.3" );
    if "PDERS" in mlist:
        factory.BookMethod( dataloader, TMVA.Types.kPDERS, "PDERS",
                            "!H:!V:NormTree=T:VolumeRangeMode=Adaptive:KernelEstimator=Gauss:GaussSigma=0.3:NEventsMin=400:NEventsMax=600" )

    if "PDERSD" in mlist:
        factory.BookMethod( dataloader, TMVA.Types.kPDERS, "PDERSD",
                            "!H:!V:VolumeRangeMode=Adaptive:KernelEstimator=Gauss:GaussSigma=0.3:NEventsMin=400:NEventsMax=600:VarTransform=Decorrelate" )

    if "PDERSPCA" in mlist:
        factory.BookMethod( dataloader, TMVA.Types.kPDERS, "PDERSPCA",
                             "!H:!V:VolumeRangeMode=Adaptive:KernelEstimator=Gauss:GaussSigma=0.3:NEventsMin=400:NEventsMax=600:VarTransform=PCA" )

   # Multi-dimensional likelihood estimator using self-adapting phase-space binning
    if "PDEFoam" in mlist:
        factory.BookMethod( dataloader, TMVA.Types.kPDEFoam, "PDEFoam",
                            "!H:!V:SigBgSeparate=F:TailCut=0.001:VolFrac=0.0666:nActiveCells=500:nSampl=2000:nBin=5:Nmin=100:Kernel=None:Compress=T" )

    if "PDEFoamBoost" in mlist:
        factory.BookMethod( dataloader, TMVA.Types.kPDEFoam, "PDEFoamBoost",
                            "!H:!V:Boost_Num=30:Boost_Transform=linear:SigBgSeparate=F:MaxDepth=4:UseYesNoCell=T:DTLogic=MisClassificationError:FillFoamWithOrigWeights=F:TailCut=0:nActiveCells=500:nBin=20:Nmin=400:Kernel=None:Compress=T" )

    # K-Nearest Neighbour classifier (KNN)
    if "KNN" in mlist:
        factory.BookMethod( dataloader, TMVA.Types.kKNN, "KNN",
                            "H:nkNN=20:ScaleFrac=0.8:SigmaFact=1.0:Kernel=Gaus:UseKernel=F:UseWeight=T:!Trim" )

    # H-Matrix (chi2-squared) method
    if "HMatrix" in mlist:
        factory.BookMethod( dataloader, TMVA.Types.kHMatrix, "HMatrix", "!H:!V" )

    # Linear discriminant (same as Fisher discriminant)
    if "LD" in mlist:
        factory.BookMethod( dataloader, TMVA.Types.kLD, "LD", "H:!V:VarTransform=None:CreateMVAPdfs:PDFInterpolMVAPdf=Spline2:NbinsMVAPdf=50:NsmoothMVAPdf=10" )

    # Fisher discriminant (same as LD)
    if "Fisher" in mlist:
        factory.BookMethod( dataloader, TMVA.Types.kFisher, "Fisher", "H:!V:Fisher:CreateMVAPdfs:PDFInterpolMVAPdf=Spline2:NbinsMVAPdf=50:NsmoothMVAPdf=10" )

    # Fisher with Gauss-transformed input variables
    if "FisherG" in mlist:
        factory.BookMethod( dataloader, TMVA.Types.kFisher, "FisherG", "H:!V:VarTransform=Gauss" )

    # Composite classifier: ensemble (tree) of boosted Fisher classifiers
    if "BoostedFisher" in mlist:
        factory.BookMethod( dataloader, TMVA.Types.kFisher, "BoostedFisher", 
                            "H:!V:Boost_Num=20:Boost_Transform=log:Boost_Type=AdaBoost:Boost_AdaBoostBeta=0.2" )

    # Function discrimination analysis (FDA) -- test of various fitters - the recommended one is Minuit (or GA or SA)
    if "FDA_MC" in mlist:
        factory.BookMethod( dataloader, TMVA.Types.kFDA, "FDA_MC",
                            "H:!V:Formula=(0)+(1)*x0+(2)*x1+(3)*x2+(4)*x3:ParRanges=(-1,1)(-10,10);(-10,10);(-10,10);(-10,10):FitMethod=MC:SampleSize=100000:Sigma=0.1" );

    if "FDA_GA" in mlist:
        factory.BookMethod( dataloader, TMVA.Types.kFDA, "FDA_GA",
                            "H:!V:Formula=(0)+(1)*x0+(2)*x1+(3)*x2+(4)*x3:ParRanges=(-1,1);(-10,10);(-10,10);(-10,10);(-10,10):FitMethod=GA:PopSize=300:Cycles=3:Steps=20:Trim=True:SaveBestGen=1" );

    if "FDA_SA" in mlist:
        factory.BookMethod( dataloader, TMVA.Types.kFDA, "FDA_SA",
                            "H:!V:Formula=(0)+(1)*x0+(2)*x1+(3)*x2+(4)*x3:ParRanges=(-1,1)(-10,10);(-10,10);(-10,10);(-10,10):FitMethod=SA:MaxCalls=15000:KernelTemp=IncAdaptive:InitialTemp=1e+6:MinTemp=1e-6:Eps=1e-10:UseDefaultScale" );

    if "FDA_MT" in mlist:
        factory.BookMethod( dataloader, TMVA.Types.kFDA, "FDA_MT",
                            "H:!V:Formula=(0)+(1)*x0+(2)*x1+(3)*x2+(4)*x3:ParRanges=(-1,1)(-10,10);(-10,10);(-10,10);(-10,10):FitMethod=MINUIT:ErrorLevel=1:PrintLevel=-1:FitStrategy=2:UseImprove:UseMinos:SetBatch" );

    if "FDA_GAMT" in mlist:
        factory.BookMethod( dataloader, TMVA.Types.kFDA, "FDA_GAMT",
                            "H:!V:Formula=(0)+(1)*x0+(2)*x1+(3)*x2+(4)*x3:ParRanges=(-1,1)(-10,10);(-10,10);(-10,10);(-10,10):FitMethod=GA:Converger=MINUIT:ErrorLevel=1:PrintLevel=-1:FitStrategy=0:!UseImprove:!UseMinos:SetBatch:Cycles=1:PopSize=5:Steps=5:Trim" );

    if "FDA_MCMT" in mlist:
        factory.BookMethod( dataloader, TMVA.Types.kFDA, "FDA_MCMT",
                            "H:!V:Formula=(0)+(1)*x0+(2)*x1+(3)*x2+(4)*x3:ParRanges=(-1,1)(-10,10);(-10,10);(-10,10);(-10,10):FitMethod=MC:Converger=MINUIT:ErrorLevel=1:PrintLevel=-1:FitStrategy=0:!UseImprove:!UseMinos:SetBatch:SampleSize=20" );

    # TMVA ANN: MLP (recommended ANN) -- all ANNs in TMVA are Multilayer Perceptrons
    if "MLP" in mlist:
        factory.BookMethod( dataloader, TMVA.Types.kMLP, "MLP", "H:!V:NeuronType=tanh:VarTransform=N:NCycles=600:HiddenLayers=N+5:TestRate=5:!UseRegulator" )

    if "MLPBFGS" in mlist:
        factory.BookMethod( dataloader, TMVA.Types.kMLP, "MLPBFGS", "H:!V:NeuronType=tanh:VarTransform=N:NCycles=600:HiddenLayers=N+5:TestRate=5:TrainingMethod=BFGS:!UseRegulator" )

    if "MLPBNN" in mlist:
        factory.BookMethod( dataloader, TMVA.Types.kMLP, "MLPBNN", "H:!V:NeuronType=tanh:VarTransform=N:NCycles=600:HiddenLayers=N+5:TestRate=5:TrainingMethod=BFGS:UseRegulator" ) # BFGS training with bayesian regulators

    # CF(Clermont-Ferrand)ANN
    if "CFMlpANN" in mlist:
        factory.BookMethod( dataloader, TMVA.Types.kCFMlpANN, "CFMlpANN", "!H:!V:NCycles=2000:HiddenLayers=N+1,N"  ) # n_cycles:#nodes:#nodes:...  

    # Tmlp(Root)ANN
    if "TMlpANN" in mlist:
        factory.BookMethod( dataloader, TMVA.Types.kTMlpANN, "TMlpANN", "!H:!V:NCycles=200:HiddenLayers=N+1,N:LearningMethod=BFGS:ValidationFraction=0.3"  ) # n_cycles:#nodes:#nodes:...

    # Support Vector Machine
    if "SVM" in mlist:
        factory.BookMethod( dataloader, TMVA.Types.kSVM, "SVM", "Gamma=0.25:Tol=0.001:VarTransform=Norm" )

    # Boosted Decision Trees
    if "BDTG" in mlist:
        factory.BookMethod( dataloader, TMVA.Types.kBDT, "BDTG",
                            "!H:!V:NTrees=1000:BoostType=Grad:Shrinkage=0.30:UseBaggedGrad:GradBaggingFraction=0.6:SeparationType=GiniIndex:nCuts=20:NNodesMax=5" )

    if "BDT" in mlist:
        factory.BookMethod( dataloader, TMVA.Types.kBDT, "BDT",
                            "!H:!V:NTrees=850:MinNodeSize=2.5%:MaxDepth=3:BoostType=AdaBoost:AdaBoostBeta=0.5:SeparationType=GiniIndex:nCuts=20:PruneMethod=NoPruning" )

    if "BDTB" in mlist:
        factory.BookMethod( dataloader, TMVA.Types.kBDT, "BDTB",
                            "!H:!V:NTrees=400:BoostType=Bagging:SeparationType=GiniIndex:nCuts=20:PruneMethod=NoPruning" )

    if "BDTD" in mlist:
        factory.BookMethod( dataloader, TMVA.Types.kBDT, "BDTD",
                            "!H:!V:NTrees=400:nEventsMin=400:MaxDepth=3:BoostType=AdaBoost:SeparationType=GiniIndex:nCuts=20:PruneMethod=NoPruning:VarTransform=Decorrelate" )

    # RuleFit -- TMVA implementation of Friedman's method
    if "RuleFit" in mlist:
        factory.BookMethod( dataloader, TMVA.Types.kRuleFit, "RuleFit",
                            "H:!V:RuleFitModule=RFTMVA:Model=ModRuleLinear:MinImp=0.001:RuleMinDist=0.001:NTrees=20:fEventsMin=0.01:fEventsMax=0.5:GDTau=-1.0:GDTauPrec=0.01:GDStep=0.01:GDNSteps=10000:GDErrScale=1.02" )

    # --------------------------------------------------------------------------------------------------
            
    # ---- Now you can tell the factory to train, test, and evaluate the MVAs. 

    # Train MVAs
    factory.TrainAllMethods()
    
    # Test MVAs
    factory.TestAllMethods()
    
    # Evaluate MVAs
    factory.EvaluateAllMethods()    
    
    # Save the output.
    outputFile.Close()
    
    print "=== wrote root file %s\n" % outfname
    print "=== TMVAClassification is done!\n"
    
    # open the GUI for the result macros
    if not gROOT.IsBatch(): TMVA.TMVAGui( outfname )
    
    # keep the ROOT thread running
    #gApplication.Run()

# ----------------------------------------------------------

if __name__ == "__main__":
    try:
        # retrive command line options
        shortopts  = "m:s:b:t:o:c:vh?"
        longopts   = ["methods=", "sigfile=", "bkgfile=", "inputtrees=", "outputfile=", "channel=", "verbose", "help", "usage"]
        opts, args = getopt.getopt( sys.argv[1:], shortopts, longopts )
        
    except getopt.GetoptError:
        # print help information and exit:
        print "ERROR: unknown options in argument %s" % sys.argv[1:]
        usage()
        sys.exit(1)

    DEFAULT_OPTNAME  = "BDTOutput"
    DEFAULT_SIGFNAME = "Lb2L0Gamma_presel_sig_mcmatch_LL_presel.root"
    DEFAULT_BKGFNAME = "Lb2L0Gamma_presel_bkg_LL_presel.root"
    DEFAULT_CHANNEL  = "LL"
    DEFAULT_TREESIG  = "Lb2L0GammaTuple/DecayTree"
    DEFAULT_TREEBKG  = "Lb2L0GammaTuple/DecayTree"
    DEFAULT_METHODS  = "BDT"

    sigfname    = DEFAULT_SIGFNAME
    bkgfname    = DEFAULT_BKGFNAME
    treeNameSig = DEFAULT_TREESIG
    treeNameBkg = DEFAULT_TREEBKG
    optname     = DEFAULT_OPTNAME
    channel     = DEFAULT_CHANNEL
    methods     = DEFAULT_METHODS
    verbose     = False
    for o, a in opts:
        if o in ("-?", "-h", "--help", "--usage"):
            usage()
            sys.exit(0)
        elif o in ("-m", "--methods"):
            methods  = a
        elif o in ("-s", "--sigfile"):
            sigfname = a
        elif o in ("-b", "--bkgfile"):
            bkgfname = a
        elif o in ("-o", "--optname"):
            optname = a
        elif o in ("-c", "--channel"):
            channel  = a
        elif o in ("-t", "--inputtrees"):
            trees = a
        elif o in ("-v", "--verbose"):
            verbose = True

    TMVAClassification(methods, sigfname, bkgfname, optname, channel, trees=trees, verbose=False)

#EOF