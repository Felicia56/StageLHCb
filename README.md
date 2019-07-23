# StageLHCb

Workflow and explication:
(sometimes special environnement needed: lb-run Urania/v7r0 bash --norc)

1) Preselection

PlotPreselection.py -> plotting for each cut step data and MC
PlotHistograms.py -> plotting after preselection for left side & right side massband for data, compared to the MC
PlotBmass.py -> plotting mass of Lb for data after all preselection

2) TMVA for reducing combinatorial bkg

copytree.C -> copy tree & apply preselection (so no error in TMVAClassification.py)
TMVAClassification.py -> optimize for best MVA method (used BDTG)
   - looking at TMVA method via TMVA::TMVAGui("output-root-file")
   - adding BDTG variable to the DecayTree for data & MC!!
     (https://github.com/pseyfert/tmva-branch-adder: do (I) and follow instructions
     (I) git clone https://github.com/pseyfert/tmva-branch-adder.git )
PlotCorrBDTGMass.py -> look for each BDTG the correlation between Lb-mass & mean BDTG value

3) Selection best BDTG cut value

FigureOfMerit.py -> fitting left and rigth sideband in data with an exponential -> Nbkg
& extracting efficiency of MC & calculates Nsig with BR-value of thesis
-> significance figure of merit & punzi figure of merit (plotting everything)
needs script simpleFit.py at same place

4) Fitting Lb-mass distribution

FitMCshape.py -> fitting naively gaussian shape to MC - is not a good fit here!
FitBMass2.py -> fitting data by using for the signal a gaussian and for the background an exponential and a gaussian & producing sweights
FitBMass.py -> fitting background by an exponential and signal by an gaussian core with two exponential tails - the exponential tails are extracted from fit to MC shape and fixed (BifurcatedCB.cxx & BifurcatedCB.h needed) & producing sweights 

5) Splot

Splot.py : produced root file from FitBMass.py needed 
-> sweighted plot of L* and cos(theta_L*) in CoM-frame of L* (measured to direction of Lb)

Bonne courage!! :)
