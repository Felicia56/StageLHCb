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
import time   # time accounting
import getopt # command line parser
import os, pickle
from ROOT import gStyle,gROOT,TLegend,TPad,TCanvas,TH1F,TH2F,TFile,TDirectory,TCut,THStack,TLatex,TTree
from /Stage/setup.py import simpleFit


#____________________________FUNCTION DEFINITION______________________________

BDTG_val_list = [-1,-0.8,-0.6,-0.4,-0.2,0,0.2,0.4,0.6,0.8,1]

#____________________________FUNCTION DEFINITION______________________________



#____________________________________MAIN_____________________________________


if __name__ == '__main__':
    import datetime
    import sys
    
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("file", action="store", type=str)
    parser.add_argument("-t", "--tree", default="DecayTree", action="store", type=str)
    #parser.add_argument("-m", "--mean", default=5620., action="store", type=float)
    parser.add_argument("-n", "--xmin", default=6120., action="store", type=float)
    parser.add_argument("-x", "--xmax", default=9120., action="store", type=float)
    #parser.add_argument("-c", "--cuts", default="", action="store", type=str)
    parser.add_argument("-v", "--version", default="1st",action="store", type=str)
    args = parser.parse_args()

    newDir = '/sps/lhcb/volle/FitResults/{}'.format(args.version)
    if not os.path.isdir(newDir):
        os.makedirs(newDir)
    FitDir = "/users/LHCb/volle/Stage/simpleFit.py"

    for i in BDTG_val_list:
        logfile = '{}/Output.log'.format(newDir)
        cmd = 'python {} -v {} -x {} -c {}'.format(FitDir,args.version,args.xmax, i)

    #!!!!!!!!!!Muss neu!!!!!!!!!!!
    from setup import partition
    sbatch_cmd = 'sbatch -p {partition} --mem-per-cpu 1024 -n 1 -t 300 -o {0} --wrap "{1}"'.format(logfile,cmd,partition=partition)
    os.system(sbatch_cmd)


    


