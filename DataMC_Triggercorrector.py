import sys
import os
tools = os.path.expandvars('$ANA/AnaTools/')
sys.path.insert(0,tools)
sys.path.insert(0,tools+'Tools/')

import argparse

import ROOT
from ROOT import gROOT, gStyle, TFile, TLorentzVector, TVector3, TRotation, TLorentzRotation, TMath, TH1D, TCanvas, TH2D, TObject, TF1, TH1F, TLegend, kTRUE
import itertools as it
import csv
import sqlite3
import math as m
from datetime import datetime
import time
import os.path
import warnings
import seaborn as sns
from numpy import array
from root_numpy import tree2array, array2root
import numpy as np
import pandas as pd
from root_pandas import read_root,to_root
import matplotlib.pyplot as plt
#from TagAndProbe_e import TagAndProbe_e
#from Presels import GenericPresel, VetoesPresel, TriggerPreselE, PIDPresel, TighterKst0Presel

from Adaptive_binning_pro import AdaptiveBinning1D
from Vocabulary import  jobsDict, type_list, channel_list, yearsRI, yearsRII, mag_list
from Tools0 import GetJob, GetList_for_ReducedChain, listdirs, listfiles

from TagAndProbe_kFold_pro import TagAndProbe_L0E, TagAndProbe_L0H, TagAndProbe_L0M, TagAndProbe_L0TIS, TagAndProbe_HLT, TagAndProbe_L0E_onlynumerator, TagAndProbe_L0H_onlynumerator, TagAndProbe_L0M_onlynumerator, TagAndProbe_L0TIS_onlynumerator, TagAndProbe_HLT_onlynumerator
from reweighting import reweighting
from Plotdf import PlotDF
import pdb
from Correct_MC_with_Data_pro import Correct_MC_with_Data, Reweighting

from UtilsTriggerCalib import CalibSelection_M,CalibSelection_E

if __name__ == "__main__" :


    '''
    Author: Michele Atzeni
    Date: June 1st, 2017

    Description:
    Takes the TightKst0 dataframes for MC (TM) and data (it should be BDT reweighted) and plots the histograms for the trigger efficiencies.

    How to run it:
    Choose the year of the desired MC/data correction of the trigger and execute!
    > python L0TriggerDataMC.py [-y 11] [--test]



  Important:                                                                                                                        
    
    Selection for MC (after TightKst0 preselection and Truth Matching):                                                               
    B02Kst0Jpsi2ee-> + B_PVandJpsiDTF_B_M in [5150., 5900.] MeV/c^2                                                                   
                     +         q^2        in [6., 11.]*10^5 MeV^2/c^4                                                                 
                                                                                                                                      
    B02Kst0Jpsi2mm-> + B_PVandJpsiDTF_B_M in [5150., 5900.] MeV/c^2                                                                   
                     +         Jpsi_M     in [2996.9., 3196.9] MeV/c^2                                                                
                                                                                                                                      
                                                                                                                                      
    #Selection for Data (after TightKst0 preselection):                                                                               
    #B02Kst0Jpsi2ee-> + B_PVandJpsiDTF_B_M in [5150., 5900.] MeV/c^2                                                                  
    #                 +         q^2        in [6., 11.]*10^5 MeV^2/c^4                                                                
    #                                                                                                                                 
    #B02Kst0Jpsi2mm-> + B_PVandJpsiDTF_B_M in [5150., 5900.] MeV/c^2                                                                  
    #                 +         Jpsi_M     in [2996.9., 3196.9] MeV/c^2                            
    '''
        
    parser = argparse.ArgumentParser(description = 'Configuration of the parameters for the SplitAndMerge')
    
    parser.add_argument("-y", "--year" , dest="y"  , required=False, help="Choose the year for the RI-RII", choices=yearsRI+yearsRII+['RunI','RunII'], default = '11')
    parser.add_argument("-l", "--leptons" , dest="l"  , required=False, help="Choose the leptons in the final state", choices= ['ee','mm'], default = 'ee')
    parser.add_argument("-v", "--versionq2" , dest="version"  , required=False, help="", choices= ['All','q2','q2_PVandBDTF'], default = 'All')
    parser.add_argument("-low", "--lowq2" , dest="low_val"  , required=False, help="", choices= ['6','7'], default = '6')
    parser.add_argument("-u", "--user" , dest="user"  , required=False, help="", choices= ['M','F'], default = 'F')

    parser.add_argument("--test", action="store_true", help="Do a test plot")
    parser.add_argument("--version", action="store_true", help="q2version")
    parser.add_argument("--TM", action="store_true", help="Truth Matched sample")

    parser.add_argument("--VERB", action="store_true", help="VERBOSE")

    

    args = parser.parse_args()
    
    # Parameters and configuration                                                                                                                             
    year= args.y
    leptons= args.l
    version= args.version
    low_val= args.low_val
    user= args.user
    TM = args.TM
    test = args.test

    VERB = args.VERB


    if (user == 'F'):
        #directory = '/disk/data12/lhcb/flionett/LHCb/Analysis/AnalysisResults/ControlMode/Clean/6To11All/EffWithSymOption/'
        directory = '/disk/data12/lhcb/flionett/LHCb/Analysis/AnalysisResults/ControlMode/Clean/6To11All/EffWithSymOption/'
        
        #version = 'All'
        #low_val = '6.'
        
        dir_L0HLT = './test/'
        print(listfiles(directory))
        for ifile in listfiles(directory):
            #ifile = 'B02Kst0Jpsi2ee-RunI-MC-central-L0E-0.95-HOP0-BkgCat0-Ctl0-51501-Mass4500To5800-Eff.h5'
            ifile_spltd = ifile.split("-")
            
            if("MC" in ifile_spltd):


                print "Opening the file: %s ...."%ifile
                
                store = pd.HDFStore(directory+ifile)
                
                if(store.keys()):
                    print "The DataFrames available are: ", store.keys()
                
                    df = store[ifile_spltd[0]]
                    
                    if( "ee" in ifile_spltd[0]):
                        df =CalibSelection_E(df, TM, version, low_val)
                        
                    elif ("mm" in ifile_spltd[0]):
                        df =CalibSelection_M(df, TM)

                    else:
                        print "Not clear what to do..."
                        exit()

                    df11 = df[df.Year == 11]
                    df12 = df[df.Year == 12]
                    dfRW11 = Reweighting(df11, ifile_spltd[4],'11', TM, version, low_val)
                    dfRW12 = Reweighting(df12, ifile_spltd[4],'12',TM, version, low_val)

                    dfRW = pd.concat([dfRW11,dfRW12],axis = 0)

                    dfRW.to_hdf(dir_L0HLT+ifile.split(".")[0]+"_L0HLTrw.h5", format='table', key="dfL0HLT")

                    store.close()
                    #exit()
                else:
                    print '========================================================================================'
                    print"WARNING! The file {} is empty, please verify that this is not a mistake.".format(ifile)
                    print '========================================================================================'

            else:
                print "Jumping..."

    elif (user == 'M'):
        #directory = '/disk/data12/lhcb/flionett/LHCb/Analysis/AnalysisResults/ControlMode/Clean/6To11All/EffWithSymOption/'
        directory = '/home/hep/matzeni/gangadir/Analysis/DFs/MC/B02Kst0Jpsi2ee/DFs/WithHOP/ParallelSplitOutput/L0HLT/'
        
        #version = 'All'
        #low_val = '6.'
        
        dir_L0HLT = './test/'
        print(listfiles(directory))
        for ifile in listfiles(directory):
            #ifile = 'B02Kst0Jpsi2ee-RunI-MC-central-L0E-0.95-HOP0-BkgCat0-Ctl0-51501-Mass4500To5800-Eff.h5'
            ifile_spltd = ifile.split("-")
            
            
            print "Opening the file: %s ...."%ifile
            
            store = pd.HDFStore(directory+ifile)
                
            if(store.keys()):
                print "The DataFrames available are: ", store.keys()
                
                df = store["dfTightKst0_noTrig"]
                
                if( "ee" in ifile_spltd[0]):
                    #df =CalibSelection_E(df, TM, version, low_val)

                    dfL0E = df[(df.L2_L0ElectronDecision_TOS == 1) | (df.L1_L0ElectronDecision_TOS == 1)]
                    dfL0H = df[(df.Kstar_L0HadronDecision_TOS==1) & (df.L2_L0ElectronDecision_TOS == 0) & (df.L1_L0ElectronDecision_TOS == 0)]
                    dfL0TIS = df[(df.B_L0Global_TIS==1) & (df.Kstar_L0HadronDecision_TOS==0) & (df.L2_L0ElectronDecision_TOS == 0) & (df.L1_L0ElectronDecision_TOS == 0)]

                elif ("mm" in ifile_spltd[0]):
                    #CalibSelection_M(df, TM)
                    
                    dfL0E = df[(df.L2_L0MuonDecision_TOS == 1) | (df.L1_L0MuonDecision_TOS == 1)]
                    dfL0H = df[(df.Kstar_L0HadronDecision_TOS==1) & (df.L2_L0MuonDecision_TOS == 0) & (df.L1_L0MuonDecision_TOS == 0)]
                    dfL0TIS = df[(df.B_L0Global_TIS==1) & (df.Kstar_L0HadronDecision_TOS==0) & (df.L2_L0MuonDecision_TOS == 0) & (df.L1_L0MuonDecision_TOS == 0)]


                else:
                    print "Not clear what to do..."
                    exit()
                '''
                df11 = df[df.Year == 11]
                df12 = df[df.Year == 12]
                
                dfRW11 = Reweighting(df11, ifile_spltd,'11', TM, version, low_val)
                dfRW12 = Reweighting(df12, ifile_spltd,'12',TM, version, low_val)
                '''
                dfRWL0E = Reweighting(dfL0E, 'L0E','11', TM, version, low_val)
                dfRWL0H = Reweighting(dfL0H, 'L0H','11',TM, version, low_val)
                dfRWL0TIS = Reweighting(dfL0TIS, 'L0TIS','11',TM, version, low_val)
                dfRW = pd.concat([dfRWL0E,dfRWL0H, dfL0TIS],axis = 0)
                
                dfRW.to_hdf(dir_L0HLT+ifile.split(".")[0]+"_L0HLTrw.h5", format='table', key="dfL0HLT")
                
                store.close()
                #exit()
            else:
                print '========================================================================================'
                print"WARNING! The file {} is empty, please verify that this is not a mistake.".format(ifile)
                print '========================================================================================'


