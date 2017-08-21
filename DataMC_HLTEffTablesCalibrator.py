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
#from Presels import GenericPresel, VetoesPresel, TriggePreselE, PIDPresel, TighterKst0Presel

from Adaptive_binning_pro import AdaptiveBinning1D




from TagAndProbe_kFold import TagAndProbe_L0E, TagAndProbe_L0H, TagAndProbe_L0M, TagAndProbe_L0TIS, TagAndProbe_HLT, TagAndProbe_L0E_onlynumerator, TagAndProbe_L0H_onlynumerator, TagAndProbe_L0M_onlynumerator, TagAndProbe_L0TIS_onlynumerator, TagAndProbe_HLT_onlynumerator
from reweighting import reweighting
from Plotdf import PlotDF
import pdb
from Correct_MC_with_Data import Correct_MC_with_Data_E_maxET, Correct_MC_with_Data_M_maxPT, Correct_MC_with_Data_KPi, Correct_MC_with_Data_TIS, Correct_MC_with_Data_E_mixedmaxET, Correct_MC_with_Data_M_mixedmaxPT, Correct_MC_with_Data_HLT


from UtilsTriggerCalib import Open_files_for_TriggerCalibration,PreselB02Kst0Jpsi2eeDF,CalibSelection_M,CalibSelection_E
from TablesCalibrators import CalibrationTables_L0, CalibrationTables_HLT


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
    
    parser.add_argument("-y", "--year" , dest="y"  , required=False, help="Choose the year for the RI-RII", choices=['11','12','15','16','RunI','RunII'], default = 'RunI')
    parser.add_argument("-l", "--leptons" , dest="l"  , required=False, help="Choose the leptons in the final state", choices= ['ee','mm'], default = 'ee')
    parser.add_argument("-v", "--versionq2" , dest="version"  , required=False, help="", choices= ['All','q2','q2_PVandBDTF'], default = 'All')
    parser.add_argument("-low", "--lowq2" , dest="low_val"  , required=False, help="", choices= ['6','7'], default = '6')
    parser.add_argument("-u", "--user" , dest="user"  , required=False, help="", choices= ['M','F'], default = 'M')

    parser.add_argument("--test", action="store_true", help="Do a test plot")
    parser.add_argument("--version", action="store_true", help="q2version")
    parser.add_argument("--TM", action="store_true", help="TruthMatchedMC")

    parser.add_argument("--plot", action="store_true", help="Do a test plot")
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
    plot = args.plot
    VERB = args.VERB

    TM = True
    #version = 'All'
    #low_val = '6'
    channelData = 'B02Kst0{}'.format(leptons)
    channelMC = 'B02Kst0Jpsi2{}'.format(leptons)
    

    if(user == "M"):
        directoryMC = '/home/hep/matzeni/gangadir/Analysis/DFs/MC/{}/DFs/WithHOP/ParallelSplitOutput/L0HLT/'.format(channelMC)
        directoryData = '/home/hep/matzeni/gangadir/Analysis/DFs/Data/{}/DFs/WithHOP/ParallelSplitOutput/L0HLT/'.format(channelData)

        from Vocabulary import  jobsDict, type_list, channel_list, yearsRI, yearsRII, mag_list
        from Tools0 import GetJob, GetList_for_ReducedChain, listdirs, listfiles
    else:
        directoryMC = '/disk/data12/lhcb/flionett/LHCb/Analysis/Ntuples/MC/{}/DFs/WithHOP/ParallelSplitOutput/Meerkat/'.format( channelMC)
        directoryData = '/disk/data12/lhcb/flionett/LHCb/Analysis/Ntuples/Data/{}/DFs/WithHOP/ParallelSplitOutput/L0HLT/'.format( channelData)
        from voc import  jobsDict, type_list, channel_list, yearsRI, yearsRII, mag_list
        from tools import GetJob, GetList_for_ReducedChain, listdirs, listfiles

    dfMC = Open_files_for_TriggerCalibration(directoryMC, jobsDict,'MC', channelMC, year, TM, version, low_val, test)
    Eff_tables_MC, Eff_root_MC = CalibrationTables_HLT(dfMC,"MC", channelMC, year, leptons, VERB)
    
    del dfMC
    #######
    os.system('mkdir -p EffTable')

    if(TM):
        Tag_sample = 'q2{}_lw{}_TM'.format(version, low_val)
    else:
        Tag_sample = 'q2{}_lw{}'.format(version, low_val)

    import pickle  
    print "Writing the efficiency tables to EffTable/EfficiencyTables_Calib_HLT_{}_{}_{}-{}.pkl".format(channelMC, year, "MC", Tag_sample)
    pickle.dump(Eff_tables_MC, open('./EffTable/EfficiencyTables_Calib_HLT_{}_{}_{}-{}.pkl'.format(channelMC, year, "MC", Tag_sample), 'wb'))
    
    
        
    print "Saving the Calibration histograms in EffTable/EffHisto_Calib_HLT_{}_{}_{}-{}.root".format(channelMC, year, "MC", Tag_sample)
    file_root_MC = TFile("EffTable/EffHisto_Calib_HLT_{}_{}_{}-{}.root".format(channelMC, year, "MC", Tag_sample),"RECREATE")

    map(lambda x:x.Write(), Eff_root_MC)
    file_root_MC.Close()
    

    dfData = Open_files_for_TriggerCalibration(directoryData, jobsDict, 'Data', channelData, year, False, version, low_val, test)
    Eff_tables_Data, Eff_root_Data = CalibrationTables_HLT(dfData,"Data", channelData, year, leptons, VERB)

    del dfData



    import pickle  
    print "Writing the efficiency tables to EffTable/EfficiencyTables_Calib_HLT_{}_{}_{}-q2{}_{}.pkl".format(channelData, year, "Data", version, low_val)
    pickle.dump(Eff_tables_Data, open('./EffTable/EfficiencyTables_Calib_HLT_{}_{}_{}-q2{}_{}.pkl'.format(channelData, year, "Data", version, low_val), 'wb'))
    
    #######
    print "Saving the Calibration histograms in EffTable/EffHisto_Calib_HLT_{}_{}_{}-q2{}_{}.root".format(channelData, year, "Data", version, low_val)
    file_root_Data = TFile("EffTable/EffHisto_Calib_HLT_{}_{}_{}-q2{}_{}.root".format(channelData, year, "Data", version, low_val),"RECREATE")
    map(lambda x:x.Write(), Eff_root_Data)
    #
    file_root_Data.Close()


    
    '''
    import pickle  
    print "Writing the efficiency tables to EffTable/EfficiencyTables_Calib_HLT_{}_{}_{}.pkl".format(channel, year, inputType)
    os.system('mkdir -p EffTable')
    pickle.dump(Eff_tables, open('EffTable/EfficiencyTables_Calib_HLT_{}_{}_{}.pkl'.format(channel, year, inputType), 'wb'))
    
    
    
    os.system('mkdir -p Plots')
    os.system('mkdir -p Plots/TriggerCalibration')
    
    print "Saving the Calibration histograms in EffTable/EffHisto_Calib_HLT_{}_{}_{}.root".format(channel, year, inputType)
    file_root = TFile("EffTable/EffHisto_Calib_HLT_{}_{}_{}.root".format(channel, year, inputType),"RECREATE")
    
    eff_HLT_L0E.Write()
    eff_HLT_L0H.Write()
    eff_HLT_L0TIS.Write()
    
    file_root.Close()
    '''
