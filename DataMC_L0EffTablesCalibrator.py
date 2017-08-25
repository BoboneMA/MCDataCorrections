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
    Date: 22 Aug 2017

    Description:

    '''
        
    parser = argparse.ArgumentParser(description = 'Configuration of the parameters for the SplitAndMerge')
    
    parser.add_argument("-y", "--year" , dest="y"  , required=False, help="Choose the year for the RI-RII", choices=['11','12','15','16','RunI','RunII'], default = '11')
    parser.add_argument("-l", "--leptons" , dest="l"  , required=False, help="Choose the leptons in the final state", choices= ['ee','mm'], default = 'ee')
    parser.add_argument("-v", "--versionq2" , dest="version"  , required=False, help="", choices= ['All','q2','q2_PVandBDTF'], default = 'All')
    parser.add_argument("-low", "--lowq2" , dest="low_val"  , required=False, help="", choices= ['6','7'], default = '6')
    parser.add_argument("-u", "--user" , dest="user"  , required=False, help="", choices= ['M','F'], default = 'M')

    parser.add_argument("--test", action="store_true", help="Do a test plot")
    parser.add_argument("--version", action="store_true", help="q2version")
    parser.add_argument("--TM", action="store_true", help="TruthMatchedMC")
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

    #TM = True
    #version = 'All'
    #low_val = '6'
    channelData = 'B02Kst0{}'.format(leptons)
    channelMC = 'B02Kst0Jpsi2{}'.format(leptons)
    

    if(user == "M"):
        directoryMC = '/home/hep/matzeni/gangadir/Analysis/DFs/MC/{}/DFs/WithHOP/ParallelSplitOutput/L0HLT/'.format(channelMC)
        directoryData = '/home/hep/matzeni/gangadir/Analysis/DFs/Data/{}/DFs/WithHOP/ParallelSplitOutput/L0HLT/'.format(channelData)

        from Vocabulary import  jobsDict, type_list, channel_list, yearsRI, yearsRII, mag_list
        from Tools0 import GetJob,  listdirs, listfiles
    else:
        directoryMC = '/disk/data12/lhcb/flionett/LHCb/Analysis/Ntuples/MC/{}/DFs/WithHOP/ParallelSplitOutput/Meerkat/'.format( channelMC)
        directoryData = '/disk/data12/lhcb/flionett/LHCb/Analysis/Ntuples/Data/{}/DFs/WithHOP/ParallelSplitOutput/L0HLT/'.format( channelData)
        from voc import  jobsDict, type_list, channel_list, yearsRI, yearsRII, mag_list
        from tools import GetJob,  listdirs, listfiles


    if(TM):
        Tag_name = 'q2{}_lw{}_TM'.format(version, low_val)
    else:
        Tag_name = 'q2{}_lw{}'.format(version, low_val)


    dfData = Open_files_for_TriggerCalibration(directoryData, jobsDict, 'Data', channelData, year, False, version, low_val, test)
    Eff_tables_Data, Eff_root_Data = CalibrationTables_L0(dfData,"Data", channelData, year, leptons,  Tag_name, VERB)

    del dfData
    ###############
    
    dfMC = Open_files_for_TriggerCalibration(directoryMC, jobsDict,'MC', channelMC, year, TM, version, low_val,  test)
    Eff_tables_MC, Eff_root_MC = CalibrationTables_L0(dfMC,"MC", channelMC, year, leptons, Tag_name,VERB)
    
    del dfMC
    #######
    os.system('mkdir -p EffTable')
    os.system('mkdir -p EffTable/{}'.format(Tag_name))
    import pickle  
    print "Writing the efficiency tables to EffTable/{}/EfficiencyTables_Calib_L0-{}_{}_{}_{}-{}.pkl".format(Tag_name, "MC", channelMC, year,"mBoth",  Tag_name)
    with open('./EffTable/{}/EfficiencyTables_Calib_L0-{}_{}_{}_{}-{}.pkl'.format(Tag_name, "MC", channelMC, year,"mBoth",  Tag_name), 'wb') as f:
        pickle.dump(Eff_tables_MC, f)
    
    
        
    print "Saving the Calibration histograms in EffTable/{}/EffHisto_Calib_L0-{}_{}_{}_{}-{}.pkl".format(Tag_name, "MC", channelMC, year,"mBoth",  Tag_name)
    file_root_MC = TFile("EffTable/{}/EffHisto_Calib_L0{}_{}_{}_{}-{}.pkl".format(Tag_name, "MC", channelMC, year,"mBoth",  Tag_name),"RECREATE")

    map(lambda x:x.Write(), Eff_root_MC)
    file_root_MC.Close()
    





    import pickle  
    print "Writing the efficiency tables to EffTable/{}/EfficiencyTables_Calib_L0-{}_{}_{}_{}-{}.pkl".format(Tag_name, "Data", channelData, year, "mBoth", Tag_name)
    with open('./EffTable/{}/EfficiencyTables_Calib_L0-{}_{}_{}_{}-{}.pkl'.format(Tag_name, "Data", channelData, year, "mBoth", Tag_name), 'wb') as f:
        pickle.dump(Eff_tables_Data, f)
    
    #######
    print "Saving the Calibration histograms in EffTable/{}/EffHisto_Calib_L0-{}_{}_{}_{}-{}.pkl".format(Tag_name, "Data", channelData, year, "mBoth", Tag_name)
    file_root_Data = TFile("EffTable/{}/EffHisto_Calib_L0-{}_{}_{}_{}-{}.pkl".format(Tag_name, "Data", channelData, year, "mBoth", Tag_name),"RECREATE")
    map(lambda x:x.Write(), Eff_root_Data)
    #
    file_root_Data.Close()

    
