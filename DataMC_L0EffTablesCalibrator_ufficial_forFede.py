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

#from TagAndProbe_kFold_pro import TagAndProbe_L0E, TagAndProbe_L0H, TagAndProbe_L0M, TagAndProbe_L0TIS, TagAndProbe_HLT, TagAndProbe_L0E_onlynumerator, TagAndProbe_L0H_onlynumerator, TagAndProbe_L0M_onlynumerator, TagAndProbe_L0TIS_onlynumerator, TagAndProbe_HLT_onlynumerator
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
    
    parser.add_argument("-y", "--year" , dest="y"  , required=False, help="Choose the year for the RI-RII", choices=yearsRI+yearsRII+['RunI','RunII'], default = '11')
    parser.add_argument("-l", "--leptons" , dest="l"  , required=False, help="Choose the leptons in the final state", choices= ['ee','mm'], default = 'ee')

    parser.add_argument("--test", action="store_true", help="Do a test plot")
    parser.add_argument("--plot", action="store_true", help="Do a test plot")
    parser.add_argument("--VERB", action="store_true", help="VERBOSE")

    

    args = parser.parse_args()
    
    # Parameters and configuration                                                                                                                             
    year= args.y
    leptons= args.l
    test = args.test
    plot = args.plot
    VERB = args.VERB

    channelData = 'B02Kst0{}'.format(leptons)
    channelMC = 'B02Kst0Jpsi2{}'.format(leptons)

    #Part 0: Obtain two dataframes, one for Data and one for MC

    if (test):
        max_files = 20
    else:
        max_files = 1000000000000
    ######################################
    dir_baseMC = '/home/hep/matzeni/gangadir/Analysis/DFs/MC/{}/DFs/WithHOP/ParallelSplitOutput/'.format(channelMC)
    dir_baseData = '/home/hep/matzeni/gangadir/Analysis/DFs/Data/{}/DFs/WithHOP/ParallelSplitOutput/'.format(channelData)

    joblistMC = []
    for imag in ['Down','Up']:
        joblistMC.append(GetJob( channelMC, year, imag, "MC"))
    print "The jobs used are {} corresponding to\n--------------\n channel: {}\n--------------\n year: {}\n--------------\n inputType: {}\n--------------\n".format(joblistMC,channelMC, year, "MC")

    filecount = 0
    MCList = []
    for ijob in joblistMC:

        for ifile in os.listdir(dir_baseMC+ijob+"/"):
            if ((".h5" in ifile) & (".h5." not in ifile)):

                print "Opening the file: %s ...."%ifile

                try:
                    store = pd.HDFStore(dir_baseMC+ijob+"/"+ifile)
                
                    if(store.keys()):
                        print "The MCFrames available are: ", store.keys()
                    
                        try:
                            data = store['dfTightKst0_noTrig']


                            if(leptons == 'mm'):
                                data = data[(data.B_PVandJpsiDTF_B_M> 5219.) & (data.B_PVandJpsiDTF_B_M <5339.) & (data.Jpsi_M > 2996.9) & (data.Jpsi_M < 3196.9 ) ]
                            if(leptons == 'ee'):
                                data = data[(data.B_PVandJpsiDTF_B_M> 5219.) & (data.B_PVandJpsiDTF_B_M <5339.) & (data.Jpsi_M*data.Jpsi_M > 6000000.) & (data.Jpsi_M*data.Jpsi_M < 11000000.)]

                            print "Events in this dataset:  ",data.index.size
                            MCList.append(data)
                            store.close()
                            filecount +=1
                            if (filecount > max_files):
                                break
                        except:
                            print "ACHTUNG!No dfTightKst0_noTrig"
                            store.close()
                    else:
                        print '========================================================================================'
                        print"WARNING! The file {} is empty, please verify that this is not a mistake.".format(ifile)
                        print '========================================================================================'
                
                except:
                    print "ACHTUNG!PROBLEM!Skipping..."
                    print "***************************"


    dfMC = pd.concat(MCList)


    print "================ Log File =========================\n"
    print '''
The MCDataCorrection for L0 trigger is usually done \n
on the control samples, for this reason the .root \n
files used are:\n\n
====================================================\n     
MC => {}\n
====================================================\n     
Data => {}\n    
====================================================\n     '''.format(dir_baseMC,dir_baseData)

    CalibrationTables_L0(dfMC,"MC", channelMC, year, leptons)
    
    del dfMC


    joblistData = []
    for imag in ['Down','Up']:
        joblistData.append(GetJob( channelData, year, imag, "Data"))

    dir_baseMC = '/home/hep/matzeni/gangadir/Analysis/DFs/MC/{}/DFs/WithHOP/ParallelSplitOutput/'.format(channelMC)
    dir_baseData = '/home/hep/matzeni/gangadir/Analysis/DFs/Data/{}/DFs/WithHOP/ParallelSplitOutput/'.format(channelData)
    #dir_base = '/home/hep/matzeni/gangadir/workspace/matzeni/LocalXML/212/0/output/'

    print "The jobs used are {} corresponding to\n--------------\n channel: {}\n--------------\n year: {}\n--------------\n inputType: {}\n--------------\n".format(joblistData,channelData, year, "Data")


    DataList = []
    for ijob in joblistData:
        filecount = 0
        for ifile in os.listdir(dir_baseData+ijob+"/"):
            if ((".h5" in ifile) & (".h5." not in ifile)):


                print "Opening the file: %s ...."%ifile

                try:
                    store = pd.HDFStore(dir_baseData+ijob+"/"+ifile)
                
                    if(store.keys()):
                        print "The DataFrames available are: ", store.keys()
                        
                        try:
                            data = store['dfTightKst0_noTrig']

                            if(leptons == 'mm'):
                                data = data[(data.B_PVandJpsiDTF_B_M> 5219.) & (data.B_PVandJpsiDTF_B_M <5339.) & (data.Jpsi_M > 2996.9) & (data.Jpsi_M < 3196.9 ) ]
                            if(leptons == 'ee'):
                                data = data[(data.B_PVandJpsiDTF_B_M> 5219.) & (data.B_PVandJpsiDTF_B_M <5339.) & (data.Jpsi_M*data.Jpsi_M > 6000000.) & (data.Jpsi_M*data.Jpsi_M < 11000000.)]

                            print "Events in this dataset:  ",data.index.size
                            
                            DataList.append(data)
                            store.close()
                            filecount +=1
                            if (filecount > max_files):
                                break       
                        except:
                            print "ACHTUNG!No dfTightKst0_noTrig"
                            store.close()
                    else:
                        print '========================================================================================'
                        print"WARNING! The file {} is empty, please verify that this is not a mistake.".format(ifile)
                        print '========================================================================================'
                except:
                    print "ACHTUNG!PROBLEM!Skipping..."
                    print "***************************"
                



    dfData = pd.concat(DataList)
                           

    CalibrationTables_L0(dfData,"Data", channelData, year, leptons)
    #MCrootfile = 'B02Kst0Jpsi2{}-MC-20{}-ForRW_AddBranchesForBDTReweighting.root'.format(leptons,year)
    #Datarootfile = 'B02Kst0Jpsi2{}-{}_sWeighted_Preselected_AddBranchesForBDTReweighting.root'.format(leptons,year)
    #MCrootfile = 'B02Kst0mm-MC-2016-Beam6500GeV-Nu1.6-MagDown-Sim09b-Reco16-11144001.root'
    #Datarootfile = 'B02Kst0mm-MC-2016-Beam6500GeV-Nu1.6-MagDown-Sim09b-Reco16-11144001.root'
