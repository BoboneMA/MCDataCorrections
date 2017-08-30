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
#from Vocabulary import  jobsDict, type_list, channel_list, yearsRI, yearsRII, mag_list
#from Tools0 import GetJob, GetList_for_ReducedChain, listdirs, listfiles

from TagAndProbe_kFold_pro import TagAndProbe_L0E, TagAndProbe_L0H, TagAndProbe_L0M, TagAndProbe_L0TIS, TagAndProbe_HLT, TagAndProbe_L0E_onlynumerator, TagAndProbe_L0H_onlynumerator, TagAndProbe_L0M_onlynumerator, TagAndProbe_L0TIS_onlynumerator, TagAndProbe_HLT_onlynumerator
from reweighting import reweighting
from Plotdf import PlotDF
import pdb
from Correct_MC_with_Data_pro import Correct_MC_with_Data, Reweighting_L0

from UtilsTriggerCalib import CalibSelection_M,CalibSelection_E,GetJob,  listdirs, listfiles
from Singlefile_L0TriggerCorrector import Singlefile_L0TriggerCorrector


def multi_helper(args):
    return Singlefile_L0TriggerCorrector(*args)



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
    
    parser.add_argument("-y", "--year" , dest="y"  , required=False, help="Choose the year for the RI-RII", choices=['11','12','15','16','RunI','RunII'], default = '11')
    parser.add_argument("-l", "--leptons" , dest="l"  , required=False, help="Choose the leptons in the final state", choices= ['ee','mm'], default = 'ee')
    parser.add_argument("-v", "--versionq2" , dest="version"  , required=False, help="", choices= ['All','q2','q2_PVandBDTF'], default = 'All')
    parser.add_argument("-low", "--lowq2" , dest="low_val"  , required=False, help="", choices= ['6','7'], default = '6')
    parser.add_argument("-u", "--user" , dest="user"  , required=False, help="", choices= ['M','F'], default = 'M')
    parser.add_argument("-c", "--cores" , dest="cores"  , required=False, help="", default = 4)

    parser.add_argument("--test", action="store_true", help="Do a test plot")
    parser.add_argument("--version", action="store_true", help="q2version")
    parser.add_argument("--TM", action="store_true", help="TruthMatchedMC")
    parser.add_argument("--VERB", action="store_true", help="VERBOSE")
    parser.add_argument("--b", action="store_true", help="BATCH")
    
    

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
    batch = args.b
    cores = args.cores

    if (test):
        max_files = 10
    else:
        max_files = -1


    if (year == 'RunI'):
        yearlist = ['11','12']
    elif (year == 'RunII'):
        yearlist = ['15','16']
    else:
        yearlist = [year]



    if(TM):
        Tag_name = 'q2{}_lw{}_TM'.format(version, low_val)
        #Tag_name = 'q2{}_lw{}_TM_PID'.format(version, low_val)                                                                                                               
        #Tag_name = 'q2{}_lw{}_TM_ECALmask'.format(version, low_val)                                                                                                             
        #Tag_name = 'q2{}_lw{}_TM_KstarPT4'.format(version, low_val)                                                                                                              
    else:
        Tag_name = 'q2{}_lw{}'.format(version, low_val)
        #Tag_name = 'q2{}_lw{}_PID'.format(version, low_val)                                                                                                                      
        #Tag_name = 'q2{}_lw{}_ECALmask'.format(version, low_val)                                                                                                                 
        #Tag_name = 'q2{}_lw{}_KstarPT4'.format(version, low_val) 






    if (user == 'M'):
        from Vocabulary import  jobsDict, type_list, channel_list, yearsRI, yearsRII, mag_list


        #directory = '/disk/data12/lhcb/flionett/LHCb/Analysis/AnalysisResults/ControlMode/Clean/6To11All/EffWithSymOption/'
        directory = '/home/hep/matzeni/gangadir/Analysis/DFs/MC/B02Kst0Jpsi2{}/DFs/WithHOP/ParallelSplitOutput/L0HLT/L0RW/'.format(leptons)
        #version = 'All'
        #low_val = '6.'

        
        channel = 'B02Kst0Jpsi2{}'.format(leptons)
        inputType = 'MC'
        joblist = []
        for iyear in yearlist:
            for imag in ['Down','Up']:
                joblist.extend(GetJob(jobsDict, channel, iyear, imag, inputType))
        print "The jobs used are {} corresponding to\n--------------\n channel: {}\n--------------\n year: {}\n--------------\n inputType: {}\n--------------\n".format(joblist,channel, year, inputType)

        filecount = 0
        List = []
        for ijob in joblist:

            for ifile in listfiles(directory):
                #print(ifile)
                #print(ifile.split("-"))
                if ((".h5" in ifile) & (".h5." not in ifile) & (ifile.split("-")[-2] == ijob)):

            
                    if(batch):
                        os.system('sbatch -p express --nodes=1  -t 02:00:00  --ntasks-per-node=2 --export=ifile={},year={},Tag_name={}   submit_HLTCorrection.slurm'.format(ifile, year,Tag_name))
                    
                        if (test):
                            print '================================================='
                            print 'As a test only the first TORQUE job was submitted'
                            print 'The file used is \n', ifile
                            print '================================================='
                            break

                    else:
                        List.append((directory+ifile, jobsDict[ijob]["year"], Tag_name))


        if(not batch):
            print(List)
            import multiprocessing
            multiprocessing.Pool(int(cores)).map(multi_helper, List)

