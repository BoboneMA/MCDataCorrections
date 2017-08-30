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



def Singlefile_L0TriggerCorrector(ifile, year, Tag_name):

    '''

    '''
    directory = ifile.rsplit("/",1)[0]+"/"
    filename = ifile.rsplit("/",1)[1]
    filename_spltd = filename.split("-")
            
            
    print "Opening the file: %s ...."%ifile
    
    store = pd.HDFStore(ifile)
    
    if(store.keys()):
        print "The DataFrames available are: ", store.keys()
        
        df = store["dfTightKst0_noTrig"]
        
        if( "ee" in filename_spltd[0]):
            #df =CalibSelection_E(df, TM, version, low_val)

            dfL0E = df[(df.L2_L0ElectronDecision_TOS == 1) | (df.L1_L0ElectronDecision_TOS == 1)]
            dfL0H = df[(df.Kstar_L0HadronDecision_TOS==1) & (df.L2_L0ElectronDecision_TOS == 0) & (df.L1_L0ElectronDecision_TOS == 0)]
            dfL0TIS = df[(df.B_L0Global_TIS==1) & (df.Kstar_L0HadronDecision_TOS==0) & (df.L2_L0ElectronDecision_TOS == 0) & (df.L1_L0ElectronDecision_TOS == 0)]

            dfRWL0L = Reweighting_L0(dfL0E, 'L0E', year, Tag_name)
            dfRWL0H = Reweighting_L0(dfL0H, 'L0H',year, Tag_name)
            dfRWL0TIS = Reweighting_L0(dfL0TIS, 'L0TIS',year,Tag_name)
            dfRW = pd.concat([dfRWL0L,dfRWL0H, dfL0TIS],axis = 0)

            os.system('mkdir -p {}'.format(directory+"L0RW/"))
            print "\n\n Saving the reweighted dataframe in {} \n\n".format(directory+"L0RW/"+filename.split(".")[0]+"_L0Hrw.h5")
            dfRW.to_hdf(directory+"L0RW/"+filename.split(".")[0]+"_L0Hrw.h5", format='table', key="dfL0HLT")

            store.close()
                        
                            #

        elif ("mm" in filename_spltd[0]):
            #CalibSelection_M(df, TM)
            
            dfL0M = df[(df.L2_L0MuonDecision_TOS == 1) | (df.L1_L0MuonDecision_TOS == 1)]
            dfL0H = df[(df.Kstar_L0HadronDecision_TOS==1) & (df.L2_L0MuonDecision_TOS == 0) & (df.L1_L0MuonDecision_TOS == 0)]
            dfL0TIS = df[(df.B_L0Global_TIS==1) & (df.Kstar_L0HadronDecision_TOS==0) & (df.L2_L0MuonDecision_TOS == 0) & (df.L1_L0MuonDecision_TOS == 0)]
            
            dfRWL0L = Reweighting_L0(dfL0M, 'L0M',year, Tag_name)
            dfRWL0H = Reweighting_L0(dfL0H, 'L0H',year,Tag_name)
            dfRWL0TIS = Reweighting_L0(dfL0TIS, 'L0TIS',year,Tag_name)
            dfRW = pd.concat([dfRWL0L,dfRWL0H, dfL0TIS],axis = 0)
            os.system('mkdir -p {}'.format(directory+"L0RW/"))            
            print "\n\n Saving the reweighted dataframe in {} \n\n".format(directory+"L0RW/"+filename.split(".")[0]+"_L0Hrw.h5")
            dfRW.to_hdf(directory+"L0RW/"+filename.split(".")[0]+"_L0Hrw.h5", format='table', key="dfL0HLT")

            store.close()
            
        else:
            print "Not clear what to do..."
            exit()
                            
    else:
        print '========================================================================================'
        print"WARNING! The file {} is empty, please verify that this is not a mistake.".format(ifile)
        print '========================================================================================'
                

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
    
    parser.add_argument("-f", "--file" , dest="ifile"  , required=False, help="Full path of the file to use")
    parser.add_argument("-y", "--year" , dest="y"  , required=False, help="Choose the year for the RI-RII", choices=['11','12','15','16','RunI','RunII'], default = '11')
    parser.add_argument("-Tag", "--Tag_name" , dest="Tag_name"  , required=False, help="Tag name for the tables to use")


    

    args = parser.parse_args()
    
    # Parameters and configuration                                                                                                                             
    year= args.y
    ifile= args.ifile
    Tag_name= args.Tag_name


    Singlefile_L0TriggerCorrector(ifile, year, Tag_name)
