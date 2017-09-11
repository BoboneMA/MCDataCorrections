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
from Correct_MC_with_Data_pro import Correct_MC_with_Data, Reweighting_L0, Reweighting_HLT

from UtilsTriggerCalib import CalibSelection_M,CalibSelection_E,GetJob,  listdirs, listfiles



def Singlefile_HLTTriggerCorrector(ifile, year, Tag_name):

    '''
    Reweights the dataframe according to the HLT Table correction.

    This is done in three steps: 1.Open the .h5 file and retrieve the dataframe 
                                 2. Reweight the full dataframe with weights that depend on the trigger category of interest;
                                 3. Save the reweighted h5
    
                                    

    '''
    directory = ifile.rsplit("/",1)[0]+"/"
    filename = ifile.rsplit("/",1)[1]
    filename_spltd = filename.split("-")
            

    print "Opening the file: %s ...."%ifile
    
    store = pd.HDFStore(ifile)
    
    if(store.keys()):
        print "The DataFrames available are: ", store.keys()
        
        df = store["dfL0"]
        
        if( "ee" in filename_spltd[0]):
            #df =CalibSelection_E(df, TM, version, low_val)

            dfRWL0E = Reweighting_HLT(df,"ee", 'L0E', year, Tag_name)
            dfRWL0H = Reweighting_HLT(dfRWL0E, "ee", 'L0H',year, Tag_name)
            dfRW = Reweighting_HLT(dfRWL0H, "ee",'L0TIS',year,Tag_name)


            os.system('mkdir -p {}'.format(directory+"L0HLTRW/"))
            print "\n\n Saving the reweighted dataframe in {} \n\n".format(directory+"L0HLTRW/"+filename.split(".")[0]+"_L0HLTrw.h5")
            dfRW.to_hdf(directory+"L0HLTRW/"+filename.split(".")[0]+"_L0Hrw.h5", format='table', key="dfL0HLT")

            store.close()
                        
                            #

        elif ("mm" in filename_spltd[0]):
            #CalibSelection_M(df, TM)
            
            dfRWL0M = Reweighting_HLT(df,"mm", 'L0M',year, Tag_name)
            dfRWL0H = Reweighting_HLT(dfRWL0M, "mm", 'L0H',year,Tag_name)
            dfRW = Reweighting_HLT(dfRWL0H, "mm", 'L0TIS',year,Tag_name)

            os.system('mkdir -p {}'.format(directory+"L0HLTRW/"))            
            print "\n\n Saving the reweighted dataframe in {} \n\n".format(directory+"L0HLTRW/"+filename.split(".")[0]+"_L0HLTrw.h5")
            dfRW.to_hdf(directory+"L0HLTRW/"+filename.split(".")[0]+"_L0HLTrw.h5", format='table', key="dfL0HLT")

            store.close()
            
        else:
            raise RuntimeError('The {} option is not available, choose between "ee,mm"'.format(leptons))
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


    Singlefile_HLTTriggerCorrector(ifile, year, Tag_name)
