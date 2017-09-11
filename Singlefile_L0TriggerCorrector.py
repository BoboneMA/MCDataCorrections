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
    Reweights the single ".h5" file according to the Reweighting_L0 module 

    Author: Michele Atzeni
    Date: June 1st, 2017

    '''


    directory = ifile.rsplit("/",1)[0]+"/"
    filename = ifile.rsplit("/",1)[1]
    filename_spltd = filename.split("-")
            
            
    print "Opening the file: %s ...."%ifile
    
    store = pd.HDFStore(ifile)
    
    if(store.keys()):
        print "The DataFrames available are: ", store.keys()

        #TO DO: get rid of this hardcoded dataframe
        df = store["dfTightKst0_noTrig"]
        sizedf = df.shape[0]
        if( "ee" in filename_spltd[0]):



            dfRWL0E = Reweighting_L0(df, "ee", 'L0E', year, Tag_name)
            dfRWL0H = Reweighting_L0(dfRWL0E, "ee", 'L0H',year, Tag_name)
            dfRW = Reweighting_L0(dfRWL0H, "ee",'L0TIS',year,Tag_name)

            #check
            if((dfRW.shape[0]>sizedf )|(dfRW.shape[0]<sizedf)):
                raise RuntimeError("Consistency check failed, error found")
                exit()


            os.system('mkdir -p {}'.format(directory+"L0RW/"))
            print "\n\n Saving the reweighted dataframe in {} \n\n".format(directory+"L0RW/"+filename.split(".")[0]+"_L0rw.h5")
            dfRW.to_hdf(directory+"L0RW/"+filename.split(".")[0]+"_L0rw.h5", format='table', key="dfL0")

            store.close()
                        
        #
        elif ("mm" in filename_spltd[0]):

            
            dfRWL0M = Reweighting_L0(df, "mm" ,'L0M',year, Tag_name)
            dfRWL0H = Reweighting_L0(dfRWL0M, "mm",'L0H',year,Tag_name)
            dfRW = Reweighting_L0(dfRWL0H, "mm",'L0TIS',year,Tag_name)

            if((dfRW.shape[0]>sizedf )|(dfRW.shape[0]<sizedf)):
                raise RuntimeError( "Consistency check failed, error found")
                exit()

            os.system('mkdir -p {}'.format(directory+"L0RW/"))            
            print "\n\n Saving the reweighted dataframe in {} \n\n".format(directory+"L0RW/"+filename.split(".")[0]+"_L0rw.h5")
            dfRW.to_hdf(directory+"L0RW/"+filename.split(".")[0]+"_L0rw.h5", format='table', key="dfL0")

            store.close()
            
        else:
            print "Not clear what to do, nor ee or mm in filename"
            exit()
                            
    else:
        print '========================================================================================'
        print"WARNING! The file {} is empty, please verify that this is not a mistake.".format(ifile)
        print '========================================================================================'
                

if __name__ == "__main__" :


    '''
    Author: Michele Atzeni
    Date: June 1st, 2017
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
