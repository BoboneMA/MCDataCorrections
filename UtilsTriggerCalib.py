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

def Open_files_for_TriggerCalibration(directory, jobDict, inputType, channel, year, TM, version, low_val, test):

    ######################################



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

    joblist = []
    for iyear in yearlist:
        for imag in ['Down','Up']:
            joblist.extend(GetJob(jobDict, channel, iyear, imag, inputType))
    print "The jobs used are {} corresponding to\n--------------\n channel: {}\n--------------\n year: {}\n--------------\n inputType: {}\n--------------\n".format(joblist,channel, year, inputType)

    filecount = 0
    List = []
    for ijob in joblist:

        for ifile in listfiles(directory):
            #print(ifile)
            #print(ifile.split("-"))
            if ((".h5" in ifile) & (".h5." not in ifile) & (ifile.split("-")[-2] == ijob)):

                print "Opening the file: {} ....".format(directory+ifile)

                try:
                    store = pd.HDFStore(directory+ifile)
                
                    if(store.keys()):
                        print "The Frames available are: ", store.keys()
                    
                        try:
                            #dataframe = store['dfL0HLT']
                            dataframe = store['dfTightKst0_noTrig']
                            
                            if('ee' in channel):
                                dataframe = CalibSelection_E(dataframe, TM, version, low_val)

                            elif ('mm' in channel):
                                dataframe = CalibSelection_M(dataframe, TM)
                            else:
                                print "Calibration channel not recognized.\n Exiting ..."
                                exit()

                            print "Events in this dataframe:  ",dataframe.index.size
                            List.append(dataframe)
                            store.close()
                            filecount +=1
                            if ((filecount > max_files) & (max_files > 0)):
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


    df = pd.concat(List)

    return df

def PreselB02Kst0Jpsi2eeDF(df,version, low_val) :
    """                                                                                                                                                                                                                                
    Preselect the B02Kst0Jpsi2ee candidates by asking:                                                                                                                                                                                 
    - q2 (calculated after constraining the B0 candidate to originate from the PV and the invariant mass of the B0 candidate to correspond to the nominal B0 mass) between 7. and 11. GeV^2;                                           
    - version ('All', 'q2', or 'q2_PVandBDTF'), according to which cuts have to be applied.                                                                                                                                            
    """
    low = float(low_val)
    high = 11.
    if (version == 'All') :
        cut = ((df.q2_PVandBDTF>low) & (df.q2_PVandBDTF<high) & (df.q2>low) & (df.q2<high))
    elif (version == 'q2') :
        cut = ((df.q2>low) & (df.q2<high))
    elif (version == 'q2_PVandBDTF') :
        cut = ((df.q2_PVandBDTF>low) & (df.q2_PVandBDTF<high))

    df = df[cut]

    return df

def CalibSelection_M(df, TM):

    df = df[(df.B_PVandJpsiDTF_B_M> 5219.) & (df.B_PVandJpsiDTF_B_M <5339.) & (df.Jpsi_M > 2996.9) & (df.Jpsi_M < 3196.9 ) ]

    if(TM):
        
        df = df[(df.B_BKGCAT == 0) | (df.B_BKGCAT == 10) | (df.B_BKGCAT == 50) | (df.B_BKGCAT == 60) ]

    return df

def CalibSelection_E(df, TM, version, low_val):


    df = df[(df.B_PVandJpsiDTF_B_M> 5219.) & (df.B_PVandJpsiDTF_B_M <5339.)]
    df = PreselB02Kst0Jpsi2eeDF(df,version, low_val)


    if(TM):

        df = df[(df.B_BKGCAT == 0) | (df.B_BKGCAT == 10) | (df.B_BKGCAT == 50) | (df.B_BKGCAT == 60) ]

    return df

#############################################################
    
def HLTTriggerSelection_E(df, year):


    if ((year == '11') or (year == '12')) :
        cutHLT1 = (df.B_Hlt1TrackAllL0Decision_TOS==1)
    elif ((year == '15') or (year == '16')) :
        cutHLT1 = ((df.B_Hlt1TrackMVADecision_TOS==1) | (df.B_Hlt1TwoTrackMVADecision_TOS==1) | (df.B_Hlt1ElectronTrackDecision_TOS==1) | (df.B_Hlt1TrackMVALooseDecision_TOS==1) | (df.B_Hlt1TwoTrackMVALooseDecision_TOS==1))
        # The HLT2 lines are different for 2011-2012 and 2015-2016.                                                                                                                                                                                                         
        
    if ((year == '11') or (year == '12')) :
        cutHLT2 = ((df.B_Hlt2Topo2BodyBBDTDecision_TOS==1) | (df.B_Hlt2Topo3BodyBBDTDecision_TOS==1) | (df.B_Hlt2Topo4BodyBBDTDecision_TOS==1) | (df.B_Hlt2TopoE2BodyBBDTDecision_TOS==1) | (df.B_Hlt2TopoE3BodyBBDTDecision_TOS==1) | (df.B_Hlt2TopoE4BodyBBDTDecision_TOS==1))
    elif ((year == '15') or (year == '16')) :
        cutHLT2 = ((df.B_Hlt2Topo2BodyDecision_TOS==1) | (df.B_Hlt2Topo3BodyDecision_TOS==1) | (df.B_Hlt2Topo4BodyDecision_TOS==1) | (df.B_Hlt2TopoE2BodyDecision_TOS==1) | (df.B_Hlt2TopoE3BodyDecision_TOS==1) | (df.B_Hlt2TopoE4BodyDecision_TOS==1) | (df.B_Hlt2TopoEE2BodyDecision_TOS==1) | (df.B_Hlt2TopoEE3BodyDecision_TOS==1) | (df.B_Hlt2TopoEE4BodyDecision_TOS==1))
    
    return df[cutHLT1 & cutHLT2]

def HLTTriggerSelection_M(df, year):
    
    cutHLT1 = ((df.B_Hlt1TrackAllL0Decision_TOS==1) | (df.B_Hlt1TrackMuonDecision_TOS==1))
    cutHLT2 = ((df.B_Hlt2Topo2BodyBBDTDecision_TOS==1) | (df.B_Hlt2Topo3BodyBBDTDecision_TOS==1) | (df.B_Hlt2Topo4BodyBBDTDecision_TOS==1) | (df.B_Hlt2TopoMu2BodyBBDTDecision_TOS==1) | (df.B_Hlt2TopoMu3BodyBBDTDecision_TOS==1) | (df.B_Hlt2TopoMu4BodyBBDTDecision_TOS==1) | (df.B_Hlt2DiMuonDetachedDecision_TOS==1))


    return df[cutHLT1 & cutHLT2]


#############################################################
def listdirs(folder):
    '''
    List only directories
    '''

    return [d for d in os.listdir(folder) if os.path.isdir(os.path.join(folder, d))]


def listfiles(folder):
    '''
    List only files in the folder
    '''

    return [f for f in os.listdir(folder) if os.path.isfile(os.path.join(folder, f))]



def GetJob(jobDict, channel, year, polarity, inputType):

    job = []
    count = 0
    for ijob in jobDict.keys():
        j = jobDict[ijob]

        if ((j['year'] == year) & (j['polarity'] == polarity) & (j['inputType'] == inputType) &  (j['channel'] == channel)):
            job.append(ijob)
            count += 1

    if (count == 0):
        print 'ACHTUNG! No job in the dictionary correponds to the one you describe!'
        
    elif(count > 1 ):
        print 'ACHTUNG! More than one job found!'
    else:
        print "The corresponding job is: ", job

    return job

