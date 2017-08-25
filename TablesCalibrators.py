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
from Vocabulary import  jobsDict, type_list, channel_list, yearsRI, yearsRII, mag_list
from Tools0 import GetJob, GetList_for_ReducedChain, listdirs, listfiles

from TagAndProbe_kFold import TagAndProbe_L0E, TagAndProbe_L0H, TagAndProbe_L0M, TagAndProbe_L0TIS, TagAndProbe_HLT, TagAndProbe_L0E_onlynumerator, TagAndProbe_L0H_onlynumerator, TagAndProbe_L0M_onlynumerator, TagAndProbe_L0TIS_onlynumerator, TagAndProbe_HLT_onlynumerator
from reweighting import reweighting
from Plotdf import PlotDF
import pdb
from Correct_MC_with_Data import Correct_MC_with_Data_E_maxET, Correct_MC_with_Data_M_maxPT, Correct_MC_with_Data_KPi, Correct_MC_with_Data_TIS, Correct_MC_with_Data_E_mixedmaxET, Correct_MC_with_Data_M_mixedmaxPT, Correct_MC_with_Data_HLT

from UtilsTriggerCalib import HLTTriggerSelection_E,HLTTriggerSelection_M



def CalibrationTables_L0(df, inputType, channel,  year, leptons, Tag_name, VERB):

    '''
    Author: Michele Atzeni
    Email: michele.atzeni@cern.ch
    Date: 22 Aug 2017


    Description:
    Given a dataframe and its main characteristics it produces the TISTOS efficiency tables and histograms for L0, thanks to the modules TagAndProbe* .


    Achtung:
    It is important to notice that Adaflag is True for Data and False for MC. This is beacuse the Data sample is usually smaller than the simulated one and the bnning should be adapted to the smallest possible sample. 
    The binning obtained for the data is then saved in Binning*.dat and used after for the corresponding simulated sample. 
    For this reason we should always calibrate first the data and then the MC.
    '''
    #This options allows you to select the TIS sampleon the basis of a hardcoded list you find in TagAndProb_kfold_pro.py

    Eff_tables = {}

    if(inputType == "Data"):
        Adaflag = True
    else:
        Adaflag = False

    print "The dataframe used for L0 efficiency calcultaion has shape: ",df.shape

    if(leptons == 'ee'):



        dfL0 = HLTTriggerSelection_E(df, year)
        print "Dataframe after the HLT selection has shape: ",dfL0.shape

        selTag = "B_L0Global_TIS"
        #TIS TOS efficiency for: L0E
        modelL0E = "maxET"
        effIn_e_maxET, effMid_e_maxET, effOut_e_maxET, dfL0Eff_maxET = TagAndProbe_L0E(dfL0,                                   #dataframe used 
                                                                                       selTag,                                 #selection used for the Tag sample
                                                                                       modelL0E,                               #model used for the efficiency calculation, i.e. the variable used to bin the efficiency
                                                                                       '{}-{}-{}'.format(inputType, Tag_name, year+modelL0E),#Tag used to distinguish the efficiency histograms
                                                                                       Adaflag,                                #Adaptive binning option
                                                                                       None,                                   #Name of the weights column to use 
                                                                                       VERB)                                   #VERBOSE
        '''
        #Histos with weights
        #effInData_e_maxETw, effMidData_e_maxETw, effOutData_e_maxETw, dfL0EffData_maxETw = TagAndProbe_L0E(dfDataL0, selTag, modelL0E, 'Data-SW-{}'.format(year+modelL0E), True, "Sig_sw", VERB)
        #effInMC_e_maxETw, effMidMC_e_maxETw, effOutMC_e_maxETw, dfL0EffMC_maxETw = TagAndProbe_L0E(df, selTag, modelL0E,'MC-BDT-{}'.format(year+modelL0E), False, 'weights_gb2', VERB)
        '''
        ######

        modelL0E = "mixedmaxET"
        effIn_e_mixedmaxET, effMid_e_mixedmaxET, effOut_e_mixedmaxET, dfL0Eff_mixedmaxET = TagAndProbe_L0E(dfL0, selTag, modelL0E, '{}-{}-{}'.format(inputType, Tag_name,year+modelL0E), Adaflag, None, VERB) 
        '''
        #histos with weights
        #effInData_e_mixedmaxETw, effMidData_e_mixedmaxETw, effOutData_e_mixedmaxETw, dfL0EffData_mixedmaxETw = TagAndProbe_L0E(dfDataL0, selTag, modelL0E, 'Data-SW-{}'.format(year+modelL0E), True, "Sig_sw", VERB) 
        #effInMC_e_mixedmaxETw, effMidMC_e_mixedmaxETw, effOutMC_e_mixedmaxETw, dfL0EffMC_mixedmaxETw = TagAndProbe_L0E(df, selTag, modelL0E, 'MC-BDT-{}'.format(year+modelL0E), False, "weights_gb2", VERB) 
        '''

        #TIS TOS efficiency for: L0H
        modelL0H = "notL0E"
        effIn_KPi_notL0E, effMid_KPi_notL0E, effOut_KPi_notL0E, dfL0Heff_notL0E = TagAndProbe_L0H(dfL0, selTag, modelL0H, '{}-{}-{}'.format(inputType, Tag_name,year+modelL0H), Adaflag, None, VERB)
        '''
        #effInMC_KPi_notL0Ew, effMidMC_KPi_notL0Ew, effOutMC_KPi_notL0Ew, dfL0HeffMC_notL0Ew = TagAndProbe_L0H(df, selTag, modelL0H, 'MC-BDT-{}'.format(year+modelL0H), False, 'weights_gb2', VERB)
        '''

        modelL0H = "alsoL0E" 
        effIn_KPi_alsoL0E, effMid_KPi_alsoL0E, effOut_KPi_alsoL0E, dfL0Heff_alsoL0E = TagAndProbe_L0H(dfL0, selTag, modelL0H, '{}-{}-{}'.format(inputType, Tag_name, year+modelL0H), Adaflag, None, VERB) 
        '''
        #effInMC_KPi_alsoL0Ew, effMidMC_KPi_alsoL0Ew, effOutMC_KPi_alsoL0Ew, dfL0HeffMC_alsoL0Ew = TagAndProbe_L0H(df, selTag, modelL0H, 'MC-BDT-{}'.format(year+modelL0H), False, "weights_gb2", VERB) 
        '''

        #TIS TOS efficiency for: L0TIS
        modelL0TIS = "notL0EH"
        selTag = "B_L0Global_TOS"
        eff_tis0_notL0EH, eff_tis1_notL0EH, eff_tis2_notL0EH, eff_tis3_notL0EH, dfL0TISeff_notL0EH = TagAndProbe_L0TIS(dfL0, selTag, modelL0TIS, '{}-{}-{}'.format(inputType, Tag_name, year+modelL0TIS), Adaflag, False, VERB) 

        ######
        modelL0TIS = "alsoL0EH"
        eff_tis0_alsoL0EH, eff_tis1_alsoL0EH, eff_tis2_alsoL0EH, eff_tis3_alsoL0EH, dfL0TISeff_alsoL0EH = TagAndProbe_L0TIS(dfL0, selTag, modelL0TIS, '{}-{}-{}'.format(inputType, Tag_name,year+modelL0TIS), Adaflag, False, VERB) 
        
        Eff_tables.update({"dfL0Eff{}_maxET".format(inputType):dfL0Eff_maxET,
                           "dfL0Eff{}_mixedmaxET".format(inputType):dfL0Eff_mixedmaxET,
                           "dfL0Heff{}_notL0E".format(inputType):dfL0Heff_notL0E,
                           "dfL0Heff{}_alsoL0E".format(inputType):dfL0Heff_alsoL0E,
                           "dfL0TISeff{}_alsoL0EH".format(inputType):dfL0TISeff_alsoL0EH,
                           "dfL0TISeff{}_notL0EH".format(inputType):dfL0TISeff_notL0EH})


        Eff_root = [effIn_e_maxET,
                    effMid_e_maxET,
                    effOut_e_maxET,
                    effIn_e_mixedmaxET,
                    effMid_e_mixedmaxET,
                    effOut_e_mixedmaxET,        
                    effIn_KPi_notL0E,
                    effMid_KPi_notL0E,
                    effOut_KPi_notL0E,
                    effIn_KPi_alsoL0E,
                    effMid_KPi_alsoL0E,
                    effOut_KPi_alsoL0E,
                    eff_tis0_notL0EH,
                    eff_tis1_notL0EH,
                    eff_tis2_notL0EH,
                    eff_tis3_notL0EH,
                    eff_tis0_alsoL0EH,
                    eff_tis1_alsoL0EH,
                    eff_tis2_alsoL0EH,
                    eff_tis3_alsoL0EH]

        return Eff_tables, Eff_root

 
    elif(leptons == 'mm'):

        selTag = "B_L0Global_TIS"
        #TIS TOS efficiency for: L0M
        df.reset_index(inplace=True, drop=True)
        dfL0 = HLTTriggerSelection_M(df, year)
        print "Dataframe after the HLT selection has shape: ",dfL0.shape

        modelL0M = "maxPT"
        #This options allows you to select two different approaches for the calculation of the TISTOS efficiency for L0L and L0H. The list is hardcoded in TagAndProb_kfold_pro.py
        #
        #For L0M there are two options: +maxPT -> The efficiency plotted is the ratio between TISTOS/ TIS as a function of the maximum PT of the muon couple,
        #                                         indipendently of which lepton triggered.
        #                                      -> 
        #
        #                               +mixedmaxPT
        eff_m_maxPT, dfL0Meff_maxPT = TagAndProbe_L0M(dfL0, selTag, modelL0M, '{}-{}-{}'.format(inputType, Tag_name,year+modelL0M), Adaflag, False, VERB) 
        '''
        #Histos with weights
        #effData_m_maxPTw, dfL0MeffData_maxPTw = TagAndProbe_L0M(dfDataL0, selTag, modelL0M, 'Data-SW-{}'.format(year+modelL0M), False, "Sig_sw", VERB) 
        #histos with weights
        #effMC_m_maxPTw, dfL0MeffMC_maxPTw = TagAndProbe_L0M(df, selTag, modelL0M, 'MC-BDT-{}'.format(year+modelL0M), False, "weights_gb1", VERB) 
        '''

        modelL0M = "mixedmaxPT"
        eff_m_mixedmaxPT, dfL0Meff_mixedmaxPT = TagAndProbe_L0M(dfL0, selTag, modelL0M, '{}-{}-{}'.format(inputType, Tag_name,year+modelL0M), Adaflag, False, VERB) 
        '''
        #Histos with weights
        #effData_m_mixedmaxPTw, dfL0MeffData_mixedmaxPTw = TagAndProbe_L0M(dfDataL0, selTag, modelL0M, 'Data-SW-{}'.format(year+modelL0M), False, "Sig_sw", VERB) 
        #histos with weights
        #effMC_m_mixedmaxPTw, dfL0Meff_MC_mixedmaxPTw = TagAndProbe_L0M(df, selTag, modelL0M, 'MC-BDT-{}'.format(year+modelL0M), False, "weights_gb1", VERB) 
        '''


        #TIS TOS efficiency for: L0H
        modelL0H = "notL0M" 
        effIn_KPi_notL0M, effMid_KPi_notL0M, effOut_KPi_notL0M, dfL0Heff_notL0M = TagAndProbe_L0H(dfL0, selTag, modelL0H, '{}-{}-{}'.format(inputType, Tag_name,year+modelL0H), Adaflag, None, VERB)

        '''
        #effInMC_KPi_notL0M, effMidMC_KPi_notL0M, effOutMC_KPi_notL0M, dfL0HeffMC_notL0M = TagAndProbe_L0H(dfL0, selTag, modelL0H, 'MC-BDT-{}'.format(year+modelL0H), False, "weights_gb1", VERB) 
        '''

        modelL0H = "alsoL0M" 
        effIn_KPi_alsoL0M, effMid_KPi_alsoL0M, effOut_KPi_alsoL0M, dfL0Heff_alsoL0M = TagAndProbe_L0H(dfL0, selTag, modelL0H, '{}-{}-{}'.format(inputType, Tag_name,year+modelL0H), Adaflag, None, VERB) 
        '''
        #effInMC_KPi_alsoL0Mw, effMidMC_KPi_alsoL0Mw, effOutMC_KPi_alsoL0Mw, dfL0HeffMC_alsoL0Mw = TagAndProbe_L0H(dfL0, selTag, modelL0H, 'MC-BDT-{}'.format(year+modelL0H), False, "weights_gb1", VERB) 
        '''


        #TIS TOS efficiency for: L0TIS
        selTag = "B_L0Global_TOS"
        modelL0TIS = "notL0MH"
        eff_tis0_notL0MH, eff_tis1_notL0MH, eff_tis2_notL0MH, eff_tis3_notL0MH, dfL0TISeff_notL0MH = TagAndProbe_L0TIS(dfL0, selTag, modelL0TIS, '{}-{}-{}'.format(inputType, Tag_name,year+modelL0TIS), Adaflag, False, VERB) 

        modelL0TIS = "alsoL0MH"
        eff_tis0_alsoL0MH, eff_tis1_alsoL0MH, eff_tis2_alsoL0MH, eff_tis3_alsoL0MH, dfL0TISeff_alsoL0MH = TagAndProbe_L0TIS(dfL0, selTag, modelL0TIS, '{}-{}-{}'.format(inputType, Tag_name,year+modelL0TIS), Adaflag, False, VERB) 
        #

        Eff_tables.update({"dfL0TISeff{}_alsoL0MH".format(inputType):dfL0TISeff_alsoL0MH,
                           "dfL0TISeff{}_notL0MH".format(inputType):dfL0TISeff_notL0MH,
                           "dfL0Heff{}_alsoL0M".format(inputType): dfL0Heff_alsoL0M,
                           "dfL0Heff{}_notL0M".format(inputType):dfL0Heff_notL0M,
                           "dfL0Meff{}_mixedmaxPT".format(inputType):dfL0Meff_mixedmaxPT,
                           "dfL0Meff{}_maxPT".format(inputType):dfL0Meff_maxPT})


        Eff_root = [eff_m_maxPT,
                    eff_m_mixedmaxPT,
                    effIn_KPi_notL0M,
                    effMid_KPi_notL0M,
                    effOut_KPi_notL0M,
                    effIn_KPi_alsoL0M,
                    effMid_KPi_alsoL0M,
                    effOut_KPi_alsoL0M,
                    eff_tis0_notL0MH,
                    eff_tis1_notL0MH,
                    eff_tis2_notL0MH,
                    eff_tis3_notL0MH,
                    eff_tis0_alsoL0MH,
                    eff_tis1_alsoL0MH,
                    eff_tis2_alsoL0MH,
                    eff_tis3_alsoL0MH]
        #################

        return Eff_tables, Eff_root

    else:
        print 'Warning!Wrong particles!'
        return 0,0






def CalibrationTables_HLT(df, inputType, channel,  year, leptons, Tag_name, VERB):

    '''
    Author: Michele Atzeni
    Email: michele.atzeni@cern.ch
    Date: 22 Aug 2017


    Description:
    Given a dataframe and its main characteristics it produces the TISTOS efficiency tables and histograms for HLT, thanks to the modules TagAndProbe* .


    Achtung:
    It is important to notice that Adaflag is True for Data and False for MC. This is beacuse the Data sample is usually smaller than the simulated one and the bnning should be adapted to the smallest possible sample. 
    The binning obtained for the data is then saved in Binning*.dat and used after for the corresponding simulated sample. 
    For this reason we should always calibrate first the data and then the MC.
    '''
    


    #This options allows you to select the TIS sampleon the basis of a hardcoded list you find in TagAndProb_kfold_pro.py
    Eff_tables = {}
    print "The dataframe used for HLT efficiency calcultaion has shape: ",df.shape

 
    if (inputType == 'Data'):
        Adaflag = True
    else:
        Adaflag = False


    if(leptons == 'mm'):

        #HLT efficiency
        df['PT_min'] = df[['K_PT','Pi_PT',"L1_PT", "L2_PT"]].min(axis=1)

        modelHLT = "M"

        dfL0M =df[(df.L1_L0MuonDecision_TOS == 1) | (df.L2_L0MuonDecision_TOS == 1)] 
        dfL0Hin = df[(df.Kstar_L0HadronDecision_TOS == 1) ] 
        dfL0TISin =df[(df.B_L0Global_TIS == 1)] 
        print "The three inclusive categories have the following shapes:\n L0M: {}\n L0H: {}\n L0I: {}\n".format(dfL0M.shape, dfL0Hin.shape, dfL0TISin.shape )
        
        eff_HLT_L0M, df_HLT_L0Meff = TagAndProbe_HLT(dfL0M, year, modelHLT, '{}-HLT-{}-{}'.format(inputType, Tag_name, year+modelHLT+"_L0M"), Adaflag, None, VERB)
        #If we want L0 reweighted
        '''
        #eff_HLT_L0M, df_HLT_L0Meff = TagAndProbe_HLT(dfL0M, year, modelHLT, '-L0Mrw-{}'.format(year+modelHLT+"_L0M"), False, "wL0M_maxPT", VERB) 
        #eff_HLT_L0Mw, df_HLT_L0Meffw = TagAndProbe_HLT(dfL0M, year, modelHLT, '-BDT-{}'.format(year+modelHLT+"L0M"), False, 'weights_gb2', VERB) 
        '''

        eff_HLT_L0H, df_HLT_L0Heff = TagAndProbe_HLT(dfL0Hin, year, modelHLT, '{}-HLT-{}-{}'.format(inputType, Tag_name,year+modelHLT+"_L0H"), Adaflag, None, VERB) 
        '''
        #eff_HLT_L0H, df_HLT_L0Heff = TagAndProbe_HLT(dfL0H, year, modelHLT, '-L0Hrw-{}'.format(year+modelHLT+"_L0H"), False, 'wL0H_alsoL0M', VERB) 
        #eff_HLT_L0Hw, df_HLT_L0Heffw = TagAndProbe_HLT(dfL0H, year, modelHLT, '-BDT-{}'.format(year+modelHLT+"L0H"), False, 'weights_gb2', VERB) 
        '''

        eff_HLT_L0TIS, df_HLT_L0TISeff = TagAndProbe_HLT(dfL0TISin, year, modelHLT, '{}-HLT-{}-{}'.format(inputType, Tag_name,year+modelHLT+"_L0TIS"), Adaflag, None, VERB) 
        '''
        #eff_HLT_L0TIS, df_HLT_L0TISeff = TagAndProbe_HLT(dfL0TIS, year, modelHLT, '-L0TISrw-{}'.format(year+modelHLT+"_L0TIS"), False, "wL0TIS_alsoL0MH", VERB) 
        #eff_HLT_L0TISw, df_HLT_L0TISeffw = TagAndProbe_HLT(dfDataL0TIS, year, modelHLT, '-BDT-{}'.format(year+modelHLT+"L0TIS"), False, 'weights_gb2', VERB) 
        '''
        
        Eff_tables.update({"df{}_HLT_L0Meff".format(inputType):df_HLT_L0Meff,
                           "df{}_HLT_L0Heff".format(inputType):df_HLT_L0Heff,
                           "df{}_HLT_L0TISeff".format(inputType):df_HLT_L0TISeff})

        Eff_root = [eff_HLT_L0M,
                    eff_HLT_L0H,
                    eff_HLT_L0TIS]

        return Eff_tables, Eff_root

    elif(leptons == 'ee'):

        # Reading the MC sample


        #HLT efficiency
        df['PT_min'] = df[['K_PT','Pi_PT',"L1_PT", "L2_PT"]].min(axis=1)

        modelHLT = "E"

        dfL0E =df[(df.L1_L0ElectronDecision_TOS == 1) | (df.L2_L0ElectronDecision_TOS == 1)] 
        dfL0Hin = df[df.Kstar_L0HadronDecision_TOS == 1] 
        dfL0TISin =df[df.B_L0Global_TIS == 1] 
        print "The three inclusive categories have the following shapes:\n L0M: {}\n L0H: {}\n L0I: {}\n".format(dfL0E.shape, dfL0Hin.shape, dfL0TISin.shape )

        eff_HLT_L0E, df_HLT_L0Eeff = TagAndProbe_HLT(dfL0E, year, modelHLT, '{}-HLT-{}-{}'.format(inputType, Tag_name, year+modelHLT+"_L0E"), Adaflag, None, VERB)
        #If we want L0 reweighted
        '''
        #eff_HLT_L0E, df_HLT_L0Eeff = TagAndProbe_HLT(dfL0E, year, modelHLT, '-L0Erw-{}'.format(year+modelHLT+"_L0E"), False, "wL0E_maxPT", VERB) 
        #eff_HLT_L0Ew, df_HLT_L0Eeffw = TagAndProbe_HLT(dfL0E, year, modelHLT, '-BDT-{}'.format(year+modelHLT+"L0E"), False, 'weights_gb2', VERB) 
        '''

        eff_HLT_L0H, df_HLT_L0Heff = TagAndProbe_HLT(dfL0Hin, year, modelHLT, '{}-HLT-{}-{}'.format(inputType, Tag_name, year+modelHLT+"_L0H"), Adaflag, None, VERB) 
        '''
        #eff_HLT_L0H, df_HLT_L0Heff = TagAndProbe_HLT(dfL0H, year, modelHLT, '-L0Hrw-{}'.format(year+modelHLT+"_L0H"), False, 'wL0H_alsoL0E', VERB) 
        #eff_HLT_L0Hw, df_HLT_L0Heffw = TagAndProbe_HLT(dfL0H, year, modelHLT, '-BDT-{}'.format(year+modelHLT+"L0H"), False, 'weights_gb2', VERB) 
        '''

        eff_HLT_L0TIS, df_HLT_L0TISeff = TagAndProbe_HLT(dfL0TISin, year, modelHLT, '{}-HLT-{}-{}'.format(inputType, Tag_name, year+modelHLT+"_L0TIS"), Adaflag, None, VERB) 
        '''
        #eff_HLT_L0TIS, df_HLT_L0TISeff = TagAndProbe_HLT(dfL0TIS, year, modelHLT, '-L0TISrw-{}'.format(year+modelHLT+"_L0TIS"), False, "wL0TIS_alsoL0EH", VERB) 
        #eff_HLT_L0TISw, df_HLT_L0TISeffw = TagAndProbe_HLT(dfDataL0TIS, year, modelHLT, '-BDT-{}'.format(year+modelHLT+"L0TIS"), False, 'weights_gb2', VERB) 
        '''
        
        Eff_tables.update({"df{}_HLT_L0Eeff".format(inputType):df_HLT_L0Eeff,
                           "df{}_HLT_L0Heff".format(inputType):df_HLT_L0Heff,
                           "df{}_HLT_L0TISeff".format(inputType):df_HLT_L0TISeff})
        Eff_root = [eff_HLT_L0E,
                    eff_HLT_L0H,
                    eff_HLT_L0TIS]
        return Eff_tables, Eff_root
    else:
        print 'Warning!Wrong particles!'
        return 0,0


