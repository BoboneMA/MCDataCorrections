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
import pdb

def Correct_MC_with_Data(df, tableData,tableMC):

    '''
    Author: Michele Atzeni
    Email: michele.atzeni@cern.ch
    Date: 22 Aug 2017

    Description:
    This script returns the weigth (and its error) necessary to correct the MC event to better emulate the Data.
    The weights are obtained per event thanks to the calibration tables {tableData,tableMC}.

    Notice that this script only works if the two tables:
    1. The MC table does not have null entries
    2. The empty entries are common between the two tables.
    '''
    #Print a table randomly to check that the weights are obtained in the correct way
    rndm = np.random.uniform()
    if((rndm > 1./500.)&(rndm <= 2./500.)):
        VERB = True
    else:
        VERB =False

    #initialize weights
    W = -999
    W_err = -999


    y_column= [i for i in tableData.columns if "Binsy" in i][0]
    var_y = y_column.split('-')[1]
    yBins = np.unique(np.reshape([i for i in tableMC.loc[:,y_column].dropna().values],2*tableMC.loc[:,y_column].dropna().values.size))
    y = np.digitize(df[var_y], yBins)#Notice that in this case is just one, it would be better to do all of them together


    x_column= [i for i in tableData.columns if "Binsx" in i][0]
    var_x = x_column.split('-')[1]
    xBins = np.unique(np.reshape([i for i in tableMC.loc[:,x_column].dropna().values],2*tableMC.loc[:,x_column].dropna().values.size))
    x = np.digitize(df[var_x], xBins)                                                




    if(VERB):
        print "Correct_MC_with_Data...."#!!!This has to be modified
        print "Row:\n",df[[var_x,var_y]].T#!!!!
        print "MC efficiency table: \n",tableMC
        print "Data efficiency table: \n",tableData


    #Reduced table
    t_MC = tableMC.drop([x_column,y_column],axis=1).dropna()
    t_Data = tableData.drop([x_column,y_column],axis=1).dropna()




    if((x<1)|(x>= len(xBins))|(y<1)|(y>= len(yBins))): #Give a 0 weight to the events that are not in the intervals
        W = 0
        W_err = 0
        if(VERB):
            print "**************************"
            print "Weight:" ,W
            print "**************************"
        
    else:
        wMC = t_MC.iloc[y-1, 2*(x-1)]
        wData = t_Data.iloc[y-1, 2*(x-1)]
        wMC_err = t_MC.iloc[y-1, 2*(x-1)+1]
        wData_err = t_Data.iloc[y-1, 2*(x-1)+1]
        if (wMC > 0):
            W = wData/wMC
        else:
            if(wData > 0):
                print "FATAL ERROR!\nThe efficiency of the MC is zero while the data is not, no reweighting possible.\n Choose a better binning."
                exit()
            else:
                W = 1
        if ((wMC_err>0) & (wData_err>0)):
            W_err = np.sqrt((wMC_err/wMC)*(wMC_err/wMC) + (wData_err/wData)*(wData_err/wData) )
        else:
            if(wData_err>0):
                W_err =  (wData_err/wData)
            elif(wMC_err>0):
                W_err = (wMC_err/wMC)
            elif( (wMC_err==0) & (wData_err==0)):
                W_err = 0
            else:
                print "FATAL ERROR! It is not possible to have negative errors."
                exit()
        if(VERB):
            print "**************************"
            print "Weight: {}/{} = {} +/- {}".format(wData,wMC, W, W_err*W)
            print "**************************"
    if(np.isnan(W)):
        print "FATAL ERROR! NaN value cannot be accecpted as a weight."
        exit()
        
    return W, W_err

def Reweighting(dftorw, trig_cat, year, TM, version, low_val):

    '''
    Author: Michele Atzeni
    Email: michele.atzeni@cern.ch
    Date: 22 Aug 2017

    Description:
    The script reweights the dataframe depending on its trigger category.


    '''

    if(trig_cat == 'L0E'):
        channelMC = 'B02Kst0Jpsi2{}'.format('ee')
        channelData = 'B02Kst0{}'.format('ee')
    elif((trig_cat == 'L0M')|(trig_cat == 'L0H')|(trig_cat == 'L0TIS')):
        channelMC = 'B02Kst0Jpsi2{}'.format('mm')
        channelData = 'B02Kst0{}'.format('mm')

    if(TM):
        Tag_sample = 'q2{}_lw{}_TM'.format(version, low_val)
    else:
        Tag_sample = 'q2{}_lw{}'.format(version, low_val)



    #Open the efficiency tables for the correction                                                                                                                             
    import pickle
    print "Obtaining the Efficiency tables from EffTable/EfficiencyTables_Calib_L0_{}_{}_{}-{}.pkl".format(channelMC, year, "MC", Tag_sample)
    Eff_MC_L0_tables = pickle.load(open("./EffTable/EfficiencyTables_Calib_L0_{}_{}_{}-{}.pkl".format(channelMC, year, "MC", Tag_sample)))
    
    print "Obtaining the Efficiency tables from EffTable/EfficiencyTables_Calib_L0_{}_{}_{}-q2{}_lw{}.pkl".format(channelData, year, "Data", version, low_val)
    Eff_Data_L0_tables = pickle.load(open('./EffTable/EfficiencyTables_Calib_L0_{}_{}_{}-q2{}_lw{}.pkl'.format(channelData, year, "Data", version, low_val)))

    print "Obtaining the Efficiency tables from EffTable/EfficiencyTables_Calib_HLT_{}_{}_{}-{}.pkl".format(channelMC, "RunI", "MC", Tag_sample)
    Eff_MC_HLT_tables = pickle.load(open('./EffTable/EfficiencyTables_Calib_HLT_{}_{}_{}-{}.pkl'.format(channelMC, "RunI", "MC", Tag_sample)))

    print "Obtaining the Efficiency tables from EffTable/EfficiencyTables_Calib_HLT_{}_{}_{}-q2{}_lw{}.pkl".format(channelData, "RunI", "Data", version, low_val)
    Eff_Data_HLT_tables = pickle.load(open('./EffTable/EfficiencyTables_Calib_HLT_{}_{}_{}-q2{}_lw{}.pkl'.format(channelData, "RunI", "Data", version, low_val)))



    dftorw['PT_min'] = dftorw[['K_PT','Pi_PT',"L1_PT", "L2_PT"]].min(axis=1)

    #weights
    if (trig_cat == 'L0M'):
        dftorw['PT_max'] = dftorw[["L1_PT", "L2_PT"]].max(axis=1)
        dftorw[['wL0M_maxPT','wL0M_maxPT_err']] = dftorw.apply(lambda x: Correct_MC_with_Data(x, Eff_Data_L0_tables['dfL0MeffData_maxPT'],Eff_MC_L0_tables['dfL0MeffMC_maxPT']), axis=1).apply(pd.Series)

        #HLT weights
        dftorw[['wHLT_L0M','wHLT_L0M_err']] = dftorw.apply(lambda x: Correct_MC_with_Data(x, Eff_Data_HLT_tables['dfData_HLT_L0Meff'], Eff_MC_HLT_tables['dfMC_HLT_L0Meff']), axis=1).apply(pd.Series)
    #
    elif (trig_cat == 'L0E'):
        dftorw['L0Calo_ECAL_max_realET'] = dftorw[["L1_L0Calo_ECAL_realET", "L2_L0Calo_ECAL_realET"]].max(axis=1)
        dftorw['L0Calo_ECAL_region_max'] = dftorw[['L2_L0Calo_ECAL_realET','L1_L0Calo_ECAL_realET']].idxmax(axis=1)
        dftorw.ix[dftorw.L0Calo_ECAL_region_max.isin(['L1_L0Calo_ECAL_realET']), 'L0Calo_ECAL_region_max'] = dftorw.ix[dftorw.L0Calo_ECAL_region_max.isin(['L1_L0Calo_ECAL_realET']), 'L1_L0Calo_ECAL_region']
        dftorw.ix[dftorw.L0Calo_ECAL_region_max.isin(['L2_L0Calo_ECAL_realET']), 'L0Calo_ECAL_region_max'] = dftorw.ix[dftorw.L0Calo_ECAL_region_max.isin(['L2_L0Calo_ECAL_realET']), 'L2_L0Calo_ECAL_region']



        dftorw[['wL0E_maxET','wL0E_maxET_err']] = dftorw.apply(lambda x: Correct_MC_with_Data(x,Eff_Data_L0_tables['dfL0EffData_maxET'],Eff_MC_L0_tables['dfL0EffMC_maxET']), axis=1).apply(pd.Series)

        dftorw[['wHLT_L0E','wHLT_L0E_err']] = dftorw.apply(lambda x: Correct_MC_with_Data(x, Eff_Data_HLT_tables['dfData_HLT_L0Eeff'], Eff_MC_HLT_tables['dfMC_HLT_L0Eeff']), axis=1).apply(pd.Series)

    #
    elif (trig_cat == 'L0H'):
        dftorw['Kstar_L0Calo_HCAL_region']= dftorw['K_L0Calo_HCAL_region']+dftorw['Pi_L0Calo_HCAL_region']        

        #notL0M -> L0H!                                                                                                                    
        dftorw[['wL0H_notL0M','wL0H_notL0M_err']] = dftorw[["Kstar_L0Calo_HCAL_region","Kstar_PT"]].apply(lambda x: Correct_MC_with_Data(x,Eff_Data_L0_tables['dfL0HeffData_notL0M'],Eff_MC_L0_tables['dfL0HeffMC_notL0M']),axis=1).apply(pd.Series)


        #also L0M -> L0H                                                                                                                                                          
        dftorw[['wL0H_alsoL0M','wL0H_notL0M_err']] = dftorw[["Kstar_L0Calo_HCAL_region","Kstar_PT"]].apply(lambda x: Correct_MC_with_Data(x,Eff_Data_L0_tables['dfL0HeffData_alsoL0M'],Eff_MC_L0_tables['dfL0HeffMC_alsoL0M']),axis=1).apply(pd.Series)


        dftorw[['wHLT_L0H','w_HLT_L0H_err']] = dftorw.apply(lambda x: Correct_MC_with_Data(x, Eff_Data_HLT_tables['dfData_HLT_L0Heff'],Eff_MC_HLT_tables['dfMC_HLT_L0Heff']),axis=1).apply(pd.Series)
    #
    elif (trig_cat == 'L0TIS'):



        #notL0MH -> L0I!                                                                                                                                                          
        dftorw[['wL0TIS_notL0MH','wL0TIS_notL0MH_err']] = dftorw[["B_PT","nSPDHits"]].apply(lambda x: Correct_MC_with_Data(x, Eff_Data_L0_tables['dfL0TISeffData_notL0MH'], Eff_MC_L0_tables['dfL0TISeffMC_notL0MH']),axis=1).apply(pd.Series)


        #alsoL0MH -> L0I                                                                                                                                                          
        dftorw[['wL0TIS_alsoL0MH','wL0TIS_alsoL0MH_err']] = dftorw[["B_PT","nSPDHits"]].apply(lambda x: Correct_MC_with_Data(x, Eff_Data_L0_tables['dfL0TISeffData_alsoL0MH'], Eff_MC_L0_tables['dfL0TISeffMC_alsoL0MH']),axis=1).apply(pd.Series)


        dftorw[['wHLT_L0TIS','wHLT_L0TIS_err']] = dftorw.apply(lambda x: Correct_MC_with_Data(x, Eff_Data_HLT_tables['dfData_HLT_L0TISeff'], Eff_MC_HLT_tables['dfMC_HLT_L0TISeff']),axis=1).apply(pd.Series)



    return dftorw


