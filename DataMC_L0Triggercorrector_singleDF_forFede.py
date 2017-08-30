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
from Vocabulary import  jobsDict, type_list, channel_list, yearsRI, yearsRII, mag_list
from Tools0 import GetJob, GetList_for_ReducedChain, listdirs, listfiles

from TagAndProbe_kFold_pro import TagAndProbe_L0E, TagAndProbe_L0H, TagAndProbe_L0M, TagAndProbe_L0TIS, TagAndProbe_HLT, TagAndProbe_L0E_onlynumerator, TagAndProbe_L0H_onlynumerator, TagAndProbe_L0M_onlynumerator, TagAndProbe_L0TIS_onlynumerator, TagAndProbe_HLT_onlynumerator
from reweighting import reweighting
from Plotdf import PlotDF
import pdb
from Correct_MC_with_Data_pro import Correct_MC_with_Data_E_maxET, Correct_MC_with_Data_M_maxPT, Correct_MC_with_Data_KPi, Correct_MC_with_Data_TIS, Correct_MC_with_Data_E_mixedmaxET, Correct_MC_with_Data_M_mixedmaxPT, Correct_MC_with_Data_HLT



def ComparisonPlots(H, title, yrange, xtitle, ytitle, legend, plot_name):

    '''
    Script that compacts the plotting procedure for two histograms
    '''
    gROOT.SetBatch(kTRUE)
    gROOT.SetStyle("Plain")
    gROOT.ForceStyle()

    #gStyle.SetOptStat(1111111)
    gStyle.SetOptStat(0)
    colorlist = [1,2,4,6,8,9]
    c = TCanvas('c','')
    #ROOT.SetOwnership(c, False)
    H[0].SetMinimum(yrange[0])
    H[0].SetMaximum(yrange[1])
    H[0].SetTitle(title)
    H[0].GetXaxis().SetTitle(xtitle)
    H[0].GetYaxis().SetTitle(ytitle)
    hc=0
    for h in H:

        h.SetMarkerStyle(20+hc)
        h.SetMarkerColor(colorlist[hc])
        h.SetMarkerSize(1.)
        h.SetLineColor(colorlist[hc])
        h.Draw("PESAME")
        hc +=1

    leg0 = TLegend(0.6,0.15,0.8,0.35)
    leg0.SetBorderSize(0)
    leg0.SetTextSize(0.045)
    

    hc =0
    for h in H:

        leg0.AddEntry(h,legend[hc],"l")
        hc+=1

    leg0.Draw()


    c.SaveAs(plot_name)
    print "Plotting  ==> ",plot_name

    del c
    del H
    del leg0


###################
def Reweighting(dftorw, ifile_spltd, year):

    if(ifile_spltd[4] == 'L0E'):
        channelMC = 'B02Kst0Jpsi2{}'.format('ee')
        channelData = 'B02Kst0{}'.format('ee')
    elif((ifile_spltd[4] == 'L0M')|(ifile_spltd[4] == 'L0H')|(ifile_spltd[4] == 'L0TIS')):
        channelMC = 'B02Kst0Jpsi2{}'.format('mm')
        channelData = 'B02Kst0{}'.format('mm')

    #Open the efficiency tables for the correction                                                                                                                             
    import pickle
    print "Obtaining the Efficiency tables from EffTable/EfficiencyTables_Calib_L0_{}_{}_MC.pkl".format(channelMC, year)
    Eff_MC_L0_tables = pickle.load(open('./EffTable/EfficiencyTables_Calib_L0_{}_{}_MC.pkl'.format(channelMC, year), 'rb'))
    
    print "Obtaining the Efficiency tables from EffTable/EfficiencyTables_Calib_L0_{}_{}_Data.pkl".format(channelData, year)
    Eff_Data_L0_tables = pickle.load(open('./EffTable/EfficiencyTables_Calib_L0_{}_{}_Data.pkl'.format(channelData, year), 'rb'))

    print "Obtaining the Efficiency tables from EffTable/EfficiencyTables_Calib_HLT_{}_RunI_MC.pkl".format(channelMC)
    Eff_MC_HLT_tables = pickle.load(open('./EffTable/EfficiencyTables_Calib_HLT_{}_RunI_MC.pkl'.format(channelMC), 'rb'))

    print "Obtaining the Efficiency tables from EffTable/EfficiencyTables_Calib_HLT_{}_RunI_Data.pkl".format(channelData)
    Eff_Data_HLT_tables = pickle.load(open('./EffTable/EfficiencyTables_Calib_HLT_{}_RunI_Data.pkl'.format(channelData), 'rb'))

    pdb.set_trace()

    dftorw['PT_min'] = dftorw[['K_PT','Pi_PT',"L1_PT", "L2_PT"]].min(axis=1)

    #weights
    if (ifile_spltd[4] == 'L0M'):
        dftorw['PT_max'] = dftorw[["L1_PT", "L2_PT"]].max(axis=1)
        dftorw[['wL0M_maxPT','wL0M_maxPT_err']] = dftorw.apply(lambda x: Correct_MC_with_Data_M_maxPT(x, Eff_Data_L0_tables['dfL0MeffData_maxPT'],Eff_MC_L0_tables['dfL0MeffMC_maxPT']), axis=1).apply(pd.Series)

        #HLT weights
        dftorw[['wHLT_L0M','wHLT_L0M_err']] = dftorw.apply(lambda x: Correct_MC_with_Data_HLT(x, Eff_Data_HLT_tables['dfData_HLT_L0Meff'], Eff_MC_HLT_tables['dfMC_HLT_L0Meff']), axis=1).apply(pd.Series)
    #
    elif (ifile_spltd[4] == 'L0E'):
        dftorw['L0Calo_ECAL_max_realET'] = dftorw[["L1_L0Calo_ECAL_realET", "L2_L0Calo_ECAL_realET"]].max(axis=1)
        dftorw['L0Calo_ECAL_region_max'] = dftorw[['L2_L0Calo_ECAL_realET','L1_L0Calo_ECAL_realET']].idxmax(axis=1)
        dftorw.ix[dftorw.L0Calo_ECAL_region_max.isin(['L1_L0Calo_ECAL_realET']), 'L0Calo_ECAL_region_max'] = dftorw.ix[dftorw.L0Calo_ECAL_region_max.isin(['L1_L0Calo_ECAL_realET']), 'L1_L0Calo_ECAL_region']
        dftorw.ix[dftorw.L0Calo_ECAL_region_max.isin(['L2_L0Calo_ECAL_realET']), 'L0Calo_ECAL_region_max'] = dftorw.ix[dftorw.L0Calo_ECAL_region_max.isin(['L2_L0Calo_ECAL_realET']), 'L2_L0Calo_ECAL_region']



        dftorw[['wL0E_maxET','wL0E_maxET_err']] = dftorw.apply(lambda x: Correct_MC_with_Data_E_maxET(x,Eff_Data_L0_tables['dfL0EffData_maxET'],Eff_MC_L0_tables['dfL0EffMC_maxET']), axis=1).apply(pd.Series)

        dftorw[['wHLT_L0E','wHLT_L0E_err']] = dftorw.apply(lambda x: Correct_MC_with_Data_HLT(x, Eff_Data_HLT_tables['dfData_HLT_L0Eeff'], Eff_MC_HLT_tables['dfMC_HLT_L0Eeff']), axis=1).apply(pd.Series)

    #
    elif (ifile_spltd[4] == 'L0H'):
        dftorw['Kstar_L0Calo_HCAL_region']= dftorw['K_L0Calo_HCAL_region']+dftorw['Pi_L0Calo_HCAL_region']        

        #notL0M -> L0H!                                                                                                                    
        dftorw[['wL0H_notL0M','wL0H_notL0M_err']] = dftorw[["Kstar_L0Calo_HCAL_region","Kstar_PT"]].apply(lambda x: Correct_MC_with_Data_E_maxET(x,Eff_Data_L0_tables['dfL0HeffData_notL0M'],Eff_MC_L0_tables['dfL0HeffMC_notL0M']),axis=1).apply(pd.Series)


        #also L0M -> L0H                                                                                                                                                          
        dftorw[['wL0H_alsoL0M','wL0H_notL0M_err']] = dftorw[["Kstar_L0Calo_HCAL_region","Kstar_PT"]].apply(lambda x: Correct_MC_with_Data_E_maxET(x,Eff_Data_L0_tables['dfL0HeffData_alsoL0M'],Eff_MC_L0_tables['dfL0HeffMC_alsoL0M']),axis=1).apply(pd.Series)


        dftorw[['wHLT_L0H','w_HLT_L0H_err']] = dftorw.apply(lambda x: Correct_MC_with_Data_HLT(x, Eff_Data_HLT_tables['dfData_HLT_L0Heff'],Eff_MC_HLT_tables['dfMC_HLT_L0Heff']),axis=1).apply(pd.Series)
    #
    elif (ifile_spltd[4] == 'L0TIS'):



        #notL0MH -> L0I!                                                                                                                                                          
        dftorw[['wL0TIS_notL0MH','wL0TIS_notL0MH_err']] = dftorw[["B_PT","nSPDHits"]].apply(lambda x: Correct_MC_with_Data_E_maxET(x, Eff_Data_L0_tables['dfL0TISeffData_notL0MH'], Eff_MC_L0_tables['dfL0TISeffMC_notL0MH']),axis=1).apply(pd.Series)


        #alsoL0MH -> L0I                                                                                                                                                          
        dftorw[['wL0TIS_alsoL0MH','wL0TIS_alsoL0MH_err']] = dftorw[["B_PT","nSPDHits"]].apply(lambda x: Correct_MC_with_Data_E_maxET(x, Eff_Data_L0_tables['dfL0TISeffData_alsoL0MH'], Eff_MC_L0_tables['dfL0TISeffMC_alsoL0MH']),axis=1).apply(pd.Series)


        dftorw[['wHLT_L0TIS','wHLT_L0TIS_err']] = dftorw.apply(lambda x: Correct_MC_with_Data_HLT(x, Eff_Data_HLT_tables['dfData_HLT_L0TISeff'], Eff_MC_HLT_tables['dfMC_HLT_L0TISeff']),axis=1).apply(pd.Series)



    return dftorw


###############

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


##################
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
    parser.add_argument("--TM", action="store_true", help="Truth Matched sample")
    parser.add_argument("--VERB", action="store_true", help="VERBOSE")

    

    args = parser.parse_args()
    
    # Parameters and configuration                                                                                                                             
    year= args.y
    leptons= args.l
    test = args.test
    TM = args.TM
    VERB = args.VERB

    directory = '/disk/data12/lhcb/flionett/LHCb/Analysis/AnalysisResults/ControlMode/Clean/6To11All/EffWithSymOption/'

    version = 'All'
    low_val = '6.'

    dir_L0HLT = './'

    for ifile in listfiles(directory):
        ifile = 'B02Kst0Jpsi2ee-RunI-MC-central-L0E-0.95-HOP0-BkgCat0-Ctl0-51501-Mass4500To5800-Eff.h5'
        ifile_spltd = ifile.split("-")

        if("MC" in ifile_spltd):


            print "Opening the file: %s ...."%ifile

            store = pd.HDFStore(directory+ifile)
            
            if(store.keys()):
                print "The DataFrames available are: ", store.keys()
                
                df = store[ifile_spltd[0]]

                if( "ee" in ifile_spltd[0]):
                     CalibSelection_E(df, TM, version, low_val)

                elif ("mm" in ifile_spltd[0]):
                     CalibSelection_M(df, TM)

                else:
                    print "Not clear what to do..."
                    exit()

                df11 = df[df.Year == 11]
                df12 = df[df.Year == 12]
                dfRW11 = Reweighting(df11, ifile_spltd,'11')
                dfRW12 = Reweighting(df12, ifile_spltd,'12')

                dfRW = pd.concat([dfRW11,dfRW12],axis = 1)
                dfRW.to_hdf(dir_L0HLT+ifile.split(".")[0]+"_L0HLTrw.h5", format='table', key="DF")

                store.close()
                exit()
            else:
                print '========================================================================================'
                print"WARNING! The file {} is empty, please verify that this is not a mistake.".format(ifile)
                print '========================================================================================'

        else:
            print "Jumping..."



