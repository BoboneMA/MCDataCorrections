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
from Tools0 import GetJob, GetList_for_ReducedChain, listdirs

#from TagAndProbe_kFold_pro import TagAndProbe_L0E, TagAndProbe_L0H, TagAndProbe_L0M, TagAndProbe_L0TIS, TagAndProbe_HLT, TagAndProbe_L0E_onlynumerator, TagAndProbe_L0H_onlynumerator, TagAndProbe_L0M_onlynumerator, TagAndProbe_L0TIS_onlynumerator, TagAndProbe_HLT_onlynumerator
from TagAndProbe_kFold import TagAndProbe_L0E, TagAndProbe_L0H, TagAndProbe_L0M, TagAndProbe_L0TIS, TagAndProbe_HLT, TagAndProbe_L0E_onlynumerator, TagAndProbe_L0H_onlynumerator, TagAndProbe_L0M_onlynumerator, TagAndProbe_L0TIS_onlynumerator, TagAndProbe_HLT_onlynumerator
from reweighting import reweighting
from Plotdf import PlotDF
import pdb
from Correct_MC_with_Data import Correct_MC_with_Data_E_maxET, Correct_MC_with_Data_M_maxPT, Correct_MC_with_Data_KPi, Correct_MC_with_Data_TIS, Correct_MC_with_Data_E_mixedmaxET, Correct_MC_with_Data_M_mixedmaxPT, Correct_MC_with_Data_HLT



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

#TRIGGER CUTS

def L0TriggerSelection_E(df):
    
    cutL0 = ((df.L1_L0ElectronDecision_TOS==1) | (df.L2_L0ElectronDecision_TOS==1) | (df.K_L0HadronDecision_TOS==1) | (df.Pi_L0HadronDecision_TOS==1) | (df.B_L0Global_TIS==1))     
    return df[cutL0]

def L0TriggerSelection_M(df):
    
    cutL0 = (df.B_L0MuonDecision_TOS==1)
    return df[cutL0]
    
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


def CalibrationTables_L0(df, inputType, channel,  year, leptons):

    #This options allows you to select the TIS sampleon the basis of a hardcoded list you find in TagAndProb_kfold_pro.py
    
    selTag = "B_L0Global_TIS"
    Eff_tables = {}
    #### Part 1: Apply some selections to dfData and dfMC necessary for our studies
    #
    #### For muons: + B_PVandJpsiDTF_B_M in [5219, 5339.] MeV/c^2
    #               + Jpsi_M in [2996.9, 3196.9] MeV/c^2
    #               + Hlt1 && Hlt2
    #
    #For electrons:  + B_PVandJpsiDTF_B_M in [5219., 5339.] MeV/c^2                                                                                       
    #               + Jpsi_M^2 in [6000000., 11000000.] MeV^2/c^4 
    #               + Hlt1 && Hlt2                                                                                                                            
    #                          



    if(inputType == "Data"):
        Adaflag = True
    else:
        Adaflag = False



    if(leptons == 'ee'):

        # Reading the MC sample

        df = df[(df.B_PVandJpsiDTF_B_M>5219.) & (df.B_PVandJpsiDTF_B_M <5339.) & (df.Jpsi_M*df.Jpsi_M  > 6000000.) & (df.Jpsi_M*df.Jpsi_M < 11000000.) ]
        print "df.shape TightKst0: ",df.shape
        dfL0 = HLTTriggerSelection_E(df, year)

        ##Part 2: Efficiency plots and tables for Data and MC

        #TIS TOS efficiency for: L0E

        modelL0E = "maxET"
        effIn_e_maxET, effMid_e_maxET, effOut_e_maxET, dfL0Eff_maxET = TagAndProbe_L0E(dfL0, selTag, modelL0E, '{}-{}'.format(inputType,year+modelL0E), Adaflag, None, VERB)
        '''
        #Histos with weights
        #effInData_e_maxETw, effMidData_e_maxETw, effOutData_e_maxETw, dfL0EffData_maxETw = TagAndProbe_L0E(dfDataL0, selTag, modelL0E, 'Data-SW-{}'.format(year+modelL0E), True, "Sig_sw", VERB)
        #effInMC_e_maxETw, effMidMC_e_maxETw, effOutMC_e_maxETw, dfL0EffMC_maxETw = TagAndProbe_L0E(df, selTag, modelL0E,'MC-BDT-{}'.format(year+modelL0E), False, 'weights_gb2', VERB)
        '''
        ######

        modelL0E = "mixedmaxET"
        effIn_e_mixedmaxET, effMid_e_mixedmaxET, effOut_e_mixedmaxET, dfL0Eff_mixedmaxET = TagAndProbe_L0E(dfL0, selTag, modelL0E, '{}-{}'.format(inputType,year+modelL0E), Adaflag, None, VERB) 
        '''
        #histos with weights
        #effInData_e_mixedmaxETw, effMidData_e_mixedmaxETw, effOutData_e_mixedmaxETw, dfL0EffData_mixedmaxETw = TagAndProbe_L0E(dfDataL0, selTag, modelL0E, 'Data-SW-{}'.format(year+modelL0E), True, "Sig_sw", VERB) 
        #effInMC_e_mixedmaxETw, effMidMC_e_mixedmaxETw, effOutMC_e_mixedmaxETw, dfL0EffMC_mixedmaxETw = TagAndProbe_L0E(df, selTag, modelL0E, 'MC-BDT-{}'.format(year+modelL0E), False, "weights_gb2", VERB) 
        '''

        #TIS TOS efficiency for: L0H

        modelL0H = "notL0E"
        effIn_KPi_notL0E, effMid_KPi_notL0E, effOut_KPi_notL0E, dfL0Heff_notL0E = TagAndProbe_L0H(dfL0, selTag, modelL0H, '{}-{}'.format(inputType,year+modelL0H), Adaflag, None, VERB)
        '''
        #effInMC_KPi_notL0Ew, effMidMC_KPi_notL0Ew, effOutMC_KPi_notL0Ew, dfL0HeffMC_notL0Ew = TagAndProbe_L0H(df, selTag, modelL0H, 'MC-BDT-{}'.format(year+modelL0H), False, 'weights_gb2', VERB)
        '''

        modelL0H = "alsoL0L" 
        effIn_KPi_alsoL0E, effMid_KPi_alsoL0E, effOut_KPi_alsoL0E, dfL0Heff_alsoL0E = TagAndProbe_L0H(dfL0, selTag, modelL0H, '{}-{}'.format(inputType, year+modelL0H), Adaflag, None, VERB) 
        '''
        #effInMC_KPi_alsoL0Ew, effMidMC_KPi_alsoL0Ew, effOutMC_KPi_alsoL0Ew, dfL0HeffMC_alsoL0Ew = TagAndProbe_L0H(df, selTag, modelL0H, 'MC-BDT-{}'.format(year+modelL0H), False, "weights_gb2", VERB) 
        '''

        #TIS TOS efficiency for: L0TIS
        
        modelL0TIS = "notL0EH"
        selTagTIS = "B_L0Global_TOS"
        eff_tis0_notL0EH, eff_tis1_notL0EH, eff_tis2_notL0EH, eff_tis3_notL0EH, dfL0TISeff_notL0EH = TagAndProbe_L0TIS(dfL0, selTagTIS, modelL0TIS, '{}-{}'.format(inputType, year+modelL0TIS), Adaflag, False, VERB) 

        #

        modelL0TIS = "alsoL0LH"
        eff_tis0_alsoL0EH, eff_tis1_alsoL0EH, eff_tis2_alsoL0EH, eff_tis3_alsoL0EH, dfL0TISeff_alsoL0EH = TagAndProbe_L0TIS(dfL0, selTagTIS, modelL0TIS, '{}-{}'.format(inputType,year+modelL0TIS), Adaflag, False, VERB) 
        
        Eff_tables.update({"dfL0Eff{}_maxET".format(inputType):dfL0Eff_maxET,
                           "dfL0Eff{}_mixedmaxET".format(inputType):dfL0Eff_mixedmaxET,
                           "dfL0Heff{}_notL0E".format(inputType):dfL0Heff_notL0E,
                           "dfL0Heff{}_alsoL0E".format(inputType):dfL0Heff_alsoL0E,
                           "dfL0TISeff{}_alsoL0EH".format(inputType):dfL0TISeff_alsoL0EH,
                           "dfL0TISeff{}_notL0EH".format(inputType):dfL0TISeff_notL0EH})

        import pickle  
        print "Writing the efficiency tables to EffTable/EfficiencyTables_Calib_L0_{}_{}_{}.pkl".format(channel, year, inputType)
        os.system('mkdir -p EffTable')
        pickle.dump(Eff_tables, open('./EffTable/EfficiencyTables_Calib_L0_{}_{}_{}.pkl'.format(channel, year, inputType), 'wb'))

        #####
        os.system('mkdir -p Plots')
        os.system('mkdir -p Plots/TriggerCalibration')

        # L0E Plots Data/MC 
        print "Saving the Calibration histograms in EffTable/EffHisto_Calib_L0_{}_{}_{}.root".format(channel, year, inputType)
        file_root = TFile("EffTable/EffHisto_Calib_L0_{}_{}_{}.root".format(channel, year, inputType),"RECREATE")

        effIn_e_maxET.Write()
        effMid_e_maxET.Write()
        effOut_e_maxET.Write()
        effIn_e_mixedmaxET.Write()
        effMid_e_mixedmaxET.Write()
        effOut_e_mixedmaxET.Write()        
        effIn_KPi_notL0E.Write()
        effMid_KPi_notL0E.Write()
        effOut_KPi_notL0E.Write()
        effIn_KPi_alsoL0E.Write()
        effMid_KPi_alsoL0E.Write()
        effOut_KPi_alsoL0E.Write()
        eff_tis0_notL0EH.Write()
        eff_tis1_notL0EH.Write()
        eff_tis2_notL0EH.Write()
        eff_tis3_notL0EH.Write()
        eff_tis0_alsoL0EH.Write()
        eff_tis1_alsoL0EH.Write()
        eff_tis2_alsoL0EH.Write()
        eff_tis3_alsoL0EH.Write()

        #
        file_root.Close()


 
    elif(leptons == 'mm'):


        #TIS TOS efficiency for: L0M

        # Reading the MC sample
        df.reset_index(inplace=True, drop=True)
        df = df[(df.B_PVandJpsiDTF_B_M>5219.) & (df.B_PVandJpsiDTF_B_M <5339.) & (df.Jpsi_M > 2996.9) & (df.Jpsi_M < 3196.9) ]
        print "df.shape TightKst0: ",df.shape
        dfL0 = HLTTriggerSelection_M(df, year)



        modelL0M = "maxPT"
        #This options allows you to select two different approaches for the calculation of the TISTOS efficiency for L0L and L0H. The list is hardcoded in TagAndProb_kfold_pro.py
        #
        #For L0M there are two options: +maxPT -> The efficiency plotted is the ratio between TISTOS/ TIS as a function of the maximum PT of the muon couple,
        #                                         indipendently of which lepton triggered.
        #                                      -> 
        #
        #                               +mixedmaxPT
        eff_m_maxPT, dfL0Meff_maxPT = TagAndProbe_L0M(dfL0, selTag, modelL0M, '{}-{}'.format(inputType,year+modelL0M), Adaflag, False, VERB) 
        '''
        #Histos with weights
        #effData_m_maxPTw, dfL0MeffData_maxPTw = TagAndProbe_L0M(dfDataL0, selTag, modelL0M, 'Data-SW-{}'.format(year+modelL0M), False, "Sig_sw", VERB) 
        #histos with weights
        #effMC_m_maxPTw, dfL0MeffMC_maxPTw = TagAndProbe_L0M(df, selTag, modelL0M, 'MC-BDT-{}'.format(year+modelL0M), False, "weights_gb1", VERB) 
        '''

        modelL0M = "mixedmaxPT"

        eff_m_mixedmaxPT, dfL0Meff_mixedmaxPT = TagAndProbe_L0M(dfL0, selTag, modelL0M, '{}-{}'.format(inputType,year+modelL0M), Adaflag, False, VERB) 
        '''
        #Histos with weights
        #effData_m_mixedmaxPTw, dfL0MeffData_mixedmaxPTw = TagAndProbe_L0M(dfDataL0, selTag, modelL0M, 'Data-SW-{}'.format(year+modelL0M), False, "Sig_sw", VERB) 
        #histos with weights
        #effMC_m_mixedmaxPTw, dfL0Meff_MC_mixedmaxPTw = TagAndProbe_L0M(df, selTag, modelL0M, 'MC-BDT-{}'.format(year+modelL0M), False, "weights_gb1", VERB) 
        '''


        #TIS TOS efficiency for: L0H
        modelL0H = "notL0M" 
        effIn_KPi_notL0M, effMid_KPi_notL0M, effOut_KPi_notL0M, dfL0Heff_notL0M = TagAndProbe_L0H(dfL0, selTag, modelL0H, '{}-{}'.format(inputType,year+modelL0H), Adaflag, None, VERB)

        '''
        #effInMC_KPi_notL0M, effMidMC_KPi_notL0M, effOutMC_KPi_notL0M, dfL0HeffMC_notL0M = TagAndProbe_L0H(dfL0, selTag, modelL0H, 'MC-BDT-{}'.format(year+modelL0H), False, "weights_gb1", VERB) 
        '''

        modelL0H = "alsoL0L" 
        effIn_KPi_alsoL0M, effMid_KPi_alsoL0M, effOut_KPi_alsoL0M, dfL0Heff_alsoL0M = TagAndProbe_L0H(dfL0, selTag, modelL0H, '{}-{}'.format(inputType,year+modelL0H), Adaflag, None, VERB) 
        '''
        #effInMC_KPi_alsoL0Mw, effMidMC_KPi_alsoL0Mw, effOutMC_KPi_alsoL0Mw, dfL0HeffMC_alsoL0Mw = TagAndProbe_L0H(dfL0, selTag, modelL0H, 'MC-BDT-{}'.format(year+modelL0H), False, "weights_gb1", VERB) 
        '''


        #TIS TOS efficiency for: L0TIS
        selTagTIS = "B_L0Global_TOS"
        modelL0TIS = "notL0MH"
        eff_tis0_notL0MH, eff_tis1_notL0MH, eff_tis2_notL0MH, eff_tis3_notL0MH, dfL0TISeff_notL0MH = TagAndProbe_L0TIS(dfL0, selTagTIS, modelL0TIS, '{}-{}'.format(inputType,year+modelL0TIS), Adaflag, False, VERB) 

        modelL0TIS = "alsoL0LH"
        eff_tis0_alsoL0MH, eff_tis1_alsoL0MH, eff_tis2_alsoL0MH, eff_tis3_alsoL0MH, dfL0TISeff_alsoL0MH = TagAndProbe_L0TIS(dfL0, selTagTIS, modelL0TIS, '{}-{}'.format(inputType,year+modelL0TIS), Adaflag, False, VERB) 
        #

        Eff_tables.update({"dfL0TISeff{}_alsoL0MH".format(inputType):dfL0TISeff_alsoL0MH,
                           "dfL0TISeff{}_notL0MH".format(inputType):dfL0TISeff_notL0MH,
                           "dfL0Heff{}_alsoL0M".format(inputType): dfL0Heff_alsoL0M,
                           "dfL0Heff{}_notL0M".format(inputType):dfL0Heff_notL0M,
                           "dfL0Meff{}_mixedmaxPT".format(inputType):dfL0Meff_mixedmaxPT,
                           "dfL0Meff{}_maxPT".format(inputType):dfL0Meff_maxPT})

        import pickle  
        print "Writing the efficiency tables to EffTable/EfficiencyTables_Calib_L0_{}_{}_{}.pkl".format(channel, year, inputType)
        os.system('mkdir -p EffTable')
        pickle.dump(Eff_tables, open('./EffTable/EfficiencyTables_Calib_L0_{}_{}_{}.pkl'.format(channel, year, inputType), 'wb'))

        #######
        os.system('mkdir -p Plots')
        os.system('mkdir -p Plots/TriggerCalibration')

        
        print "Saving the Calibration histograms in EffTable/EffHisto_Calib_L0_{}_{}_{}.root".format(channel, year, inputType)
        file_root = TFile("EffTable/EffHisto_Calib_L0_{}_{}_{}.root".format(channel, year, inputType),"RECREATE")

        eff_m_maxPT.Write()
        eff_m_mixedmaxPT.Write()
        effIn_KPi_notL0M.Write()
        effMid_KPi_notL0M.Write()
        effOut_KPi_notL0M.Write()
        effIn_KPi_alsoL0M.Write()
        effMid_KPi_alsoL0M.Write()
        effOut_KPi_alsoL0M.Write()
        eff_tis0_notL0MH.Write()
        eff_tis1_notL0MH.Write()
        eff_tis2_notL0MH.Write()
        eff_tis3_notL0MH.Write()
        eff_tis0_alsoL0MH.Write()
        eff_tis1_alsoL0MH.Write()
        eff_tis2_alsoL0MH.Write()
        eff_tis3_alsoL0MH.Write()

        #
        file_root.Close()

        #################

    else:
        print 'Warning!Wrong particles!'





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
