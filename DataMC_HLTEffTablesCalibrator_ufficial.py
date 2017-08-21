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


def CalibrationTables_HLT(df, inputType, channel,  year, leptons):

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


 
    if (inputType == 'Data'):
        Adaflag = True
    else:
        Adaflag = False





    if(leptons == 'mm'):

        #df = df[(df.B_PVandJpsiDTF_B_M>5219.) & (df.B_PVandJpsiDTF_B_M <5339.) & (data.Jpsi_M > 2996.9) & (data.Jpsi_M < 3196.9 ) ]
        print "df.shape TightKst0: ",df.shape
        Eff_tables = {}


        #HLT efficiency
        '''
        import pickle  
        df = pickle.load(open('df_{}_{}.pkl'.format(channelMC, year), 'rb'))
        print "The dataframe reweighted has the following shape: ", df.shape
        '''
        df['PT_min'] = df[['K_PT','Pi_PT',"L1_PT", "L2_PT"]].min(axis=1)

        modelHLT = "M"

        dfL0M =df[(df.L1_L0MuonDecision_TOS == 1) | (df.L2_L0MuonDecision_TOS == 1)] 
        dfL0Hin = df[(df.Kstar_L0HadronDecision_TOS == 1) ] 
        dfL0TISin =df[(df.B_L0Global_TIS == 1) ] 

        eff_HLT_L0M, df_HLT_L0Meff = TagAndProbe_HLT(dfL0M, year, modelHLT, '{}-HLT-{}'.format(inputType,year+modelHLT+"_L0M"), Adaflag, None, VERB)
        #If we want L0 reweighted
        '''
        #eff_HLT_L0M, df_HLT_L0Meff = TagAndProbe_HLT(dfL0M, year, modelHLT, '-L0Mrw-{}'.format(year+modelHLT+"_L0M"), False, "wL0M_maxPT", VERB) 
        #eff_HLT_L0Mw, df_HLT_L0Meffw = TagAndProbe_HLT(dfL0M, year, modelHLT, '-BDT-{}'.format(year+modelHLT+"L0M"), False, 'weights_gb2', VERB) 
        '''

        eff_HLT_L0H, df_HLT_L0Heff = TagAndProbe_HLT(dfL0Hin, year, modelHLT, '{}-HLT-{}'.format(inputType,year+modelHLT+"_L0H"), Adaflag, None, VERB) 
        '''
        #eff_HLT_L0H, df_HLT_L0Heff = TagAndProbe_HLT(dfL0H, year, modelHLT, '-L0Hrw-{}'.format(year+modelHLT+"_L0H"), False, 'wL0H_alsoL0M', VERB) 
        #eff_HLT_L0Hw, df_HLT_L0Heffw = TagAndProbe_HLT(dfL0H, year, modelHLT, '-BDT-{}'.format(year+modelHLT+"L0H"), False, 'weights_gb2', VERB) 
        '''

        eff_HLT_L0TIS, df_HLT_L0TISeff = TagAndProbe_HLT(dfL0TISin, year, modelHLT, '{}-HLT-{}'.format(inputType,year+modelHLT+"_L0TIS"), Adaflag, None, VERB) 
        '''
        #eff_HLT_L0TIS, df_HLT_L0TISeff = TagAndProbe_HLT(dfL0TIS, year, modelHLT, '-L0TISrw-{}'.format(year+modelHLT+"_L0TIS"), False, "wL0TIS_alsoL0MH", VERB) 
        #eff_HLT_L0TISw, df_HLT_L0TISeffw = TagAndProbe_HLT(dfDataL0TIS, year, modelHLT, '-BDT-{}'.format(year+modelHLT+"L0TIS"), False, 'weights_gb2', VERB) 
        '''
        
        Eff_tables.update({"df{}_HLT_L0Meff".format(inputType):df_HLT_L0Meff,
                           "df{}_HLT_L0Heff".format(inputType):df_HLT_L0Heff,
                           "df{}_HLT_L0TISeff".format(inputType):df_HLT_L0TISeff})

        import pickle  
        print "Writing the efficiency tables to EffTable/EfficiencyTables_Calib_HLT_{}_{}_{}.pkl".format(channel, year, inputType)
        os.system('mkdir -p EffTable')
        pickle.dump(Eff_tables, open('EffTable/EfficiencyTables_Calib_HLT_{}_{}_{}.pkl'.format(channel, year, inputType), 'wb'))



        os.system('mkdir -p Plots')
        os.system('mkdir -p Plots/TriggerCalibration')

        print "Saving the Calibration histograms in EffTable/EffHisto_Calib_HLT_{}_{}_{}.root".format(channel, year, inputType)
        file_root = TFile("EffTable/EffHisto_Calib_HLT_{}_{}_{}.root".format(channel, year, inputType),"RECREATE")
        
        eff_HLT_L0M.Write()
        eff_HLT_L0H.Write()
        eff_HLT_L0TIS.Write()


        file_root.Close()

    elif(leptons == 'ee'):
        # Reading the MC sample

        #df = df[(df.B_PVandJpsiDTF_B_M>5219.) & (df.B_PVandJpsiDTF_B_M <5339.) & (df.Jpsi_M*df.Jpsi_M  > 6000000.) & (df.Jpsi_M*df.Jpsi_M < 11000000.) ]
        print "df.shape TightKst0: ",df.shape

        #HLT efficiency
        '''
        import pickle  
        df = pickle.load(open('df_{}_{}.pkl'.format(channelMC, year), 'rb'))
        print "The dataframe reweighted has the following shape: ", df.shape
        '''
        df['PT_min'] = df[['K_PT','Pi_PT',"L1_PT", "L2_PT"]].min(axis=1)

        modelHLT = "E"

        dfL0E =df[(df.L1_L0ElectronDecision_TOS == 1) | (df.L2_L0ElectronDecision_TOS == 1)] 
        dfL0Hin = df[df.Kstar_L0HadronDecision_TOS == 1] 
        dfL0TISin =df[df.B_L0Global_TIS == 1] 
        

        eff_HLT_L0E, df_HLT_L0Eeff = TagAndProbe_HLT(dfL0E, year, modelHLT, '{}-HLT-{}'.format(inputType,year+modelHLT+"_L0E"), Adaflag, None, VERB)
        #If we want L0 reweighted
        '''
        #eff_HLT_L0E, df_HLT_L0Eeff = TagAndProbe_HLT(dfL0E, year, modelHLT, '-L0Erw-{}'.format(year+modelHLT+"_L0E"), False, "wL0E_maxPT", VERB) 
        #eff_HLT_L0Ew, df_HLT_L0Eeffw = TagAndProbe_HLT(dfL0E, year, modelHLT, '-BDT-{}'.format(year+modelHLT+"L0E"), False, 'weights_gb2', VERB) 
        '''

        eff_HLT_L0H, df_HLT_L0Heff = TagAndProbe_HLT(dfL0Hin, year, modelHLT, '{}-HLT-{}'.format(inputType,year+modelHLT+"_L0H"), Adaflag, None, VERB) 
        '''
        #eff_HLT_L0H, df_HLT_L0Heff = TagAndProbe_HLT(dfL0H, year, modelHLT, '-L0Hrw-{}'.format(year+modelHLT+"_L0H"), False, 'wL0H_alsoL0E', VERB) 
        #eff_HLT_L0Hw, df_HLT_L0Heffw = TagAndProbe_HLT(dfL0H, year, modelHLT, '-BDT-{}'.format(year+modelHLT+"L0H"), False, 'weights_gb2', VERB) 
        '''

        eff_HLT_L0TIS, df_HLT_L0TISeff = TagAndProbe_HLT(dfL0TISin, year, modelHLT, '{}-HLT-{}'.format(inputType,year+modelHLT+"_L0TIS"), Adaflag, None, VERB) 
        '''
        #eff_HLT_L0TIS, df_HLT_L0TISeff = TagAndProbe_HLT(dfL0TIS, year, modelHLT, '-L0TISrw-{}'.format(year+modelHLT+"_L0TIS"), False, "wL0TIS_alsoL0EH", VERB) 
        #eff_HLT_L0TISw, df_HLT_L0TISeffw = TagAndProbe_HLT(dfDataL0TIS, year, modelHLT, '-BDT-{}'.format(year+modelHLT+"L0TIS"), False, 'weights_gb2', VERB) 
        '''
        
        Eff_tables.update({"df{}_HLT_L0Eeff".format(inputType):df_HLT_L0Eeff,
                           "df{}_HLT_L0Heff".format(inputType):df_HLT_L0Heff,
                           "df{}_HLT_L0TISeff".format(inputType):df_HLT_L0TISeff})
        

        import pickle  
        print "Writing the efficiency tables to EffTable/EfficiencyTables_Calib_HLT_{}_{}_{}.pkl".format(channel, year, inputType)
        os.system('mkdir -p EffTable')
        pickle.dump(Eff_tables, open('EffTable/EfficiencyTables_Calib_HLT_{}_{}_{}.pkl'.format(channel, year, inputType), 'wb'))



        os.system('mkdir -p Plots')
        os.system('mkdir -p Plots/TriggerCalibration')

        print "Saving the Calibration histograms in EffTable/EffHisto_Calib_HLT_{}_{}_{}.root".format(channel, year, inputType)
        file_root = TFile("EffTable/EffHisto_Calib_HLT_{}_{}_{}.root".format(channel, year, inputType),"RECREATE")
        
        eff_HLT_L0E.Write()
        eff_HLT_L0H.Write()
        eff_HLT_L0TIS.Write()

        file_root.Close()

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
    
    parser.add_argument("-y", "--year" , dest="y"  , required=False, help="Choose the year for the RI-RII", choices=yearsRI+yearsRII+['RunI','RunII'], default = 'RunI')
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


    if (year == 'RunI'):
        yearlist = ['11','12']
    elif (year == 'RunII'):
        yearlist = ['15','16']
    else:
        yearlist = [year]

    joblistData = []
    for iyear in yearlist:
        for imag in ['Down','Up']:
            joblistData.append(GetJob( channelData, iyear, imag, "Data"))

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
    ######################################
    CalibrationTables_HLT(dfData,"Data", channelData, year, leptons)
    del dfData

    #MCrootfile = 'B02Kst0Jpsi2{}-MC-20{}-ForRW_AddBranchesForBDTReweighting.root'.format(leptons,year)
    #Datarootfile = 'B02Kst0Jpsi2{}-{}_sWeighted_Preselected_AddBranchesForBDTReweighting.root'.format(leptons,year)
    #MCrootfile = 'B02Kst0mm-MC-2016-Beam6500GeV-Nu1.6-MagDown-Sim09b-Reco16-11144001.root'
    #Datarootfile = 'B02Kst0mm-MC-2016-Beam6500GeV-Nu1.6-MagDown-Sim09b-Reco16-11144001.root'

    joblistMC = []
    for iyear in yearlist:
        for imag in ['Down','Up']:
            joblistMC.append(GetJob( channelMC, iyear, imag, "MC"))
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

                                branches = ['L1_PT','L2_PT','K_PT','Pi_PT',"B_PT","B_Hlt1TrackAllL0Decision_TOS","B_Hlt1TrackMuonDecision_TOS","B_Hlt2Topo2BodyBBDTDecision_TOS","B_Hlt2Topo3BodyBBDTDecision_TOS","B_Hlt2Topo4BodyBBDTDecision_TOS","B_Hlt2TopoMu2BodyBBDTDecision_TOS","B_Hlt2TopoMu3BodyBBDTDecision_TOS","B_Hlt2TopoMu4BodyBBDTDecision_TOS","B_Hlt2DiMuonDetachedDecision_TOS","B_Hlt1Phys_TIS","B_Hlt2Phys_TIS","L1_L0MuonDecision_TOS","L2_L0MuonDecision_TOS","Kstar_L0HadronDecision_TOS","B_L0Global_TIS","Kstar_PT"]
                                # Reading the MC sample
                                data = data[branches]

                            if(leptons == 'ee'):
                                data = data[(data.B_PVandJpsiDTF_B_M> 5219.) & (data.B_PVandJpsiDTF_B_M <5339.) & (data.Jpsi_M*data.Jpsi_M > 6000000.) & (data.Jpsi_M*data.Jpsi_M < 11000000.)]
 
                                if ((year == '11') or (year == '12') or (year == 'RunI')) :
                                    branches = ['L1_PT','L2_PT','K_PT','Pi_PT',"B_PT","B_Hlt1TrackAllL0Decision_TOS","B_Hlt2Topo2BodyBBDTDecision_TOS","B_Hlt2Topo3BodyBBDTDecision_TOS","B_Hlt2Topo4BodyBBDTDecision_TOS","B_Hlt2TopoE2BodyBBDTDecision_TOS","B_Hlt2TopoE3BodyBBDTDecision_TOS","B_Hlt2TopoE4BodyBBDTDecision_TOS","B_Hlt1Phys_TIS","B_Hlt2Phys_TIS","L1_L0ElectronDecision_TOS","L2_L0ElectronDecision_TOS","Kstar_L0HadronDecision_TOS","B_L0Global_TIS","Kstar_PT"]
                                elif ((year == '15') or (year == '16')or (year == 'RunII')) :
                                    branches = ['L1_PT','L2_PT','K_PT','Pi_PT',"B_PT","B_Hlt1TrackMVADecision_TOS","B_Hlt1TwoTrackMVADecision_TOS","B_Hlt1ElectronTrackDecision_TOS","B_Hlt1TrackMVALooseDecision_TOS","B_Hlt1TwoTrackMVALooseDecision_TOS","B_Hlt2Topo2BodyDecision_TOS","B_Hlt2Topo3BodyDecision_TOS","B_Hlt2Topo4BodyDecision_TOS","B_Hlt2TopoE2BodyDecision_TOS","B_Hlt2TopoE3BodyDecision_TOS","B_Hlt2TopoE4BodyDecision_TOS","B_Hlt2TopoEE2BodyDecision_TOS","B_Hlt2TopoEE3BodyDecision_TOS","B_Hlt2TopoEE4BodyDecision_TOS","B_Hlt1Phys_TIS","B_Hlt2Phys_TIS","L1_L0ElectronDecision_TOS","L2_L0ElectronDecision_TOS","Kstar_L0HadronDecision_TOS","B_L0Global_TIS","Kstar_PT"]

                                data = data[branches]

                    
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


    CalibrationTables_HLT(dfMC,"MC", channelMC, year, leptons)
    
