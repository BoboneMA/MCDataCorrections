
"""
Author: R. Coutinho
Date: 10/2016

Description:
 Simple method to estimate the L0 corrections using a tag & probe approach using B0 -> J/psi K*(Kpi), J/psi->mumu 
 using mu MUTOS as tag and examining the corresponding hadron pair. 

   
How to run it:
Called by e.g. python L0Tables.py -t K [K, Pi] -y 2011 [2011, 2012] -m MagDown [MagDown, MagUp, MagAll] -c Posit [Posit, Negat] -p inner outer
"""

print(__doc__)

import sys
import argparse


import ROOT
from ROOT import gROOT, gStyle, TFile, TH2F , TH1F,  TCanvas , TLegend
from itertools import repeat
import itertools as it
import csv
import re

from numpy import array
from root_numpy import  array2root, tree2array 
import numpy as np

from array import array 
import pandas as pd
import matplotlib.mlab as mlab
import matplotlib.pyplot as plt
from Adaptive_binning_pro import AdaptiveBinning1D
import pdb
#Before using this script be  sure all the variables used here are in the dataframe
'''

'''



def TAP_L0E2_mixedmaxET(df, selTag, VERB):




   # Defining the triggered sample (the one we are interested in), the Tag sample and the intersection of the two
   #selTrig = (df.K_L0HadronDecision_TOS==1) & (df.Pi_L0HadronDecision_TOS==1)
   selTrig11 = ((df.L1_L0ElectronDecision_TOS==1) & (df.L2_L0ElectronDecision_TOS==1))
   selTrig10 = ((df.L1_L0ElectronDecision_TOS==1) & (df.L2_L0ElectronDecision_TOS==0))
   selTrig01 = ((df.L1_L0ElectronDecision_TOS==0) & (df.L2_L0ElectronDecision_TOS==1))

   #selTag = (df.B_L0Global_TIS == 1)
   selProb11    = selTag &  selTrig11
   selProb10    = selTag &  selTrig10
   selProb01    = selTag &  selTrig01
   
   ## L1_TOS==1 & L2_TOS==1
   dfTag11 = df[selTag][['L1_L0Calo_ECAL_region','L1_L0Calo_ECAL_realET','L2_L0Calo_ECAL_region','L2_L0Calo_ECAL_realET','weights']]
   dfTag11['L0Calo_ECAL_region_max'] = dfTag11[['L2_L0Calo_ECAL_realET','L1_L0Calo_ECAL_realET']].idxmax(axis=1)
   dfTag11['L0Calo_ECAL_max_realET'] = dfTag11[["L1_L0Calo_ECAL_realET", "L2_L0Calo_ECAL_realET"]].max(axis=1)
   dfTag11.ix[dfTag11.L0Calo_ECAL_region_max.isin(['L1_L0Calo_ECAL_realET']), 'L0Calo_ECAL_region_max'] = dfTag11.ix[dfTag11.L0Calo_ECAL_region_max.isin(['L1_L0Calo_ECAL_realET']), 'L1_L0Calo_ECAL_region']
   dfTag11.ix[dfTag11.L0Calo_ECAL_region_max.isin(['L2_L0Calo_ECAL_realET']), 'L0Calo_ECAL_region_max'] = dfTag11.ix[dfTag11.L0Calo_ECAL_region_max.isin(['L2_L0Calo_ECAL_realET']), 'L2_L0Calo_ECAL_region']
   arrayTag11 = dfTag11[['L0Calo_ECAL_region_max','L0Calo_ECAL_max_realET','weights']].values


   dfTagAndProb11 = df[selProb11][['L1_L0Calo_ECAL_region','L1_L0Calo_ECAL_realET','L2_L0Calo_ECAL_region','L2_L0Calo_ECAL_realET','weights']]
   dfTagAndProb11['L0Calo_ECAL_max_realET'] = dfTagAndProb11[["L1_L0Calo_ECAL_realET", "L2_L0Calo_ECAL_realET"]].max(axis=1)
   dfTagAndProb11['L0Calo_ECAL_region_max'] = dfTagAndProb11[['L2_L0Calo_ECAL_realET','L1_L0Calo_ECAL_realET']].idxmax(axis=1)
   dfTagAndProb11.ix[dfTagAndProb11.L0Calo_ECAL_region_max.isin(['L1_L0Calo_ECAL_realET']), 'L0Calo_ECAL_region_max'] = dfTagAndProb11.ix[dfTagAndProb11.L0Calo_ECAL_region_max.isin(['L1_L0Calo_ECAL_realET']), 'L1_L0Calo_ECAL_region']
   dfTagAndProb11.ix[dfTagAndProb11.L0Calo_ECAL_region_max.isin(['L2_L0Calo_ECAL_realET']), 'L0Calo_ECAL_region_max'] = dfTagAndProb11.ix[dfTagAndProb11.L0Calo_ECAL_region_max.isin(['L2_L0Calo_ECAL_realET']), 'L2_L0Calo_ECAL_region']


   arrayTagAndProb11 = dfTagAndProb11[['L0Calo_ECAL_region_max','L0Calo_ECAL_max_realET','weights']].values

   ## L1_TOS==1 & L2_TOS==0
   dfTag10 = df[selTag]
   arrayTag10 = dfTag10[['L1_L0Calo_ECAL_region','L1_L0Calo_ECAL_realET','weights']].values
   dfTagAndProb10 = df[selProb10]
   arrayTagAndProb10 = df[selProb10][['L1_L0Calo_ECAL_region','L1_L0Calo_ECAL_realET','weights']].values

   dfTag01 = df[selTag]
   arrayTag01 = dfTag01[['L2_L0Calo_ECAL_region','L2_L0Calo_ECAL_realET','weights']].values
   dfTagAndProb01 = df[selProb01]
   arrayTagAndProb01 = df[selProb01][['L2_L0Calo_ECAL_region','L2_L0Calo_ECAL_realET','weights']].values
   

   print "=================== TAP_L0E2_mixedmaxET ============================\n\n"

   print "TIS & {},{},{}\\\\".format(dfTag11.shape[0],dfTag01.shape[0],dfTag10.shape[0])
   print "L0E1 and L0E2 & {}\\\\".format(dfTagAndProb11.shape[0])
   print "L0E1 and !L0E2 & {}\\\\".format(dfTagAndProb10.shape[0])
   print "!L0E1 and L0E2 & {}\\\\".format(dfTagAndProb01.shape[0])

   if VERB:

      print '========================== Appending the maximum ET column  ==============================='
      print "dfTag:\n", dfTag11[["L1_L0Calo_ECAL_realET", "L2_L0Calo_ECAL_realET",'L0Calo_ECAL_max_realET']].head()
      print "dfTagAndProb11:\n", dfTagAndProb11[["L1_L0Calo_ECAL_realET", "L2_L0Calo_ECAL_realET",'L0Calo_ECAL_max_realET']].head()
      print "dfTagAndProb11:\n", dfTagAndProb11[["L1_L0Calo_ECAL_region", "L2_L0Calo_ECAL_region",'L0Calo_ECAL_region_max']].head()
      #print "dfTagAndProb01:\n", dfTagAndProb11[["L1_L0Calo_ECAL_realET", "L2_L0Calo_ECAL_realET",'L0Calo_ECAL_max_realET']].head()
      #print "dfTagAndProb10:\n", dfTagAndProb11[["L1_L0Calo_ECAL_realET", "L2_L0Calo_ECAL_realET",'L0Calo_ECAL_max_realET']].head()
      #print "dfTagAndProb01:\n", dfTagAndProb01[["L1_L0Calo_ECAL_region", "L2_L0Calo_ECAL_region",'L0Calo_ECAL_region_max']].head()
      #print "dfTagAndProb10:\n", dfTagAndProb10[["L1_L0Calo_ECAL_region", "L2_L0Calo_ECAL_region",'L0Calo_ECAL_region_max']].head()


   
   
   #arrayTag = np.concatenate((arrayTag11, arrayTag10, arrayTag01), axis=0)
   arrayTag = np.concatenate((arrayTag11, arrayTag10 ,arrayTag01), axis=0)
   arrayTagAndProb = np.concatenate((arrayTagAndProb11, arrayTagAndProb10 ,arrayTagAndProb01), axis=0)

   return arrayTag, arrayTagAndProb



def TAP_L0E_maxET(df, selTag, VERB):
   '''
   selTagDict = {"B_L0Global_TIS" : (df.B_L0Global_TIS == 1),
                 "All" : ((df.B_L0Global_TIS == 1) | (df.B_L0Global_TIS == 0)),
                 "test" : ((df.L1_L0ElectronDecision_TOS == 1) | (df.L2_L0ElectronDecision_TOS == 0)),
}

   selTag = selTagDict[option]
   '''
   selTrig = ((df.L1_L0ElectronDecision_TOS==1) | (df.L2_L0ElectronDecision_TOS==1))
   #selTag = (df.B_L0Global_TIS == 1)
   selProb    = selTag &  selTrig

   
   ## L1_TOS==1 & L2_TOS==1
   dfTag = df[selTag][['L1_L0Calo_ECAL_region','L1_L0Calo_ECAL_realET','L2_L0Calo_ECAL_region','L2_L0Calo_ECAL_realET','weights']]
   dfTag['L0Calo_ECAL_region_max'] = dfTag[['L2_L0Calo_ECAL_realET','L1_L0Calo_ECAL_realET']].idxmax(axis=1)
   dfTag['L0Calo_ECAL_max_realET'] = dfTag[["L1_L0Calo_ECAL_realET", "L2_L0Calo_ECAL_realET"]].max(axis=1)
   dfTag.ix[dfTag.L0Calo_ECAL_region_max.isin(['L1_L0Calo_ECAL_realET']), 'L0Calo_ECAL_region_max'] = dfTag.ix[dfTag.L0Calo_ECAL_region_max.isin(['L1_L0Calo_ECAL_realET']), 'L1_L0Calo_ECAL_region']
   dfTag.ix[dfTag.L0Calo_ECAL_region_max.isin(['L2_L0Calo_ECAL_realET']), 'L0Calo_ECAL_region_max'] = dfTag.ix[dfTag.L0Calo_ECAL_region_max.isin(['L2_L0Calo_ECAL_realET']), 'L2_L0Calo_ECAL_region']



   dfTagAndProb = df[selProb][['L1_L0Calo_ECAL_region','L1_L0Calo_ECAL_realET','L2_L0Calo_ECAL_region','L2_L0Calo_ECAL_realET','weights']]
   dfTagAndProb['L0Calo_ECAL_max_realET'] = dfTagAndProb[["L1_L0Calo_ECAL_realET", "L2_L0Calo_ECAL_realET"]].max(axis=1)
   dfTagAndProb['L0Calo_ECAL_region_max'] = dfTagAndProb[['L2_L0Calo_ECAL_realET','L1_L0Calo_ECAL_realET']].idxmax(axis=1)
   dfTagAndProb.ix[dfTagAndProb.L0Calo_ECAL_region_max.isin(['L1_L0Calo_ECAL_realET']), 'L0Calo_ECAL_region_max'] = dfTagAndProb.ix[dfTagAndProb.L0Calo_ECAL_region_max.isin(['L1_L0Calo_ECAL_realET']), 'L1_L0Calo_ECAL_region']
   dfTagAndProb.ix[dfTagAndProb.L0Calo_ECAL_region_max.isin(['L2_L0Calo_ECAL_realET']), 'L0Calo_ECAL_region_max'] = dfTagAndProb.ix[dfTagAndProb.L0Calo_ECAL_region_max.isin(['L2_L0Calo_ECAL_realET']), 'L2_L0Calo_ECAL_region']

   arrayTag = dfTag[['L0Calo_ECAL_region_max','L0Calo_ECAL_max_realET','weights']].values
   arrayTagAndProb = dfTagAndProb[['L0Calo_ECAL_region_max','L0Calo_ECAL_max_realET','weights']].values


   print "=================== TAP_L0E_maxET ============================\n\n"
   print "TIS & {}\\\\".format(dfTag.shape[0])
   print "L0E1 or L0E2 & {}\\\\".format(dfTagAndProb.shape[0])


   if VERB:

      print '========================== Appending the maximum ET column  ==============================='  
      print "dfTag:\n", dfTag[["L1_L0Calo_ECAL_realET", "L2_L0Calo_ECAL_realET",'L0Calo_ECAL_max_realET']].head()                           
      print "dfTag:\n", dfTag[["L1_L0Calo_ECAL_region", "L2_L0Calo_ECAL_region",'L0Calo_ECAL_region_max']].head()  
      print "dfTagAndProb:\n", dfTagAndProb[["L1_L0Calo_ECAL_realET", "L2_L0Calo_ECAL_realET",'L0Calo_ECAL_max_realET']].head() 
      print "dfTagAndProb:\n", dfTagAndProb[["L1_L0Calo_ECAL_region", "L2_L0Calo_ECAL_region",'L0Calo_ECAL_region_max']].head()      

   return arrayTag, arrayTagAndProb

def TAP_L0E_indipE(df, selTag, VERB, particle='both'):

   
   if (particle == 'e+'):    particle_list = ['e+']
   if (particle == 'e-'):    particle_list = ['e-']
   if (particle == 'both'):    particle_list = ['e+','e-']

   
   arrayTag =  np.empty([0,3])
   arrayTagAndProb =  np.empty([0,3])



   for ipart in particle_list:
               
      #Distinguish between e+ and e-
      if ipart == 'e-' :  
      

         selcharge1 = (df.L1_ID <0)
         selcharge2 = (df.L2_ID <0)
      
      
      elif ipart == 'e+':
         selcharge1 = (df.L1_ID >0)
         selcharge2 = (df.L2_ID >0)
      
      
      #Number of L1 corresponding to i.e. e- that are triggered by L0Electron   
      selProb1 = (df.L1_L0ElectronDecision_TOS==1) & selcharge1
      #Number of i.e. e- that are selected thanks of B_L0Global_TIS
      selTag1 = selTag & selcharge1
      #Number of i.e. e- that are selected thanks of B_L0Global_TIS and L0Electron
      selTagAndProb1    = selTag1 &  selProb1
         
      #Number of L2 corresponding to i.e. e- that are triggered by L0Electron   
      selProb2 = (df.L2_L0ElectronDecision_TOS==1) & selcharge2
      #Number of i.e. e- that are selected thanks of B_L0Global_TIS
      selTag2 = selTag & selcharge2
      #Number of i.e. e- that are selected thanks of B_L0Global_TIS and L0Electron
      selTagAndProb2    = selTag2 &  selProb2
      
      #Create np.array to fill TH1F
      #dfselTag1 =  df[selTag1][['L1_L0Calo_ECAL_region','L1_L0Calo_ECAL_realET','wL01',"L1_L0ElectronDecision_TOS","B_L0Global_TIS"]]
      #dfselTagAndProb1 = df[selTagAndProb1][['L1_L0Calo_ECAL_region','L1_L0Calo_ECAL_realET','wL01',"L1_L0ElectronDecision_TOS","B_L0Global_TIS"]]
      arrayTag1 = df[selTag1][['L1_L0Calo_ECAL_region','L1_L0Calo_ECAL_realET','weights']].values
      arrayTagAndProb1 = df[selTagAndProb1][['L1_L0Calo_ECAL_region','L1_L0Calo_ECAL_realET','weights']].values
      
      #dfselTag2 = df[selTag2][['L2_L0Calo_ECAL_region','L2_L0Calo_ECAL_realET','wL02',"L2_L0ElectronDecision_TOS","B_L0Global_TIS"]]
      #dfselTagAndProb2 = df[selTagAndProb2][['L2_L0Calo_ECAL_region','L2_L0Calo_ECAL_realET','wL02',"L2_L0ElectronDecision_TOS","B_L0Global_TIS"]]
      arrayTag2 = df[selTag2][['L2_L0Calo_ECAL_region','L2_L0Calo_ECAL_realET','weights']].values
      arrayTagAndProb2 = df[selTagAndProb2][['L2_L0Calo_ECAL_region','L2_L0Calo_ECAL_realET','weights']].values



      arrayTag_t = np.concatenate((arrayTag1, arrayTag2), axis=0)
      arrayTagAndProb_t = np.concatenate((arrayTagAndProb1, arrayTagAndProb2), axis=0)

      #It is strange to have two electrons in the same decay if only the signal was selected
      # Verify that the L1 L2 samples have zero intersection
      if (arrayTag_t.shape[0] == df[selTag].shape[0]):
         print "Consistency check passed..."
      else:
         print "WARNING! Consistency check NOT PASSED!"
         quit


      #Concatenate e+ and e- if the option both was chosen, otherwise in just renaming
      arrayTag = np.concatenate((arrayTag, arrayTag_t), axis=0)
      arrayTagAndProb = np.concatenate((arrayTagAndProb, arrayTagAndProb_t),axis=0)




   return arrayTag, arrayTagAndProb
#################

def TAP_L0M_indipM(df, selTag, VERB, particle='both'):

   
   if (particle == 'mu+'):    particle_list = ['mu+']
   if (particle == 'mu-'):    particle_list = ['mu-']
   if (particle == 'both'):    particle_list = ['mu+','mu-']

   
   arrayTag =  np.empty([0,2])
   arrayTagAndProb =  np.empty([0,2])



   for ipart in particle_list:
               
      #Distinguish between e+ and e-
      if ipart == 'mu-' :  
      

         selcharge1 = (df.L1_ID <0)
         selcharge2 = (df.L2_ID <0)
      
      
      elif ipart == 'mu+':
         selcharge1 = (df.L1_ID >0)
         selcharge2 = (df.L2_ID >0)
      
      
      #Number of L1 corresponding to i.e. e- that are triggered by L0Muon   
      selProb1 = (df.L1_L0MuonDecision_TOS==1) & selcharge1
      #Number of i.e. e- that are selected thanks of B_L0Global_TIS
      selTag1 = selTag & selcharge1
      #Number of i.e. e- that are selected thanks of B_L0Global_TIS and L0Muon
      selTagAndProb1    = selTag1 &  selProb1
         
      #Number of L2 corresponding to i.e. e- that are triggered by L0Muon   
      selProb2 = (df.L2_L0MuonDecision_TOS==1) & selcharge2
      #Number of i.e. e- that are selected thanks of B_L0Global_TIS
      selTag2 = selTag & selcharge2
      #Number of i.e. e- that are selected thanks of B_L0Global_TIS and L0Muon
      selTagAndProb2    = selTag2 &  selProb2


      #Create np.array to fill TH1F
      #dfselTag1 =  df[selTag1][['L1_L0Calo_ECAL_region','L1_PT','wL01',"L1_L0MuonDecision_TOS","B_L0Global_TIS"]]
      #dfselTagAndProb1 = df[selTagAndProb1][['L1_L0Calo_ECAL_region','L1_PT','wL01',"L1_L0MuonDecision_TOS","B_L0Global_TIS"]]
      arrayTag1 = df[selTag1][['L1_PT','weights']].values
      arrayTagAndProb1 = df[selTagAndProb1][['L1_PT','weights']].values
      
      #dfselTag2 = df[selTag2][['L2_L0Calo_ECAL_region','L2_PT','wL02',"L2_L0MuonDecision_TOS","B_L0Global_TIS"]]
      #dfselTagAndProb2 = df[selTagAndProb2][['L2_L0Calo_ECAL_region','L2_PT','wL02',"L2_L0MuonDecision_TOS","B_L0Global_TIS"]]
      arrayTag2 = df[selTag2][['L2_PT','weights']].values
      arrayTagAndProb2 = df[selTagAndProb2][['L2_PT','weights']].values



      arrayTag_t = np.concatenate((arrayTag1, arrayTag2), axis=0)
      arrayTagAndProb_t = np.concatenate((arrayTagAndProb1, arrayTagAndProb2), axis=0)

      #It is strange to have two electrons in the same decay if only the signal was selected
      # Verify that the L1 L2 samples have zero intersection
      if (arrayTag_t.shape[0] == df[selTag].shape[0]):
         print "Consistency check passed..."
      else:
         print "WARNING! Consistency check NOT PASSED!"
         quit


      #Concatenate e+ and e- if the option both was chosen, otherwise in just renaming
      arrayTag = np.concatenate((arrayTag, arrayTag_t), axis=0)
      arrayTagAndProb = np.concatenate((arrayTagAndProb, arrayTagAndProb_t),axis=0)




   return arrayTag, arrayTagAndProb
############


def TAP_L0E_indipendentE(df, selTag, VERB, particle='both'):

   
   if (particle == 'e+'):    particle_list = ['e+']
   if (particle == 'e-'):    particle_list = ['e-']
   if (particle == 'both'):    particle_list = ['e+','e-']

   
   arrayTag =  np.empty([0,3])
   arrayTagAndProb =  np.empty([0,3])



   for ipart in particle_list:
               
      #Distinguish between e+ and e-
      if ipart == 'e-' :  
      

         selcharge1 = (df.L1_ID <0)
         selcharge2 = (df.L2_ID <0)
      
      
      elif ipart == 'e+':
         selcharge1 = (df.L1_ID >0)
         selcharge2 = (df.L2_ID >0)
      
      
      #Number of L1 corresponding to i.e. e- that are triggered by L0Electron   
      selProb1 = (df.L1_L0ElectronDecision_TOS==1) & selcharge1
      #Number of i.e. e- that are selected thanks of B_L0Global_TIS
      selTag1 = selTag & selcharge1
      #Number of i.e. e- that are selected thanks of B_L0Global_TIS and L0Electron
      selTagAndProb1    = selTag1 &  selProb1
         
      #Number of L2 corresponding to i.e. e- that are triggered by L0Electron   
      selProb2 = (df.L2_L0ElectronDecision_TOS==1) & selcharge2
      #Number of i.e. e- that are selected thanks of B_L0Global_TIS
      selTag2 = selTag & selcharge2
      #Number of i.e. e- that are selected thanks of B_L0Global_TIS and L0Electron
      selTagAndProb2    = selTag2 &  selProb2
      
      #Create np.array to fill TH1F
      #dfselTag1 =  df[selTag1][['L1_L0Calo_ECAL_region','L1_L0Calo_ECAL_realET','wL01',"L1_L0ElectronDecision_TOS","B_L0Global_TIS"]]
      #dfselTagAndProb1 = df[selTagAndProb1][['L1_L0Calo_ECAL_region','L1_L0Calo_ECAL_realET','wL01',"L1_L0ElectronDecision_TOS","B_L0Global_TIS"]]
      arrayTag1 = df[selTag1][['L1_L0Calo_ECAL_region','L1_L0Calo_ECAL_realET','wL01']].values
      arrayTagAndProb1 = df[selTagAndProb1][['L1_L0Calo_ECAL_region','L1_L0Calo_ECAL_realET','wL01']].values
      
      #dfselTag2 = df[selTag2][['L2_L0Calo_ECAL_region','L2_L0Calo_ECAL_realET','wL02',"L2_L0ElectronDecision_TOS","B_L0Global_TIS"]]
      #dfselTagAndProb2 = df[selTagAndProb2][['L2_L0Calo_ECAL_region','L2_L0Calo_ECAL_realET','wL02',"L2_L0ElectronDecision_TOS","B_L0Global_TIS"]]
      arrayTag2 = df[selTag2][['L2_L0Calo_ECAL_region','L2_L0Calo_ECAL_realET','wL02']].values
      arrayTagAndProb2 = df[selTagAndProb2][['L2_L0Calo_ECAL_region','L2_L0Calo_ECAL_realET','wL02']].values



      arrayTag_t = np.concatenate((arrayTag1, arrayTag2), axis=0)
      arrayTagAndProb_t = np.concatenate((arrayTagAndProb1, arrayTagAndProb2), axis=0)

      #It is strange to have two electrons in the same decay if only the signal was selected
      # Verify that the L1 L2 samples have zero intersection
      if (arrayTag_t.shape[0] == df[df.B_L0Global_TIS==1].shape[0]):
         print "Consistency check passed..."
      else:
         print "WARNING! Consistency check NOT PASSED!"
         quit


      #Concatenate e+ and e- if the option both was chosen, otherwise in just renaming
      arrayTag = np.concatenate((arrayTag, arrayTag_t), axis=0)
      arrayTagAndProb = np.concatenate((arrayTagAndProb, arrayTagAndProb_t),axis=0)




   return arrayTag, arrayTagAndProb


def TAP_L0M_maxpt(df, selTag) :    # Open calibration dataset    

   '''
   Author: Michele Atzeni
   Email Address: michele.atzeni@cern.ch
   Date: June 7th, 2017
   Note: script based on T. Humair's B+ -> J/psi K+ algorithm and R.Coutinho

   Description:
   This script computes the trigger-on-signal efficiency for mu+ and/or mu- as a function of the trasversal momentum PT using a Tag&Probe approach.

   Skeleton: 1)Divide the sample in two: i.e. L1_ID<0 and L2_ID<0
             2) Select the events respect to the Tag criteria B_L0Global_TIS==1
             3) Select the events respect to the Tag&Probe criteria (B_L0Global_TIS==1) & (L#_L0Muon_TOS==1)
             4)Plot the ratio as a function of the PT of the muon

   '''

   # Defining the triggered sample (the one we are interested in), the Tag sample and the intersection of the two
   
   selTrig = (df.L1_L0MuonDecision_TOS==1) | (df.L2_L0MuonDecision_TOS==1)
   #selTag = (df.B_L0Global_TIS == 1)
   selProb    = selTag &  selTrig

   
   ## L1_TOS==1 & L2_TOS==1
   dfTag = df[selTag][['L1_PT','L2_PT','weights']]
   dfTag['PT_max'] = dfTag[["L1_PT", "L2_PT"]].max(axis=1)


   dfTagAndProb = df[selProb][['L1_PT','L2_PT','weights']]
   dfTagAndProb['PT_max'] = dfTagAndProb[["L1_PT", "L2_PT"]].max(axis=1)

   arrayTag = dfTag[['PT_max','weights']].values
   arrayTagAndProb = dfTagAndProb[['PT_max','weights']].values


   return arrayTag, arrayTagAndProb
   ########################


def TAP_L0M2_mixedmaxPT(df, selTag, VERB):




   # Defining the triggered sample (the one we are interested in), the Tag sample and the intersection of the two
   #selTrig = (df.K_L0HadronDecision_TOS==1) & (df.Pi_L0HadronDecision_TOS==1)
   selTrig11 = ((df.L1_L0MuonDecision_TOS==1) & (df.L2_L0MuonDecision_TOS==1))
   selTrig10 = ((df.L1_L0MuonDecision_TOS==1) & (df.L2_L0MuonDecision_TOS==0))
   selTrig01 = ((df.L1_L0MuonDecision_TOS==0) & (df.L2_L0MuonDecision_TOS==1))

   #selTag = (df.B_L0Global_TIS == 1)
   selProb11    = selTag &  selTrig11
   selProb10    = selTag &  selTrig10
   selProb01    = selTag &  selTrig01
   
   ## L1_TOS==1 & L2_TOS==1
   dfTag11 = df[selTag][['L1_PT','L2_PT','weights']]
   dfTag11['PT_max'] = dfTag11[["L1_PT", "L2_PT"]].max(axis=1)

   dfTagAndProb11 = df[selProb11][['L1_PT','L2_PT','weights']]
   dfTagAndProb11['PT_max'] = dfTagAndProb11[["L1_PT", "L2_PT"]].max(axis=1)

   arrayTag11 = dfTag11[['PT_max','weights']].values
   arrayTagAndProb11 = dfTagAndProb11[['PT_max','weights']].values

   ## L1_TOS==1 & L2_TOS==0
   dfTag10 = df[selTag]
   arrayTag10 = dfTag10[['L1_PT','weights']].values
   dfTagAndProb10 = df[selProb10]
   arrayTagAndProb10 = df[selProb10][['L1_PT','weights']].values

   dfTag01 = df[selTag]
   arrayTag01 = dfTag01[['L2_PT','weights']].values
   dfTagAndProb01 = df[selProb01]
   arrayTagAndProb01 = df[selProb01][['L2_PT','weights']].values
   

   print "=================== TAP_L0M2_mixedmaxPT ============================\n\n"
   print "TIS & {},{},{}\\\\".format(dfTag11.shape[0],dfTag01.shape[0],dfTag10.shape[0])
   print "L0M1 and L0M2 & {}\\\\".format(dfTagAndProb11.shape[0])
   print "L0M1 and !L0M2 & {}\\\\".format(dfTagAndProb10.shape[0])
   print "!L0M1 and L0M2 & {}\\\\".format(dfTagAndProb01.shape[0])

   if VERB:

      print '========================== Appending the maximum PT column  ==============================='
      print "dfTag:\n", dfTag11[["L1_PT", "L2_PT",'PT_max']].head()
      print "dfTagAndProb11:\n", dfTagAndProb11[["L1_PT", "L2_PT",'PT_max']].head()
      print "dfTagAndProb11:\n", dfTagAndProb11[["L1_PT", "L2_PT",'PT_max']].head()
      #print "dfTagAndProb01:\n", dfTagAndProb11[["L1_L0Calo_ECAL_realET", "L2_L0Calo_ECAL_realET",'L0Calo_ECAL_max_realET']].head()
      #print "dfTagAndProb10:\n", dfTagAndProb11[["L1_L0Calo_ECAL_realET", "L2_L0Calo_ECAL_realET",'L0Calo_ECAL_max_realET']].head()
      #print "dfTagAndProb01:\n", dfTagAndProb01[["L1_L0Calo_ECAL_region", "L2_L0Calo_ECAL_region",'L0Calo_ECAL_region_max']].head()
      #print "dfTagAndProb10:\n", dfTagAndProb10[["L1_L0Calo_ECAL_region", "L2_L0Calo_ECAL_region",'L0Calo_ECAL_region_max']].head()


   
   
   #arrayTag = np.concatenate((arrayTag11, arrayTag10, arrayTag01), axis=0)
   arrayTag = np.concatenate((arrayTag11, arrayTag10 ,arrayTag01), axis=0)
   arrayTagAndProb = np.concatenate((arrayTagAndProb11, arrayTagAndProb10 ,arrayTagAndProb01), axis=0)

   return arrayTag, arrayTagAndProb



def TAP_L0H_KstarPT_notL0E(df, selTag):    # Open calibration dataset    

   '''
   Author: Michele Atzeni
   Email Address: michele.atzeni@cern.ch
   '''
   # Defining the triggered sample (the one we are interested in), the Tag sample and the intersection of the two

   #selTrig = ((df.Kstar_L0HadronDecision_TOS==1) & (df.L1_L0ElectronDecision_TOS==0) & (df.L2_L0ElectronDecision_TOS==0) & (df.Kstar_PT > 4000.))
   selTrig = ((df.Kstar_L0HadronDecision_TOS==1) & (df.L1_L0ElectronDecision_TOS==0) & (df.L2_L0ElectronDecision_TOS==0))

   selProb    = selTag &  selTrig
   
   
   arrayTag = df[selTag][['K_L0Calo_HCAL_region','Pi_L0Calo_HCAL_region','Kstar_PT',"weights"]].values
   arrayTagAndProb = df[selProb][['K_L0Calo_HCAL_region','Pi_L0Calo_HCAL_region','Kstar_PT',"weights"]].values

   return  arrayTag, arrayTagAndProb

   #########################

def TAP_L0H_KstarPT_notL0M(df, selTag):    # Open calibration dataset    

   '''
   Author: Michele Atzeni
   Email Address: michele.atzeni@cern.ch
   '''
   # Defining the triggered sample (the one we are interested in), the Tag sample and the intersection of the two
   
   #selTrig = ((df.Kstar_L0HadronDecision_TOS==1) & (df.L1_L0MuonDecision_TOS==0) & (df.L2_L0MuonDecision_TOS==0) & (df.Kstar_PT > 4000.))
   selTrig = ((df.Kstar_L0HadronDecision_TOS==1) & (df.L1_L0MuonDecision_TOS==0) & (df.L2_L0MuonDecision_TOS==0))
   
   selProb    = selTag &  selTrig
   
   
   arrayTag = df[selTag][['K_L0Calo_HCAL_region','Pi_L0Calo_HCAL_region','Kstar_PT',"weights"]].values
   arrayTagAndProb = df[selProb][['K_L0Calo_HCAL_region','Pi_L0Calo_HCAL_region','Kstar_PT',"weights"]].values

   return  arrayTag, arrayTagAndProb

   ######################

def TAP_L0H_KstarPT_alsoL0L(df, selTag):    # Open calibration dataset    

   '''
   Author: Michele Atzeni
   Email Address: michele.atzeni@cern.ch
   '''
   # Defining the triggered sample (the one we are interested in), the Tag sample and the intersection of the two
   
   selTrig = (df.Kstar_L0HadronDecision_TOS==1)
   #selTrig = ((df.Kstar_L0HadronDecision_TOS==1) & (df.Kstar_PT > 4000.))
   selProb    = selTag &  selTrig
   
   
   arrayTag = df[selTag][['K_L0Calo_HCAL_region','Pi_L0Calo_HCAL_region','Kstar_PT',"weights"]].values
   arrayTagAndProb = df[selProb][['K_L0Calo_HCAL_region','Pi_L0Calo_HCAL_region','Kstar_PT',"weights"]].values

   return  arrayTag, arrayTagAndProb


   ############################

def TAP_L0TIS_alsoL0LH(df, selTag):    # Open calibration dataset    

   '''
   
   '''

   # Defining the triggered sample (the one we are interested in), the Tag sample and the intersection of the two
   selTrig = (df.B_L0Global_TIS==1) #& ((df.Kstar_L0HadronDecision_TOS==0) & (df.L1_L0ElectronDecision_TOS==0) & (df.L2_L0ElectronDecision_TOS==0))
   #selTrig = (df.B_L0Global_TIS==1) & ((df.Kstar_L0HadronDecision_TOS==0) & (df.L1_L0ElectronDecision_TOS==0) & (df.L2_L0ElectronDecision_TOS==0))
   #selTag = (df.B_L0Global_TOS == 1)
   selProb    = selTag &  selTrig
   
   
   arrayTag = df[selTag][['nSPDHits','B_PT','weights']].values
   arrayTagAndProb = df[selProb][['nSPDHits','B_PT','weights']].values

   return arrayTag, arrayTagAndProb

def TAP_L0TIS_notL0MH(df, selTag):    # Open calibration dataset    

   '''
   
   '''

   # Defining the triggered sample (the one we are interested in), the Tag sample and the intersection of the two
   selTrig = (df.B_L0Global_TIS==1) & ((df.Kstar_L0HadronDecision_TOS==0) & (df.L1_L0MuonDecision_TOS==0) & (df.L2_L0MuonDecision_TOS==0))
   selProb    = selTag &  selTrig
   
   
   arrayTag = df[selTag][['nSPDHits','B_PT','weights']].values
   arrayTagAndProb = df[selProb][['nSPDHits','B_PT','weights']].values

   return arrayTag, arrayTagAndProb


def TAP_L0TIS_notL0EH(df, selTag):    # Open calibration dataset    

   '''
   
   '''

   # Defining the triggered sample (the one we are interested in), the Tag sample and the intersection of the two
   
   selTrig = (df.B_L0Global_TIS==1) & ((df.Kstar_L0HadronDecision_TOS==0) & (df.L1_L0ElectronDecision_TOS==0) & (df.L2_L0ElectronDecision_TOS==0))
   selProb    = selTag &  selTrig
   
   
   arrayTag = df[selTag][['nSPDHits','B_PT','weights']].values
   arrayTagAndProb = df[selProb][['nSPDHits','B_PT','weights']].values

   return arrayTag, arrayTagAndProb


def TAP_HLT_E(df, year, VERB) :    # Open calibration dataset    

   '''
   Author: Michele Atzeni
   Email Address: michele.atzeni@cern.ch
   Date: June 7th, 2017
   Note: script based on T. Humair's B+ -> J/psi K+ algorithm and R.Coutinho

   Description:
   This script computes the trigger-on-signal efficiency for mu+ and/or mu- as a function of the trasversal momentum PT using a Tag&Probe approach.

   Skeleton  1)Divide the sample in two: i.e. L1_ID<0 and L2_ID<0
             2) Select the events respect to the Tag criteria B_L0Global_TIS==1
             3) Select the events respect to the Tag&Probe criteria (B_L0Global_TIS==1) & (L#_L0Muon_TOS==1)
             4)Plot the ratio as a function of the PT of the muon

   '''

   # Defining the triggered sample (the one we are interested in), the Tag sample and the intersection of the two
   selTag = (df.B_Hlt1Phys_TIS ==1 ) & (df.B_Hlt2Phys_TIS ==1)
   #selTag = (df.B_L0Global_TIS==1)
   # The HLT1 lines are different for 2011-2012 and 2015-2016.                                                                                   
   if ((year == '11') or (year == '12') or (year == 'RunI')) :
      cutHLT1 = (df.B_Hlt1TrackAllL0Decision_TOS==1)
   elif ((year == '15') or (year == '16')or (year == 'RunII')) :
      cutHLT1 = ((df.B_Hlt1TrackMVADecision_TOS==1) | (df.B_Hlt1TwoTrackMVADecision_TOS==1) | (df.B_Hlt1ElectronTrackDecision_TOS==1) | (df.B_Hlt1TrackMVALooseDecision_TOS==1) | (df.B_Hlt1TwoTrackMVALooseDecision_TOS==1))
   # The HLT2 lines are different for 2011-2012 and 2015-2016.                                                                                        
   if ((year == '11') or (year == '12') or (year == 'RunI')) :
      cutHLT2 = ((df.B_Hlt2Topo2BodyBBDTDecision_TOS==1) | (df.B_Hlt2Topo3BodyBBDTDecision_TOS==1) | (df.B_Hlt2Topo4BodyBBDTDecision_TOS==1) | (df.B_Hlt2TopoE2BodyBBDTDecision_TOS==1) | (df.B_Hlt2TopoE3BodyBBDTDecision_TOS==1) | (df.B_Hlt2TopoE4BodyBBDTDecision_TOS==1))
   elif ((year == '15') or (year == '16')or (year == 'RunII')) :
      cutHLT2 = ((df.B_Hlt2Topo2BodyDecision_TOS==1) | (df.B_Hlt2Topo3BodyDecision_TOS==1) | (df.B_Hlt2Topo4BodyDecision_TOS==1) | (df.B_Hlt2TopoE2BodyDecision_TOS==1) | (df.B_Hlt2TopoE3BodyDecision_TOS==1) | (df.B_Hlt2TopoE4BodyDecision_TOS==1) | (df.B_Hlt2TopoEE2BodyDecision_TOS==1) | (df.B_Hlt2TopoEE3BodyDecision_TOS==1) | (df.B_Hlt2TopoEE4BodyDecision_TOS==1))

   selTrig = cutHLT1 & cutHLT2
   selProb    = selTag &  selTrig

   
   ## L1_TOS==1 & L2_TOS==1
   dfTag = df[selTag][['L1_PT','L2_PT','K_PT','Pi_PT','weights']]
   dfTag['PT_min'] = dfTag[['K_PT','Pi_PT',"L1_PT", "L2_PT"]].min(axis=1)


   dfTagAndProb = df[selProb][['L1_PT','L2_PT','K_PT','Pi_PT','weights']]
   dfTagAndProb['PT_min'] = dfTagAndProb[["L1_PT", "L2_PT",'K_PT','Pi_PT']].min(axis=1)

   arrayTag = dfTag[['PT_min','weights']].values
   arrayTagAndProb = dfTagAndProb[['PT_min','weights']].values


   return arrayTag, arrayTagAndProb
   ########################


def TAP_HLT_M(df, year, VERB) :    # Open calibration dataset    

   '''
   Author: Michele Atzeni
   Email Address: michele.atzeni@cern.ch
   Date: June 7th, 2017
   Note: script based on T. Humair's B+ -> J/psi K+ algorithm and R.Coutinho

   Description:
   This script computes the trigger-on-signal efficiency for mu+ and/or mu- as a function of the trasversal momentum PT using a Tag&Probe approach.

   Skeleton  1)Divide the sample in two: i.e. L1_ID<0 and L2_ID<0
             2) Select the events respect to the Tag criteria B_L0Global_TIS==1
             3) Select the events respect to the Tag&Probe criteria (B_L0Global_TIS==1) & (L#_L0Muon_TOS==1)
             4)Plot the ratio as a function of the PT of the muon

   '''

   # Defining the triggered sample (the one we are interested in), the Tag sample and the intersection of the two
#   selTag = (df.B_Hlt1Phys_TIS ==1 ) & (df.B_Hlt2Phys_TIS ==1)

   selTag = (df.B_Hlt1Phys_TIS ==1 ) & (df.B_Hlt2Phys_TIS ==1)
   cutHLT1 = ((df.B_Hlt1TrackAllL0Decision_TOS==1) | (df.B_Hlt1TrackMuonDecision_TOS==1))
   cutHLT2 = ((df.B_Hlt2Topo2BodyBBDTDecision_TOS==1) | (df.B_Hlt2Topo3BodyBBDTDecision_TOS==1) | (df.B_Hlt2Topo4BodyBBDTDecision_TOS==1) | (df.B_Hlt2TopoMu2BodyBBDTDecision_TOS==1) | (df.B_Hlt2TopoMu3BodyBBDTDecision_TOS==1) | (df.B_Hlt2TopoMu4BodyBBDTDecision_TOS==1) | (df.B_Hlt2DiMuonDetachedDecision_TOS==1))

   selTrig = cutHLT1 & cutHLT2
   selProb    = selTag &  selTrig

   
   ## L1_TOS==1 & L2_TOS==1
   dfTag = df[selTag][['L1_PT','L2_PT','K_PT','Pi_PT','weights']]
   dfTag['PT_min'] = dfTag[['K_PT','Pi_PT',"L1_PT", "L2_PT"]].min(axis=1)


   dfTagAndProb = df[selProb][['L1_PT','L2_PT','K_PT','Pi_PT','weights']]
   dfTagAndProb['PT_min'] = dfTagAndProb[["L1_PT", "L2_PT",'K_PT','Pi_PT']].min(axis=1)

   arrayTag = dfTag[['PT_min','weights']].values
   arrayTagAndProb = dfTagAndProb[['PT_min','weights']].values


   return arrayTag, arrayTagAndProb
   ########################



if __name__ == "__main__" :

   parser = argparse.ArgumentParser(description = 'Configuration of the parameters for the MoM calculation')
   

   parser.add_argument("-df", "--DataFrame" , dest="df"  , required=True,  help="Choose the DF to use")
   parser.add_argument("-p", "--particleName" , dest="particle"  , required=True,  help="Choose the name of the particle to be evaluated (i.e. K or pi)")
   parser.add_argument("-t", "--Taginfo" , dest="Tag"  , required=True,  help="Tag used to distinguish the efficiency plot")


   args = parser.parse_args()

   # Parameters and configuration
   df      = args.df
   particle  = args.particle 
   Tag = args.Tag

   

   TagAndProbe_e(df, particle, Tag)
   
