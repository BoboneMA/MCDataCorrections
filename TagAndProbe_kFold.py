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
import os



import argparse

import ROOT
from ROOT import gROOT, gStyle, TFile, TH2F , TH1F,  TCanvas , TLegend, kTRUE
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
import seaborn as sns
from Adaptive_binning_pro import AdaptiveBinning1D
import pdb
from TagAndProbe_core import TAP_L0E_maxET, TAP_L0E2_mixedmaxET, TAP_L0E_indipE, TAP_L0M_indipM, TAP_L0H_KstarPT_notL0E, TAP_L0H_KstarPT_alsoL0L, TAP_L0M_maxpt, TAP_L0M2_mixedmaxPT, TAP_L0H_KstarPT_notL0M, TAP_L0TIS_alsoL0LH, TAP_L0TIS_notL0EH, TAP_L0TIS_notL0MH, TAP_HLT_E, TAP_HLT_M



#Before using this script be  sure all the variables used here are in the dataframe
'''

'''

def TagAndProbe_L0E(df1, selTag, Model,  Tag, Ada, weight, VERB) :    # Open calibration dataset    

   '''
   Author: Michele Atzeni
   Email Address: michele.atzeni@cern.ch
   Date: June 7th, 2017


   Description:
   This script computes the trigger-on-signal efficiency for electrons and/or positrons as a function of the trasversal energy deposited in the ECAL using a Tag&Probe approach.
   The electrons are divided depending on the region of the calorimeter hit: inner, medium or outer.

   Concept:
   Ideally to calculate a trigger efficiency you would just do the following:
   
   N_pt = Number of electrons that passed the trigger
   N_t  = Number of total electrons

   e_trigger = N_pt/N_t

   However the quantity N_t is not accessible experimentally.
   We can have a datadriven estimate of the trigger efficiency using a  Tag&Probe approach.

   The idea is very simple: we look at the ratio e_trigger for all those events that passed a specific trigger condition, called TAGGING. If the TAGGING trigger  samples similarly the set of all the electrons and the set of triggered ones the two efficiencies must be the same. 

   e_trigger_TAG = N_pt_tag/ N_t_tag =>e_trigger

   Where:  + N_pt_tag = Number of electrons that are TAGged with B_L0Global_TIS and that are triggered as L#_L0ElectronDecision_TOS(PROBE)
           + N_t_tag = Number of electrons that are TAGged with B_L0Global_TIS 
   

   Scheme:
   1)Get the full dataframe
   2)Choose the charge of the particle: charge =(L1_ID >0) 
                                                (L2_ID >0)

   3)Define the Probe and Tag selections: TAG    => B_L0Global_TIS==1 & L1_ID >0
                                                 => B_L0Global_TIS==1 & L2_ID >0

                                          PROBE  => L1_L0ElectronDecision_TOS==1 & L1_ID >0
                                                 => L2_L0ElectronDecision_TOS==1 & L2_ID >0

   4)Select the dataframe and put together the L1 and L2 
   5)Get the arrays for the TAGged  and the TAG&Probe sample
   6)Calculate efficiency

   '''


   print "================ Tag&Probe algorithm: L0E ==================="
   print "=================for {}===============".format(Tag)

   if(not weight): 


      branchesTag = ["B_PVandJpsiDTF_B_M","L1_L0ElectronDecision_TOS","B_L0Global_TIS","L1_ID","L2_ID","L2_L0ElectronDecision_TOS","L1_L0Calo_ECAL_region","L2_L0Calo_ECAL_region","L1_L0Calo_ECAL_realET","L2_L0Calo_ECAL_realET"]
      #Reduce the dataframe to speed up
      df = df1[branchesTag]
      df.insert(1, 'weights',1)
      
      if VERB:
         print "********* In the following a uniform weight (1) is used for all the calculations *********"
         print(df[["L1_L0ElectronDecision_TOS","weights"]].head())
         print "******************************************************************************************"

   else:
      if(weight not in df1.columns):
      
         branchesTag = ["B_PVandJpsiDTF_B_M","L1_L0ElectronDecision_TOS","B_L0Global_TIS","L1_ID","L2_ID","L2_L0ElectronDecision_TOS","L1_L0Calo_ECAL_region","L2_L0Calo_ECAL_region","L1_L0Calo_ECAL_realET","L2_L0Calo_ECAL_realET"]
         #Reduce the dataframe to speed up
         df = df1[branchesTag]
         df.insert(1, 'weights',1)
         print "****************************************=WARNING******************************************************="
         print "********= Was requested a weighted histogram, however no column found with the name 'weights' ********="
         print "************** In the following a uniform weight (1) is used for all the calculations ****************="
         print(df[["L1_L0ElectronDecision_TOS","weights"]].head())
         print "******************************************************************************************************="
         exit()
      else:

         branchesTag = ["B_PVandJpsiDTF_B_M","L1_L0ElectronDecision_TOS","B_L0Global_TIS","L1_ID","L2_ID","L2_L0ElectronDecision_TOS","L1_L0Calo_ECAL_region","L2_L0Calo_ECAL_region","L1_L0Calo_ECAL_realET","L2_L0Calo_ECAL_realET"]
         branchesTag.append(weight)
         #Reduce the dataframe to speed up
         df = df1[branchesTag]
         df = df.rename(columns={weight:'weights'})
         if VERB:
            print "************ Column of weights successfully loaded **************"
            print(df[["L1_L0ElectronDecision_TOS","weights"]].head())
            print "****************************************************************="
         
   # Defining the triggered sample (the one we are interested in), the Tag sample and the intersection of the two
   
   print "********************************************************************************"
   print "****** Number of events in the dataframe used: {}********".format(df.shape[0])
   print "********************************************************************************"


   selTagDict = {"B_L0Global_TIS" : (df.B_L0Global_TIS == 1),
                 "All" : ((df.B_L0Global_TIS == 1) | (df.B_L0Global_TIS == 0)),
                 "test" : ((df.L1_L0ElectronDecision_TOS == 1) | (df.L2_L0ElectronDecision_TOS == 1)),
                 #"All_TIS" : ((df.B_L0MuonDecision_TIS == 1) | (df.B_L0DiMuonDecision_TIS == 1) | (df.B_L0ElectronDecision_TIS == 1) | (df.B_L0PhotonDecision_TIS == 1) | (df.B_L0HadronDecision_TIS == 1)),
                 #"EH_TIS" : ((df.B_L0ElectronDecision_TIS == 1) | (df.B_L0HadronDecision_TIS == 1)),
                 #"H_TIS" : (df.B_L0HadronDecision_TIS == 1),

}

   print "********** TIS sample used: {} ***********".format(selTag)

   if (Model == "maxET"):
      arrayTag, arrayTagAndProb = TAP_L0E_maxET(df, selTagDict[selTag], VERB)
   elif (Model == "mixedmaxET"):
      arrayTag, arrayTagAndProb = TAP_L0E2_mixedmaxET(df, selTagDict[selTag],VERB) 
   elif(Model == "indipE"):
      arrayTag, arrayTagAndProb = TAP_L0E_indipE(df, selTagDict[selTag], VERB, )
   else:
      print "WARNING!No model Found!"
      exit()
   #arrayTag, arrayTagAndProb = TAP_L0E_indipE(df,'both')



   print "********************************************************************************"
   print "****** Number of events in the TIS Sample: {}********".format(len(arrayTag))
   print "********************************************************************************"
   print "********************************************************************************"
   print "****** Number of events in the TOS Sample: {}********".format(len(arrayTagAndProb))
   print "********************************************************************************"
 
  ######################################################################################            


   #print 'The ProbandTag efficiency for the L0Electron trigger in Data is ', float(len(arrayTagAndProb))/float(len(arrayTag))

   #
   locationBinTab = [-1.5, -0.5, 0.5, 1.5, 2.5]
   regRange = np.empty(len(locationBinTab),dtype=np.float32)
   for index, iRegion in enumerate(locationBinTab): 
      regRange[index] = iRegion
   #
   weights_a = np.arange(0,4,0.05)
   weightsRange = np.empty(len(weights_a),dtype=np.float32)
   for index, iRegion in enumerate(weights_a): 
      weightsRange[index] = iRegion
   
   #

   #Choice of the binning for the efficiency histograms: either we use an adaptive binning or we pickup the binning from a file
   if (Ada):
      
      dfAda = pd.DataFrame(arrayTagAndProb, columns=['L0Calo_ECAL_region_max','L0Calo_ECAL_max_realET','weights'])
      binET = AdaptiveBinning1D(dfAda, 'L0Calo_ECAL_max_realET', 8, VERB)
      print (binET)
      os.system('mkdir -p Bins')
      os.system('mkdir -p Bins/{}'.format(Tag.split("-")[-2]))
      with open('./Bins/{}/BinBoundaries_TAP_E_{}.dat'.format(Tag.split("-")[-2],Tag.split("-")[-2]+"-"+Tag.split("-")[-1]), 'wb') as f:
         for a in binET:
            f.write(str(a)+'\n')
         print "Writing the binning scheme to be adopted in ./Bins/{}/BinBoundaries_TAP_E_{}.dat".format(Tag.split("-")[-2],Tag.split("-")[-2]+"-"+Tag.split("-")[-1])

   else:
      with open('./Bins/{}/BinBoundaries_TAP_E_{}.dat'.format(Tag.split("-")[-2],Tag.split("-")[-2]+"-"+Tag.split("-")[-1]), 'rb') as f:
        fills = f.read().splitlines()
        binET = np.asarray(fills)
      print "Reading the binning scheme to be adopted in ./Bins/{}/BinBoundaries_TAP_E_{}.dat".format(Tag.split("-")[-2],Tag.split("-")[-2]+"-"+Tag.split("-")[-1])

   #binET = [0., 2500., 3500, 4000., 4500., 5000., 6000.,  8000., 20000.] #RKstar
   etRange = np.empty(len(binET),dtype=np.float32)
   for index, iET in enumerate(binET): 
      etRange[index] = iET
   

   w_a = np.arange(-10.,10.)
   ###################Weighted Histogram


   totalHist  = TH2F("totalHist_E"+Tag , "totalHist " , len(regRange)-1, regRange, len(etRange)-1, etRange)
   passedHist = TH2F("passedHist_E"+Tag, "passedHist ", len(regRange)-1, regRange, len(etRange)-1, etRange)
   effHist    = TH2F("effHist_E"+Tag   , "efficiency ", len(regRange)-1, regRange, len(etRange)-1, etRange)
   totalHist.Sumw2(), passedHist.Sumw2(), effHist.Sumw2()
   totalHist_NoW  = TH2F("totalHist_E_NoW"+Tag , "totalHist " , len(regRange)-1, regRange, len(etRange)-1, etRange)
   passedHist_NoW = TH2F("passedHist_E_NoW"+Tag, "passedHist ", len(regRange)-1, regRange, len(etRange)-1, etRange)
   effHist_NoW    = TH2F("effHist_E_NoW"+Tag   , "efficiency ", len(regRange)-1, regRange, len(etRange)-1, etRange)
   totalHist_NoW.Sumw2(), passedHist_NoW.Sumw2(), effHist_NoW.Sumw2()
   total_weights  = TH2F("total_weights_E_NoW"+Tag , "totalHist " ,  len(etRange)-1, etRange,len(weightsRange)-1, weightsRange)
   passed_weights = TH2F("passed_weights_E_NoW"+Tag, "passedHist ",  len(etRange)-1, etRange,len(weightsRange)-1, weightsRange)
   eff_weights = TH2F("eff_weights_E_NoW"+Tag, "passedHist ",  len(etRange)-1, etRange,len(weightsRange)-1, weightsRange)
   total_weights.Sumw2(), passed_weights.Sumw2(), eff_weights.Sumw2()

 # Loop over the data tree and fill the histos                                                                                                          
   for iEntry in xrange(len(arrayTag)):
      totalHist.Fill(arrayTag[iEntry][0], arrayTag[iEntry][1], arrayTag[iEntry][2])
      totalHist_NoW.Fill(arrayTag[iEntry][0], arrayTag[iEntry][1])
      if(arrayTag[iEntry][0]==2):
         total_weights.Fill(arrayTag[iEntry][1], arrayTag[iEntry][2])
   for iEntry in xrange(len(arrayTagAndProb)):
      passedHist.Fill(arrayTagAndProb[iEntry][0], arrayTagAndProb[iEntry][1], arrayTagAndProb[iEntry][2])
      passedHist_NoW.Fill(arrayTagAndProb[iEntry][0], arrayTagAndProb[iEntry][1])
      if(arrayTagAndProb[iEntry][0]==2):
         passed_weights.Fill(arrayTagAndProb[iEntry][1], arrayTagAndProb[iEntry][2])

   effHist.Divide(passedHist, totalHist, 1, 1, "B")
   eff_weights.Divide(passed_weights, total_weights, 1, 1, "B")
   effHist_NoW.Divide(passedHist_NoW, totalHist_NoW, 1, 1, "B")
   effIn  = effHist.ProjectionY("effInE_{}".format(Tag) , 4, 4)
   effMid  = effHist.ProjectionY("effMidE_{}".format(Tag) , 3, 3)
   effOut  = effHist.ProjectionY("effOutE_{}".format(Tag) , 2, 2)


   effIn_NoW  = effHist_NoW.ProjectionY("effInE_NoW_{}".format(Tag) , 4, 4)
   effMid_NoW  = effHist_NoW.ProjectionY("effMidE_NoW_{}".format(Tag) , 3, 3)
   effOut_NoW  = effHist_NoW.ProjectionY("effOutE_NoW_{}".format(Tag) , 2, 2)

   ###################Non Weighted Histogram

   gROOT.SetBatch(kTRUE)
   gROOT.SetStyle("Plain")      
   gROOT.ForceStyle()
   gStyle.SetPalette(1)
   gStyle.SetOptStat(1111111)
   '''
   c = TCanvas('c','')       
   total_weights.Draw("COLZ text")
   c.SaveAs("WeightsDistr-L0E-{}Intotal.pdf".format(Tag))
   c1 = TCanvas('c1','')       
   passed_weights.Draw("COLZ")
   c1.SaveAs("WeightsDistr-L0E-{}Inpassed.pdf".format(Tag))

   c2 = TCanvas('c2','')       
   eff_weights.Draw("COLZ")
   c2.SaveAs("Plots/Checks/WeightsDistr-L0E-{}InEff.pdf".format(Tag))

   c3 = TCanvas('c3','')       
   totalHist.Draw("COLZ text")
   c3.SaveAs("Plots/Checks/totalDistr-L0E-{}.pdf".format(Tag))

   c4 = TCanvas('c4','')       
   passedHist.Draw("COLZ text ")
   c4.SaveAs("Plots/Checks/passedDistr-L0E-{}.pdf".format(Tag))
   '''   





   #The projection is done on the bin number
   '''From the code
   inside = false;

   if (!inside)
   return -1;
   else if (inner)
   return 2;
   else if (middle)
   return 1;
   else if (outer)
   return 0;
   else
   return -999;
   '''
   xBins=[[-0.5, 0.5],[0.5, 1.5],[1.5, 2.5]]
   yBins =[]
   
   effTable_In_NoW = []
   effTable_Mid_NoW = []
   effTable_Out_NoW = []
   ErreffTable_In_NoW = []
   ErreffTable_Mid_NoW = []
   ErreffTable_Out_NoW = []

   for i in range(1,len(etRange)):
      effTable_In_NoW.append(effIn_NoW.GetBinContent(i))
      effTable_Mid_NoW.append(effMid_NoW.GetBinContent(i))
      effTable_Out_NoW.append(effOut_NoW.GetBinContent(i))
      ErreffTable_In_NoW.append(effIn_NoW.GetBinError(i))
      ErreffTable_Mid_NoW.append(effMid_NoW.GetBinError(i))
      ErreffTable_Out_NoW.append(effOut_NoW.GetBinError(i))
      yBins.append([etRange[i-1],etRange[i]])





   dict_tempNoW = {'yBins':yBins, 
                   'L0E-Eff-In':effTable_In_NoW,
                   'L0E-Eff-Mid':effTable_Mid_NoW,
                   'L0E-Eff-Out':effTable_Out_NoW,
                   'L0E-Eff-In-Err':ErreffTable_In_NoW,
                   'L0E-Eff-Mid-Err':ErreffTable_Mid_NoW,
                   'L0E-Eff-Out-Err':ErreffTable_Out_NoW}
   dfEffNoW = pd.DataFrame(dict_tempNoW)

   '''
   print "effTable_In_NoW: ",effTable_In_NoW
   print "effTable_Mid_NoW: ",effTable_Mid_NoW
   print "effTable_Out_NoW: ",effTable_Out_NoW
   print "yBins: ",yBins
   '''



   effTable_In = []
   effTable_Mid = []
   effTable_Out = []
   ErreffTable_In = []
   ErreffTable_Mid = []
   ErreffTable_Out = []

   for i in range(1,len(etRange)):

      effTable_In.append(effIn.GetBinContent(i))
      effTable_Mid.append(effMid.GetBinContent(i))
      effTable_Out.append(effOut.GetBinContent(i))
      ErreffTable_In.append(effIn.GetBinError(i))
      ErreffTable_Mid.append(effMid.GetBinError(i))
      ErreffTable_Out.append(effOut.GetBinError(i))

   '''
   print "effTable_In: ",effTable_In
   print "effTable_Mid: ",effTable_Mid
   print "effTable_Out: ",effTable_Out
   print "yBins: ",yBins
   '''
   


   dict_temp1 = {'Binsy-L0Calo_ECAL_max_realET':yBins,
                 'L0E-Eff-2':effTable_In,
                 'L0E-Eff-1':effTable_Mid,
                 'L0E-Eff-0':effTable_Out,
                 'L0E-Eff-2-Err':ErreffTable_In,
                 'L0E-Eff-1-Err':ErreffTable_Mid,
                 'L0E-Eff-0-Err':ErreffTable_Out}
   
   dfEff1 = pd.DataFrame(dict_temp1)
   
   dict_temp2 = {'Binsx-L0Calo_ECAL_region_max':xBins}
   dfEff2 = pd.DataFrame(dict_temp2)
   
   dfEff = pd.concat([dfEff1,dfEff2],axis=1)


   if VERB:
      print "======Efficiency table - No Weights ======="
      print(dfEffNoW)
   print "======Efficiency table -  Weights ======="
   print(dfEff)
   print "=========================== END of ===================================="
   print "================ Tag&Probe algorithm for L0Electron ==================="
   print "=========================== END of ====================================\n"

   return  effIn, effMid, effOut, dfEff#, passedHist, totalHist
##########################################################################


def TagAndProbe_L0E_single(df1, particle, Tag, Ada, weight, VERB) :    # Open calibration dataset    

   '''
   Author: Michele Atzeni
   Email Address: michele.atzeni@cern.ch
   Date: June 7th, 2017
   Note: script based on T. Humair's B+ -> J/psi K+ algorithm and R.Coutinho

   Description:
   This script computes the trigger-on-signal efficiency for electrons and/or positrons as a function of the trasversal energy deposited in the ECAL using a Tag&Probe approach.
   The electrons are divided depending on the region of the calorimeter hit: inner, medium or outer.

   Concept:
   Ideally to calculate a trigger efficiency you would just do the following:
   
   N_pt = Number of electrons that passed the trigger
   N_t  = Number of total electrons

   e_trigger = N_pt/N_t

   However the quantity N_t is not accessible experimentally.
   We can have a datadriven estimate of the trigger efficiency using a  Tag&Probe approach.

   The idea is very simple: we look at the ratio e_trigger for all those events that passed a specific trigger condition, called TAGGING. If the TAGGING trigger  samples similarly the set of all the electrons and the set of triggered ones the two efficiencies must be the same. 

   e_trigger_TAG = N_pt_tag/ N_t_tag =>e_trigger

   Where:  + N_pt_tag = Number of electrons that are TAGged with B_L0Global_TIS and that are triggered as L#_L0ElectronDecision_TOS(PROBE)
           + N_t_tag = Number of electrons that are TAGged with B_L0Global_TIS 
   

   Scheme:
   1)Get the full dataframe
   2)Choose the charge of the particle: charge =(L1_ID >0) 
                                                (L2_ID >0)

   3)Define the Probe and Tag selections: TAG    => B_L0Global_TIS==1 & L1_ID >0
                                                 => B_L0Global_TIS==1 & L2_ID >0

                                          PROBE  => L1_L0ElectronDecision_TOS==1 & L1_ID >0
                                                 => L2_L0ElectronDecision_TOS==1 & L2_ID >0

   4)Select the dataframe and put together the L1 and L2 
   5)Get the arrays for the TAGged  and the TAG&Probe sample
   6)Calculate efficiency

   '''


   print "================ Tag&Probe algorithm for L0Electron ==================="
   print "=================for {}===============".format(Tag)

   if(not weight):

      branchesTag = ["B_PVandJpsiDTF_B_M","L1_L0ElectronDecision_TOS","B_L0Global_TIS","L1_ID","L2_ID","L2_L0ElectronDecision_TOS","L1_L0Calo_ECAL_region","L2_L0Calo_ECAL_region","L1_L0Calo_ECAL_realET","L2_L0Calo_ECAL_realET"]
      #Reduce the dataframe to speed up
      df = df1[branchesTag]
      df.insert(1, 'wL01',1)
      df.insert(1, 'wL02',1)
      
      if VERB:
         print "========= In the following a uniform weight (1) is used for all the calculations ========="
         print(df[["L1_L0ElectronDecision_TOS","wL01","wL02"]].head())
         print "=========================================================================================="

   else:
      if((weight[0] not in df1.columns) | (weight[1] not in df1.columns)):
      
         branchesTag = ["B_PVandJpsiDTF_B_M","L1_L0ElectronDecision_TOS","B_L0Global_TIS","L1_ID","L2_ID","L2_L0ElectronDecision_TOS","L1_L0Calo_ECAL_region","L2_L0Calo_ECAL_region","L1_L0Calo_ECAL_realET","L2_L0Calo_ECAL_realET"]
         #Reduce the dataframe to speed up
         df = df1[branchesTag]
         df.insert(1, 'wL01',1)
         df.insert(1, 'wL02',1)
         print "=========================================WARNING======================================================="
         print "========= Was requested a weighted histogram, however no column found with the name 'weights' ========="
         print "============== In the following a uniform weight (1) is used for all the calculations ================="
         print(df[["L1_L0ElectronDecision_TOS","wL01","wL02"]].head())
         print "======================================================================================================="
         pdb.set_trace()
      else:

         branchesTag = ["B_PVandJpsiDTF_B_M","L1_L0ElectronDecision_TOS","B_L0Global_TIS","L1_ID","L2_ID","L2_L0ElectronDecision_TOS","L1_L0Calo_ECAL_region","L2_L0Calo_ECAL_region","L1_L0Calo_ECAL_realET","L2_L0Calo_ECAL_realET"]
         branchesTag.append(weight)
         #Reduce the dataframe to speed up
         df = df1[branchesTag]
         df = df.rename(columns={weight[0]:'wL01'})
         df = df.rename(columns={weight[1]:'wL02'})
         if VERB:
            print "============ Column of weights successfully loaded =============="
            print(df[["L1_L0ElectronDecision_TOS","wL01","wL02"]].head())
            print "================================================================="
         
   # Defining the triggered sample (the one we are interested in), the Tag sample and the intersection of the two
   
   print "================================================================================"
   print "====== Number of events in the dataframe used: {}========".format(df.shape[0])
   print "================================================================================"

   #arrayTag, arrayTagAndProb = TAP_L0E_maxET(df, (df.B_L0Global_TIS == 1), VERB)
   #arrayTag, arrayTagAndProb = TAP_L0E2_mixedmaxET(df, (df.B_L0Global_TIS == 1),VERB) 



   selTagDict = {"B_L0Global_TIS" : (df.B_L0Global_TIS == 1),
                 "All" : ((df.B_L0Global_TIS == 1) | (df.B_L0Global_TIS == 0)),
                 "test" : ((df.L1_L0ElectronDecision_TOS == 1) | (df.L2_L0ElectronDecision_TOS == 1)),
}


   print "********** TIS sample used: {} ***********".format(selTag)

   arrayTag, arrayTagAndProb = TAP_L0E_indipE(df, selTagDict[selTag], particle)



   print "================================================================================"
   print "====== Number of events in the TIS Sample: {}========".format(len(arrayTag))
   print "================================================================================"
   print "================================================================================"
   print "====== Number of events in the TOS Sample: {}========".format(len(arrayTagAndProb))
   print "================================================================================"
 
  ######################################################################################            

   print '=================================='
   #print 'The ProbandTag efficiency for the L0Electron trigger in Data is ', float(len(arrayTagAndProb))/float(len(arrayTag))

   #
   locationBinTab = [-1.5, -0.5, 0.5, 1.5, 2.5]
   regRange = np.empty(len(locationBinTab),dtype=np.float32)
   for index, iRegion in enumerate(locationBinTab): 
      regRange[index] = iRegion
   #
   weights_a = np.arange(0,4,0.05)
   weightsRange = np.empty(len(weights_a),dtype=np.float32)
   for index, iRegion in enumerate(weights_a): 
      weightsRange[index] = iRegion
   
   #

   #Choice of the binning for the efficiency histograms: either we use an adaptive binning or we pickup the binning from a file
   if (Ada):
      
      dfAda = pd.DataFrame(arrayTagAndProb, columns=['L0Calo_ECAL_region_max','L0Calo_ECAL_max_realET','weights'])
      binET = AdaptiveBinning1D(dfAda, 'L0Calo_ECAL_max_realET', 8, VERB)
      print (binET)
      os.system('mkdir -p Bins')
      os.system('mkdir -p Bins/{}'.format(Tag.split("-")[-2]))
      with open('./Bins/{}/BinBoundaries_TAP_singleE_{}.dat'.format(Tag.split("-")[-2],Tag.split("-")[-2]+"-"+Tag.split("-")[-1]), 'wb') as f:
         for a in binET:
            f.write(str(a)+'\n')
         print "Writing the binning scheme to be adopted in ./Bins/{}/BinBoundaries_TAP_singleE_{}.dat".format(Tag.split("-")[-2],Tag.split("-")[-2]+"-"+Tag.split("-")[-1])

   else:
      with open('./Bins/{}/BinBoundaries_TAP_singleE_{}.dat'.format(Tag.split("-")[-2],Tag.split("-")[-2]+"-"+Tag.split("-")[-1]), 'rb') as f:
        fills = f.read().splitlines()
        binET = np.asarray(fills)
      print "Reading the binning scheme to be adopted in ./Bins/{}/BinBoundaries_TAP_singleE_{}.dat".format(Tag.split("-")[-2],Tag.split("-")[-2]+"-"+Tag.split("-")[-1])


   #binET = [0., 2500., 3500, 4000., 4500., 5000., 6000.,  8000., 20000.] #RKstar
   etRange = np.empty(len(binET),dtype=np.float32)
   for index, iET in enumerate(binET): 
      etRange[index] = iET
   

   w_a = np.arange(-10.,10.)
   ###################Weighted Histogram


   totalHist  = TH2F("totalHist_E"+Tag , "totalHist " , len(regRange)-1, regRange, len(etRange)-1, etRange)
   passedHist = TH2F("passedHist_E"+Tag, "passedHist ", len(regRange)-1, regRange, len(etRange)-1, etRange)
   effHist    = TH2F("effHist_E"+Tag   , "efficiency ", len(regRange)-1, regRange, len(etRange)-1, etRange)
   totalHist.Sumw2(), passedHist.Sumw2(), effHist.Sumw2()
   totalHist_NoW  = TH2F("totalHist_E_NoW"+Tag , "totalHist " , len(regRange)-1, regRange, len(etRange)-1, etRange)
   passedHist_NoW = TH2F("passedHist_E_NoW"+Tag, "passedHist ", len(regRange)-1, regRange, len(etRange)-1, etRange)
   effHist_NoW    = TH2F("effHist_E_NoW"+Tag   , "efficiency ", len(regRange)-1, regRange, len(etRange)-1, etRange)
   totalHist_NoW.Sumw2(), passedHist_NoW.Sumw2(), effHist_NoW.Sumw2()
   total_weights  = TH2F("total_weights_E_NoW"+Tag , "totalHist " ,  len(etRange)-1, etRange,len(weightsRange)-1, weightsRange)
   passed_weights = TH2F("passed_weights_E_NoW"+Tag, "passedHist ",  len(etRange)-1, etRange,len(weightsRange)-1, weightsRange)
   eff_weights = TH2F("eff_weights_E_NoW"+Tag, "passedHist ",  len(etRange)-1, etRange,len(weightsRange)-1, weightsRange)
   total_weights.Sumw2(), passed_weights.Sumw2(), eff_weights.Sumw2()

 # Loop over the data tree and fill the histos                                                                                                          
   for iEntry in xrange(len(arrayTag)):
      totalHist.Fill(arrayTag[iEntry][0], arrayTag[iEntry][1], arrayTag[iEntry][2])
      totalHist_NoW.Fill(arrayTag[iEntry][0], arrayTag[iEntry][1])
      if(arrayTag[iEntry][0]==2):
         total_weights.Fill(arrayTag[iEntry][1], arrayTag[iEntry][2])
   for iEntry in xrange(len(arrayTagAndProb)):
      passedHist.Fill(arrayTagAndProb[iEntry][0], arrayTagAndProb[iEntry][1], arrayTagAndProb[iEntry][2])
      passedHist_NoW.Fill(arrayTagAndProb[iEntry][0], arrayTagAndProb[iEntry][1])
      if(arrayTagAndProb[iEntry][0]==2):
         passed_weights.Fill(arrayTagAndProb[iEntry][1], arrayTagAndProb[iEntry][2])

   effHist.Divide(passedHist, totalHist, 1, 1, "B")
   eff_weights.Divide(passed_weights, total_weights, 1, 1, "B")
   effHist_NoW.Divide(passedHist_NoW, totalHist_NoW, 1, 1, "B")
   effIn  = effHist.ProjectionY("effInE_{}".format(Tag) , 4, 4)
   effMid  = effHist.ProjectionY("effMidE_{}".format(Tag) , 3, 3)
   effOut  = effHist.ProjectionY("effOutE_{}".format(Tag) , 2, 2)


   effIn_NoW  = effHist_NoW.ProjectionY("effInE_NoW_{}".format(Tag) , 4, 4)
   effMid_NoW  = effHist_NoW.ProjectionY("effMidE_NoW_{}".format(Tag) , 3, 3)
   effOut_NoW  = effHist_NoW.ProjectionY("effOutE_NoW_{}".format(Tag) , 2, 2)

   ###################Non Weighted Histogram

   gROOT.SetBatch(kTRUE)
   gROOT.SetStyle("Plain")      
   gROOT.ForceStyle()
   gStyle.SetPalette(1)
   gStyle.SetOptStat(1111111)
   '''
   c = TCanvas('c','')       
   total_weights.Draw("COLZ text")
   c.SaveAs("WeightsDistr-L0E-{}Intotal.pdf".format(Tag))
   c1 = TCanvas('c1','')       
   passed_weights.Draw("COLZ")
   c1.SaveAs("WeightsDistr-L0E-{}Inpassed.pdf".format(Tag))

   c2 = TCanvas('c2','')       
   eff_weights.Draw("COLZ")
   c2.SaveAs("Plots/Checks/WeightsDistr-L0E-{}InEff.pdf".format(Tag))

   c3 = TCanvas('c3','')       
   totalHist.Draw("COLZ text")
   c3.SaveAs("Plots/Checks/totalDistr-L0E-{}.pdf".format(Tag))

   c4 = TCanvas('c4','')       
   passedHist.Draw("COLZ text ")
   c4.SaveAs("Plots/Checks/passedDistr-L0E-{}.pdf".format(Tag))

   c41 = TCanvas('c41','')       
   effHist.Draw("COLZ text ")
   c41.SaveAs("Plots/Checks/EffDistr-L0E-{}.pdf".format(Tag))
   '''

   #The projection is done on the bin number
   '''From the code
   inside = false;

   if (!inside)
   return -1;
   else if (inner)
   return 2;
   else if (middle)
   return 1;
   else if (outer)
   return 0;
   else
   return -999;
   '''
   xBins=[[-0.5, 0.5],[0.5, 1.5],[1.5, 2.5]]
   yBins =[]
   
   effTable_In_NoW = []
   effTable_Mid_NoW = []
   effTable_Out_NoW = []
   ErreffTable_In_NoW = []
   ErreffTable_Mid_NoW = []
   ErreffTable_Out_NoW = []

   for i in range(1,len(etRange)):
      effTable_In_NoW.append(effIn_NoW.GetBinContent(i))
      effTable_Mid_NoW.append(effMid_NoW.GetBinContent(i))
      effTable_Out_NoW.append(effOut_NoW.GetBinContent(i))
      ErreffTable_In_NoW.append(effIn_NoW.GetBinError(i))
      ErreffTable_Mid_NoW.append(effMid_NoW.GetBinError(i))
      ErreffTable_Out_NoW.append(effOut_NoW.GetBinError(i))
      yBins.append([etRange[i-1],etRange[i]])





   dict_tempNoW = {'yBins':yBins, 
                   'L0E-Eff-In':effTable_In_NoW,
                   'L0E-Eff-Mid':effTable_Mid_NoW,
                   'L0E-Eff-Out':effTable_Out_NoW,
                   'L0E-Eff-In-Err':ErreffTable_In_NoW,
                   'L0E-Eff-Mid-Err':ErreffTable_Mid_NoW,
                   'L0E-Eff-Out-Err':ErreffTable_Out_NoW}
   dfEffNoW = pd.DataFrame(dict_tempNoW)

   '''
   print "effTable_In_NoW: ",effTable_In_NoW
   print "effTable_Mid_NoW: ",effTable_Mid_NoW
   print "effTable_Out_NoW: ",effTable_Out_NoW
   print "yBins: ",yBins
   '''



   effTable_In = []
   effTable_Mid = []
   effTable_Out = []
   ErreffTable_In = []
   ErreffTable_Mid = []
   ErreffTable_Out = []

   for i in range(1,len(etRange)):

      effTable_In.append(effIn.GetBinContent(i))
      effTable_Mid.append(effMid.GetBinContent(i))
      effTable_Out.append(effOut.GetBinContent(i))
      ErreffTable_In.append(effIn.GetBinError(i))
      ErreffTable_Mid.append(effMid.GetBinError(i))
      ErreffTable_Out.append(effOut.GetBinError(i))

   '''
   print "effTable_In: ",effTable_In
   print "effTable_Mid: ",effTable_Mid
   print "effTable_Out: ",effTable_Out
   print "yBins: ",yBins
   '''

   dict_temp1 = {'Binsy-L0Calo_ECAL_max_realET':yBins,
                'L0E-Eff-2':effTable_In,
                'L0E-Eff-1':effTable_Mid,
                'L0E-Eff-0':effTable_Out,
                'L0E-Eff-2-Err':ErreffTable_In,
                'L0E-Eff-1-Err':ErreffTable_Mid,
                'L0E-Eff-0-Err':ErreffTable_Out}

   dfEff1 = pd.DataFrame(dict_temp1)

   dict_temp2 = {'Binsx-L0Calo_ECAL_region_max':xBins}
   dfEff2 = pd.DataFrame(dict_temp2)

   dfEff = pd.concat([dfEff1,dfEff2],axis=1)

   if VERB:
      print "======Efficiency table - No Weights ======="
      print(dfEffNoW)
   print "======Efficiency table -  Weights ======="
   print(dfEff)
   print "=========================== END of ===================================="
   print "================ Tag&Probe algorithm for L0Electron ==================="
   print "=========================== END of ====================================\n"

   return  effIn, effMid, effOut, dfEff#, passedHist, totalHist
 

#########################################################################



def TagAndProbe_L0M(df1, selTag, Model, Tag, Ada, weight, VERB) :    # Open calibration dataset    

   '''
   Author: Michele Atzeni
   Email Address: michele.atzeni@cern.ch
   Date: June 7th, 2017
   Note: script based on T. Humair's B+ -> J/psi K+ algorithm and R.Coutinho

   Description:
   This script computes the trigger-on-signal efficiency for mu+ and/or mu- as a function of the trasversal momentum PT using a Tag&Probe approach.

   Skeleton: 1) Select the events respect to the Tag criteria B_L0Global_TIS==1
             2) Select the events respect to the Tag&Probe criteria (B_L0Global_TIS==1) & (L1_L0MuonDecision_TOS==1 | L2_L0MuonDecision_TOS==1)
             3)Plot the ratio as a function of the PT of the muon (maxPT) or MixedMaxPT.

   '''


   print "================ Tag&Probe algorithm: L0M {}===================".format(Model)
   print "=================for {}===============".format(Tag)

   if(not weight):

      branchesTag = ["B_PVandJpsiDTF_B_M","L1_L0MuonDecision_TOS","B_L0Global_TIS","L1_ID","L2_ID","L2_L0MuonDecision_TOS","L1_PT","L2_PT"]
      #Reduce the dataframe to speed up
      df = df1[branchesTag]
      df.insert(1, 'weights',1)
      if(VERB):
         print "********* In the following a uniform weight (1) is used for all the calculations *********"
         print(df[["L1_L0MuonDecision_TOS","weights"]].head())
         print "******************************************************************************************"

   else:
      if(weight not in df1.columns):
      
         branchesTag = ["B_PVandJpsiDTF_B_M","L1_L0MuonDecision_TOS","B_L0Global_TIS","L1_ID","L2_ID","L2_L0MuonDecision_TOS","L1_PT","L2_PT"]
         #Reduce the dataframe to speed up
         df = df1[branchesTag]
         df.insert(1, 'weights',1)
         
         print "***************************************==WARNING******************************************************="
         print "********* Was requested a weighted histogram, however no column found with the name 'weights' *********"
         print "************== In the following a uniform weight (1) is used for all the calculations ***************=="
         print(df[["L1_L0MuonDecision_TOS","weights"]].head())
         print "******************************************************************************************************="
         exit()
      else:

         branchesTag = ["B_PVandJpsiDTF_B_M","L1_L0MuonDecision_TOS","B_L0Global_TIS","L1_ID","L2_ID","L2_L0MuonDecision_TOS","L1_PT","L2_PT"]
         branchesTag.append(weight)
         #Reduce the dataframe to speed up
         df = df1[branchesTag]
         df = df.rename(columns={weight:'weights'})
         if(VERB):
            print "************ Column of weights successfully loaded ************=="
            print(df[["L1_L0MuonDecision_TOS","weights"]].head())
            print "***************************************************************=="
         
         
   '''
   '''
   # Defining the triggered sample (the one we are interested in), the Tag sample and the intersection of the two
   selTagDict = {"B_L0Global_TIS" : (df.B_L0Global_TIS == 1),
                 "All" : (df.B_L0Global_TIS == 1) | (df.B_L0Global_TIS == 0),
                 "test" : ((df.L1_L0MuonDecision_TOS == 1) | (df.L2_L0MuonDecision_TOS == 1)),
                 #"All_TIS" : ((df.B_L0MuonDecision_TIS == 1) | (df.B_L0DiMuonDecision_TIS == 1) | (df.B_L0ElectronDecision_TIS == 1) | (df.B_L0PhotonDecision_TIS == 1) | (df.B_L0HadronDecision_TIS == 1)),
                 #"All_TIS_exceptE" : ((df.B_L0MuonDecision_TIS == 1) | (df.B_L0DiMuonDecision_TIS == 1) | (df.B_L0PhotonDecision_TIS == 1) | (df.B_L0HadronDecision_TIS == 1)),
              }

   print "TIS sample used: {}".format(selTag)

   if (Model == "maxPT"):
      arrayTag, arrayTagAndProb = TAP_L0M_maxpt(df, selTagDict[selTag])
   elif(Model ==  "mixedmaxPT"):
      arrayTag, arrayTagAndProb = TAP_L0M2_mixedmaxPT(df, selTagDict[selTag], VERB)

   elif(Model ==  "indipM"):
      arrayTag, arrayTagAndProb = TAP_L0M_indipM(df, selTagDict[selTag], VERB)
   else:
      print "WARNING!Wrong Model used!"
      exit()
   #
   if (VERB):
      print '************************== Appending the maximum ET column  ***************************===='
      print "df:\n", df[['PT_max','weights']].head()
      print '*********************************='
      print 'The ProbandTag efficiency for the L0Muon trigger in Data is ', float(len(arrayTagAndProb))/float(len(arrayTag))


   if (Ada):
      dfAda = pd.DataFrame(arrayTagAndProb, columns=['PT_max','weights'])
      binPT = AdaptiveBinning1D(dfAda, 'PT_max', 8, VERB)
      print (binPT)
      os.system('mkdir -p Bins')
      os.system('mkdir -p Bins/{}'.format(Tag.split("-")[-2]))
      with open('./Bins/{}/BinBoundaries_TAP_M_{}.dat'.format(Tag.split("-")[-2],Tag.split("-")[-2]+"-"+Tag.split("-")[-1]), 'wb') as f:
         for a in binPT:
            f.write(str(a)+'\n')
         print "Writing the binning scheme to be adopted in ./Bins/{}/BinBoundaries_TAP_M_{}.dat".format(Tag.split("-")[-2],Tag.split("-")[-2]+"-"+Tag.split("-")[-1])
   else:
      with open('./Bins/{}/BinBoundaries_TAP_M_{}.dat'.format(Tag.split("-")[-2],Tag.split("-")[-2]+"-"+Tag.split("-")[-1]), 'rb') as f:
        fills = f.read().splitlines()
        binPT = np.asarray(fills)
      print "Reading the binning scheme to be adopted in ./Bins/{}/BinBoundaries_TAP_M_{}.dat".format(Tag.split("-")[-2],Tag.split("-")[-2]+"-"+Tag.split("-")[-1])


   #binET = [0., 2500., 3500, 4000., 4500., 5000., 6000.,  8000., 20000.]
   ptRange = np.empty(len(binPT),dtype=np.float32)
   for index, iPT in enumerate(binPT): 
      ptRange[index] = iPT
         
   # Define histograms 
               
   totalHist  = TH1F("totalHist_m_{}".format(Tag) , "totalHist"  , len(ptRange)-1, ptRange)
   passedHist = TH1F("passedHist_m_{}".format(Tag), "passedHist", len(ptRange)-1, ptRange) 
   effHist    = TH1F("effHist_m_{}".format(Tag)   , "L0Muon efficiency (Tag&Probe)", len(ptRange)-1, ptRange)
   totalHist.Sumw2(), passedHist.Sumw2(), effHist.Sumw2() 

   # Loop over the data tree and fill the histos
   for iEntry in xrange(len(arrayTag)): 
      totalHist.Fill(arrayTag[iEntry][0], arrayTag[iEntry][1])

   for iEntry in xrange(len(arrayTagAndProb)): 
      passedHist.Fill(arrayTagAndProb[iEntry][0], arrayTagAndProb[iEntry][1])


   effHist.Divide(passedHist, totalHist, 1, 1, "B")
   #Define the efficiency plots

   # Projections 
   effTable = [] 
   ErreffTable = []
   yBins = []

   for i in range(1,len(ptRange)):

      effTable.append(effHist.GetBinContent(i))
      ErreffTable.append(effHist.GetBinError(i))
      yBins.append([ptRange[i-1],ptRange[i]])

   dict_temp1 = {'Binsy-PT_max':yBins, 
                'L0M-Eff':effTable,
                'L0M-Eff-Err':ErreffTable}
                
   dfEff1 = pd.DataFrame(dict_temp1)

   dict_temp2 = {"Binsx-Year":[[10,13]]}
   dfEff2 = pd.DataFrame(dict_temp2)

   dfEff = pd.concat([dfEff1, dfEff2], axis=1)
   print "======Efficiency table -  Weights ======="
   print(dfEff)

   print "=========================== End of ============================"
   print "================ Tag&Probe algorithm: L0M {}===================".format(Model)
   print "=================for {}===============".format(Tag)

   return  effHist, dfEff
##########################################################################



def TagAndProbe_L0H(df1, selTag, Model, Tag, Ada, weight, VERB):    # Open calibration dataset    

   '''
   Author: Michele Atzeni
   Email Address: michele.atzeni@cern.ch
   Date: June 7th, 2017
   Note: script based on T. Humair's B+ -> J/psi K+ algorithm and R.Coutinho

   Description:
   This script computes the trigger-on-signal efficiency for K_L0HadronDecision_TOS==1 && Pi_L0HadronDecision_TOS==1 as a function of the trasversal momentum PT of the Kstar using a Tag&Probe approach.


   Skeleton: 2) Select the events respect to the Tag criteria B_L0Global_TIS==1
             3) Select the events respect to the Tag&Probe criteria (B_L0Global_TIS==1) & (Kstar_L0Hadron_TOS==1)
             4)Plot the ratio as a function of the PT of the Kstar
   
   '''
   print "=========================== Beginning of ============================"
   print "================ Tag&Probe algorithm for L0Hadron {}===================".format(Model)
   print "=================for {}===============".format(Tag)

   if(not weight):

      branchesTag = ["Pi_L0HadronDecision_TOS","B_L0Global_TIS","K_ID","Pi_ID","K_L0HadronDecision_TOS","K_L0Calo_HCAL_region","Pi_L0Calo_HCAL_region","Pi_L0Calo_HCAL_realET","K_L0Calo_HCAL_realET","Kstar_PT","Kstar_L0HadronDecision_TOS","L1_L0MuonDecision_TOS","L2_L0MuonDecision_TOS","L1_L0ElectronDecision_TOS","L2_L0ElectronDecision_TOS"]
      #Reduce the dataframe to speed up
      df = df1[branchesTag]
      df.insert(1, 'weights',1)
      if VERB:
         print "********* In the following a uniform weight (1) is used for all the calculations *********"
         print(df[["K_ID","weights"]].head())
         print "******************************************************************************************"

   else:
      if(weight not in df1.columns):
      
         branchesTag = ["Pi_L0HadronDecision_TOS","B_L0Global_TIS","K_ID","Pi_ID","K_L0HadronDecision_TOS","K_L0Calo_HCAL_region","Pi_L0Calo_HCAL_region","Pi_L0Calo_HCAL_realET","K_L0Calo_HCAL_realET","Kstar_PT","Kstar_L0HadronDecision_TOS","L1_L0MuonDecision_TOS","L2_L0MuonDecision_TOS","L1_L0ElectronDecision_TOS","L2_L0ElectronDecision_TOS"]
         #Reduce the dataframe to speed up
         df = df1[branchesTag]
         df.insert(1, 'weights',1)
         
         print "***************************************==WARNING******************************************************="
         print "********* Was requested a weighted histogram, however no column found with the name 'weights' *********"
         print "************== In the following a uniform weight (1) is used for all the calculations ***************=="
         print(df[["K_ID","weights"]].head())
         print "******************************************************************************************************="
      
      else:

         branchesTag = ["Pi_L0HadronDecision_TOS","B_L0Global_TIS","K_ID","Pi_ID","K_L0HadronDecision_TOS","K_L0Calo_HCAL_region","Pi_L0Calo_HCAL_region","Pi_L0Calo_HCAL_realET","K_L0Calo_HCAL_realET","Kstar_PT","Kstar_L0HadronDecision_TOS","L1_L0MuonDecision_TOS","L2_L0MuonDecision_TOS","L1_L0ElectronDecision_TOS","L2_L0ElectronDecision_TOS"]
         branchesTag.append(weight)
         #Reduce the dataframe to speed up
         df = df1[branchesTag]
         df = df.rename(columns={weight:'weights'})
         if VERB:
            print "************ Column of weights successfully loaded ************=="
            print(df[["K_ID","weights"]].head())
            print "***************************************************************=="
         
         
   print "******************************************************************************=="
   print "****** Number of events in the dataframe used: {}******==".format(df.shape[0])
   print "******************************************************************************=="
   selTagDict = {"B_L0Global_TIS" : (df.B_L0Global_TIS == 1),
                 "All" : (df.B_L0Global_TIS == 1) | (df.B_L0Global_TIS == 0),
                 "testM" : ((df.Kstar_L0HadronDecision_TOS == 1) & ((df.L1_L0MuonDecision_TOS == 0)&(df.L2_L0MuonDecision_TOS == 0))),
                 "testE" : ((df.Kstar_L0HadronDecision_TOS == 1) & ((df.L1_L0ElectronDecision_TOS == 0)&(df.L2_L0ElectronDecision_TOS == 0))),
                 #"All_TIS" : ((df.B_L0MuonDecision_TIS == 1) | (df.B_L0DiMuonDecision_TIS == 1) | (df.B_L0ElectronDecision_TIS == 1) | (df.B_L0PhotonDecision_TIS == 1) | (df.B_L0HadronDecision_TIS == 1)),
                 #"All_TIS_exceptE" : ((df.B_L0MuonDecision_TIS == 1) | (df.B_L0DiMuonDecision_TIS == 1) | (df.B_L0PhotonDecision_TIS == 1) | (df.B_L0HadronDecision_TIS == 1)),
}

   print "TIS sample used: {}".format(selTag)

   if (Model == "notL0E"):
      nbins = 3
      arrayTag, arrayTagAndProb = TAP_L0H_KstarPT_notL0E(df, selTagDict[selTag]) 

   elif (Model == "notL0M"):
      nbins = 3
      arrayTag, arrayTagAndProb = TAP_L0H_KstarPT_notL0M(df, selTagDict[selTag]) 

   elif (Model == "alsoL0E"):
      arrayTag, arrayTagAndProb = TAP_L0H_KstarPT_alsoL0L(df, selTagDict[selTag]) 
      nbins = 5

   elif (Model == "alsoL0M"):
      arrayTag, arrayTagAndProb = TAP_L0H_KstarPT_alsoL0L(df, selTagDict[selTag]) 
      nbins = 5
   else:
      print "WARNING!No model Found!"
      exit()
                                                                                  


   print "******************************************************************************=="
   print "****** Number of events in the TIS Sample: {}******==".format(len(arrayTag))
   print "******************************************************************************=="
   print "******************************************************************************=="
   print "****** Number of events in the TOS Sample: {}******==".format(len(arrayTagAndProb))
   print "******************************************************************************=="
 
   ########################
   
   locationBinTab = [ -0.5, 0.5, 1.5, 2.5]
   regRange = np.empty(len(locationBinTab),dtype=np.float32)
   for index, iRegion in enumerate(locationBinTab): 
      regRange[index] = iRegion


   if (Ada):
      dfAda = pd.DataFrame(arrayTagAndProb, columns=['K_L0Calo_HCAL_region','Pi_L0Calo_HCAL_region','Kstar_PT',"weights"])
      binET = AdaptiveBinning1D(dfAda, 'Kstar_PT', nbins, VERB)
      print (binET)
      os.system('mkdir -p Bins')
      os.system('mkdir -p Bins/{}'.format(Tag.split("-")[-2]))
      with open('./Bins/{}/BinBoundaries_TAP_L0H_{}.dat'.format(Tag.split("-")[-2],Tag.split("-")[-2]+"-"+Tag.split("-")[-1]), 'wb') as f:
         for a in binET:
            f.write(str(a)+'\n')
         print "Writing the binning scheme to be adopted in ./Bins/{}/BinBoundaries_TAP_L0H_{}.dat".format(Tag.split("-")[-2],Tag.split("-")[-2]+"-"+Tag.split("-")[-1])
   else:
      with open('./Bins/{}/BinBoundaries_TAP_L0H_{}.dat'.format(Tag.split("-")[-2],Tag.split("-")[-2]+"-"+Tag.split("-")[-1]), 'rb') as f:
        fills = f.read().splitlines()
        binET = np.asarray(fills)
      print "Reading the binning scheme to be adopted in ./Bins/{}/BinBoundaries_TAP_L0H_{}.dat".format(Tag.split("-")[-2],Tag.split("-")[-2]+"-"+Tag.split("-")[-1])

   #binET = [0., 4000, 7000, 8000, 9500, 12000, 30000]#RKstar
   etRange = np.empty(len(binET),dtype=np.float32)
   for index, iET in enumerate(binET): 
      etRange[index] = iET
         
   # Define histograms 
   
   totalHist  = TH2F("totalHist_KPi"+Tag , "totalHist " , len(regRange)-1, regRange, len(etRange)-1, etRange)
   passedHist = TH2F("passedHist_KPi"+Tag, "passedHist ", len(regRange)-1, regRange, len(etRange)-1, etRange) 
   effHist    = TH2F("effHist_KPi"+Tag   , "efficiency ", len(regRange)-1, regRange, len(etRange)-1, etRange)
   totalHist.Sumw2(), passedHist.Sumw2(), effHist.Sumw2() 
   totalHist_NoW  = TH2F("totalHist_KPi_NoW"+Tag , "totalHist " , len(regRange)-1, regRange, len(etRange)-1, etRange)
   passedHist_NoW = TH2F("passedHist_KPi_NoW"+Tag, "passedHist ", len(regRange)-1, regRange, len(etRange)-1, etRange) 
   effHist_NoW    = TH2F("effHist_KPi_NoW"+Tag   , "efficiency ", len(regRange)-1, regRange, len(etRange)-1, etRange)
   totalHist_NoW.Sumw2(), passedHist_NoW.Sumw2(), effHist_NoW.Sumw2() 
   
   # Loop over the data tree and fill the histos
   for iEntry in xrange(len(arrayTag)): 
      totalHist_NoW.Fill(arrayTag[iEntry][0]+arrayTag[iEntry][1], arrayTag[iEntry][2])
      totalHist.Fill(arrayTag[iEntry][0]+arrayTag[iEntry][1], arrayTag[iEntry][2], arrayTag[iEntry][3])

   for iEntry in xrange(len(arrayTagAndProb)): 
      passedHist_NoW.Fill(arrayTagAndProb[iEntry][0]+arrayTagAndProb[iEntry][1], arrayTagAndProb[iEntry][2])
      passedHist.Fill(arrayTagAndProb[iEntry][0]+arrayTagAndProb[iEntry][1], arrayTagAndProb[iEntry][2], arrayTagAndProb[iEntry][3])

   effHist_NoW.Divide(passedHist_NoW, totalHist_NoW, 1, 1, "B")
   effHist.Divide(passedHist, totalHist, 1, 1, "B")

   '''
   c3 = TCanvas('c3','')       
   totalHist.Draw("COLZ text")
   c3.SaveAs("Plots/Checks/totalDistr-L0H-{}InEff.pdf".format(Tag))

   c4 = TCanvas('c4','')       
   passedHist.Draw("COLZ text ")
   c4.SaveAs("Plots/Checks/passedDistr-L0H-{}InEff.pdf".format(Tag))
   '''

   # Projections 

   effIn_NoW  = effHist_NoW.ProjectionY("effIn_KPi_NoW{}".format(Tag) , 3, 3)
   effMid_NoW  = effHist_NoW.ProjectionY("effMid_KPi_NoW{}".format(Tag) , 2, 2)
   effOut_NoW = effHist_NoW.ProjectionY("effOut_KPi_NoW{}".format(Tag), 1, 1)
   effIn  = effHist.ProjectionY("effIn_KPi_{}".format(Tag) , 3, 3)
   effMid  = effHist.ProjectionY("effMid_KPi_{}".format(Tag) , 2, 2)
   effOut = effHist.ProjectionY("effOut_KPi_{}".format(Tag), 1, 1)

   yBins =[]
   effTable_In = []
   effTable_Mid = []
   effTable_Out = []
   ErreffTable_In = []
   ErreffTable_Mid = []
   ErreffTable_Out = []

   effTable_In_NoW = []
   effTable_Mid_NoW = []
   effTable_Out_NoW = []
   ErreffTable_In_NoW = []
   ErreffTable_Mid_NoW = []
   ErreffTable_Out_NoW = []
   for i in range(1,len(etRange)):

      effTable_In.append(effIn.GetBinContent(i))
      effTable_Mid.append(effMid.GetBinContent(i))
      effTable_Out.append(effOut.GetBinContent(i))
      ErreffTable_In.append(effIn.GetBinError(i))
      ErreffTable_Mid.append(effMid.GetBinError(i))
      ErreffTable_Out.append(effOut.GetBinError(i))
      yBins.append([etRange[i-1],etRange[i]])
      effTable_In_NoW.append(effIn_NoW.GetBinContent(i))
      effTable_Mid_NoW.append(effMid_NoW.GetBinContent(i))
      effTable_Out_NoW.append(effOut_NoW.GetBinContent(i))
      ErreffTable_In_NoW.append(effIn_NoW.GetBinError(i))
      ErreffTable_Mid_NoW.append(effMid_NoW.GetBinError(i))
      ErreffTable_Out_NoW.append(effOut_NoW.GetBinError(i))

   '''
   print "effTable_In: ",effTable_In
   print "effTable_Mid: ",effTable_Mid
   print "effTable_Out: ",effTable_Out
   print "yBins: ",yBins
   '''
   if(VERB):
      print "====== L0H Efficiency table - No  Weights ======="   
      dict_temp_NoW = {'yBins':yBins, 
                       'L0H-Eff-InIn':effTable_In_NoW,
                       'L0H-Eff-InOut':effTable_Mid_NoW,
                       'L0H-Eff-OutOut':effTable_Out_NoW,
                       'L0H-Eff-InIn-Err':ErreffTable_In_NoW,
                       'L0H-Eff-InOut-Err':ErreffTable_Mid_NoW,
                       'L0H-Eff-OutOut-Err':ErreffTable_Out_NoW}
      dfEff_NoW = pd.DataFrame(dict_temp_NoW)
      print(dfEff_NoW)

   xBins = [[-0.5, 0.5],[0.5, 1.5],[1.5, 2.5]]

   print "====== L0H Efficiency table -  Weights ========="
   dict_temp1 = {'Binsy-Kstar_PT':yBins, 
                'L0H-Eff-2':effTable_In,
                'L0H-Eff-1':effTable_Mid,
                'L0H-Eff-0':effTable_Out,
                'L0H-Eff-2-Err':ErreffTable_In,
                'L0H-Eff-1-Err':ErreffTable_Mid,
                'L0H-Eff-0-Err':ErreffTable_Out}
   dfEff1 = pd.DataFrame(dict_temp1)

   dict_temp2 = { 'Binsx-Kstar_L0Calo_HCAL_region':xBins}
   dfEff2 = pd.DataFrame(dict_temp2)

   dfEff = pd.concat([dfEff1,dfEff2],axis=1)

   print(dfEff)
   print "=========================== END of ===================================="
   print "================ Tag&Probe algorithm for L0Hadron =====================\n\n"
   
   return  effIn, effMid, effOut, dfEff


   ############################

def TagAndProbe_L0TIS(df1, selTag, Model, Tag, Ada, weight, VERB):    # Open calibration dataset    

   '''
   
   '''
   print "=========================== Beginning of ============================"
   print "================ Tag&Probe algorithm for L0TIS {} ===================".format(Model)
   print "=================for {}===============".format(Tag)
   if(not weight):
      branchesTag = ["B_L0Global_TIS","nSPDHits","B_PT","K_L0HadronDecision_TOS","Kstar_L0HadronDecision_TOS","Pi_L0HadronDecision_TOS","L1_L0MuonDecision_TOS","L2_L0MuonDecision_TOS","L1_L0ElectronDecision_TOS","L2_L0ElectronDecision_TOS","B_L0Global_TOS"]
      #Reduce the dataframe to speed up
      df = df1[branchesTag]
      df.insert(1, 'weights',1)
      if VERB:
         print "********* In the following a uniform weight (1) is used for all the calculations *********"
         print(df[["B_L0Global_TIS","weights"]].head())
         print "******************************************************************************************"

   else:
      if(weight not in df1.columns):
         branchesTag = ["B_L0Global_TIS","nSPDHits","B_PT","K_L0HadronDecision_TOS","Kstar_L0HadronDecision_TOS","Pi_L0HadronDecision_TOS","L1_L0MuonDecision_TOS","L2_L0MuonDecision_TOS","L1_L0ElectronDecision_TOS","L2_L0ElectronDecision_TOS","B_L0Global_TOS"]
      #Reduce the dataframe to speed up
         df = df1[branchesTag]
         df.insert(1, 'weights',1)
         
         print "***************************************==WARNING******************************************************="
         print "********* Was requested a weighted histogram, however no column found with the name 'weights' *********"
         print "************== In the following a uniform weight (1) is used for all the calculations ***************=="
         print(df[["B_L0Global_TIS","weights"]].head())
         print "******************************************************************************************************="
         exit()
      else:

         branchesTag = ["B_L0Global_TIS","nSPDHits","B_PT","K_L0HadronDecision_TOS","Pi_L0HadronDecision_TOS","Kstar_L0HadronDecision_TOS","L1_L0MuonDecision_TOS","L2_L0MuonDecision_TOS","L1_L0ElectronDecision_TOS","L2_L0ElectronDecision_TOS","B_L0Global_TOS"]
         branchesTag.append(weight)
         #Reduce the dataframe to speed up
         df = df1[branchesTag]
         df = df.rename(columns={weight:'weights'})
         if VERB:
            print "************ Column of weights successfully loaded ************=="
            print(df[["B_L0Global_TIS","weights"]].head())
            print "***************************************************************=="
         

         
   print "******************************************************************************=="
   print "****** Number of events in the dataframe used: {}******==".format(df.shape[0])
   print "******************************************************************************=="
   selTagDict = {"temp_TIS" : ((df.K_L0HadronDecision_TOS == 1) | (df.Pi_L0HadronDecision_TOS == 1)),
                 "temp_ETIS" : ((df.L1_L0ElectronDecision_TOS == 1) | (df.L2_L0ElectronDecision_TOS == 1)),
                 "All" : ((df.B_L0Global_TIS == 1) | (df.B_L0Global_TIS == 0)),
                 "test" : (df.B_L0Global_TIS == 1),
                 "B_L0Global_TOS" : (df.B_L0Global_TOS == 1),
                 #"All_TIS" : ((df.B_L0MuonDecision_TIS == 1) | (df.B_L0DiMuonDecision_TIS == 1) | (df.B_L0MuonDecision_TIS == 1) | (df.B_L0PhotonDecision_TIS == 1) | (df.B_L0HadronDecision_TIS == 1)),
                 #"All_TIS_exceptE" : ((df.B_L0MuonDecision_TIS == 1) | (df.B_L0DiMuonDecision_TIS == 1) | (df.B_L0PhotonDecision_TIS == 1) | (df.B_L0HadronDecision_TIS == 1)),
}

   print "TIS sample used: {}".format(selTag)

   if (Model == "alsoL0EH"):
      nbins = 11
      arrayTag, arrayTagAndProb = TAP_L0TIS_alsoL0LH(df, selTagDict[selTag]) 
   elif (Model == "alsoL0MH"):
      nbins = 11
      arrayTag, arrayTagAndProb = TAP_L0TIS_alsoL0LH(df, selTagDict[selTag]) 
   elif (Model == "notL0EH"):
      arrayTag, arrayTagAndProb = TAP_L0TIS_notL0EH(df, selTagDict[selTag]) 
      nbins = 5
   elif (Model == "notL0MH"):
      arrayTag, arrayTagAndProb = TAP_L0TIS_notL0MH(df, selTagDict[selTag]) 
      nbins = 5
   else:
      print "WARNING!No model Found!"
      pdb.set_trace()
      exit()



   print "******************************************************************************=="
   print "****** Number of events in the TIS Sample: {}******==".format(len(arrayTag))
   print "******************************************************************************=="
   print "******************************************************************************=="
   print "****** Number of events in the TOS Sample: {}******==".format(len(arrayTagAndProb))
   print "******************************************************************************=="

   ########################
   
   locationBinTab = [ 0., 250., 350., 450., 600.]
   regRange = np.empty(len(locationBinTab),dtype=np.float32)
   for index, iRegion in enumerate(locationBinTab): 
      regRange[index] = iRegion


   if (Ada):
      dfAda = pd.DataFrame(arrayTagAndProb, columns=['nSPDHits','B_PT','weights'])
      binET = AdaptiveBinning1D(dfAda, 'B_PT', nbins, VERB)
      print (binET)
      os.system('mkdir -p Bins')
      os.system('mkdir -p Bins/{}'.format(Tag.split("-")[-2]))
      with open('./Bins/{}/BinBoundaries_TAP_GlobalTIS_{}.dat'.format(Tag.split("-")[-2],Tag.split("-")[-2]+"-"+Tag.split("-")[-1]), 'wb') as f:
         for a in binET:
            f.write(str(a)+'\n')
         print "Writing the binning scheme to be adopted in ./Bins/{}/BinBoundaries_TAP_GlobalTIS_{}.dat".format(Tag.split("-")[-2],Tag.split("-")[-2]+"-"+Tag.split("-")[-1])
   else:
      with open('./Bins/{}/BinBoundaries_TAP_GlobalTIS_{}.dat'.format(Tag.split("-")[-2],Tag.split("-")[-2]+"-"+Tag.split("-")[-1]), 'rb') as f:
        fills = f.read().splitlines()
        binET = np.asarray(fills)
      print "Reading the binning scheme to be adopted in ./Bins/{}/BinBoundaries_TAP_GlobalTIS_{}.dat".format(Tag.split("-")[-2],Tag.split("-")[-2]+"-"+Tag.split("-")[-1])
      
#   binET = [0., 4000, 8000, 9500, 12000, 30000]
   etRange = np.empty(len(binET),dtype=np.float32)
   for index, iET in enumerate(binET): 
      etRange[index] = iET
         
   # Define histograms 
   
   totalHist  = TH2F("totalHist" , "totalHist " , len(regRange)-1, regRange, len(etRange)-1, etRange)
   passedHist = TH2F("passedHist", "passedHist ", len(regRange)-1, regRange, len(etRange)-1, etRange) 
   effHist    = TH2F("effHist"   , "efficiency ", len(regRange)-1, regRange, len(etRange)-1, etRange)
   totalHist.Sumw2(), passedHist.Sumw2(), effHist.Sumw2() 
   
   # Loop over the data tree and fill the histos
   for iEntry in xrange(len(arrayTag)): 
      totalHist.Fill(arrayTag[iEntry][0],arrayTag[iEntry][1], arrayTag[iEntry][2])#, arrayMuonTag[iEntry][2])

   for iEntry in xrange(len(arrayTagAndProb)): 
      passedHist.Fill(arrayTagAndProb[iEntry][0],arrayTagAndProb[iEntry][1], arrayTagAndProb[iEntry][2])#, arrayTagAndProb[iEntry][2])

   effHist.Divide(passedHist, totalHist, 1, 1, "B")



   # Projections 
   eff3= effHist.ProjectionY("eff0_Global_{}".format(Tag), 4, 4)  
   eff2  = effHist.ProjectionY("eff3_Global_{}".format(Tag) , 3, 3)
   eff1  = effHist.ProjectionY("eff2_Global_{}".format(Tag) , 2, 2)
   eff0= effHist.ProjectionY("eff1_Global_{}".format(Tag), 1, 1)


   xBins = [[0.,200.],[200.,350.],[350.,450.],[450.,600.]]
   yBins =[]
   effTable_0 = []
   effTable_1 = []
   effTable_2 = []
   effTable_3 = []
   ErreffTable_0 = []
   ErreffTable_1 = []
   ErreffTable_2 = []
   ErreffTable_3 = []

   for i in range(1,len(etRange)):

      effTable_0.append(eff0.GetBinContent(i))
      effTable_1.append(eff1.GetBinContent(i))
      effTable_2.append(eff2.GetBinContent(i))
      effTable_3.append(eff3.GetBinContent(i))
      ErreffTable_0.append(eff0.GetBinError(i))
      ErreffTable_1.append(eff1.GetBinError(i))
      ErreffTable_2.append(eff2.GetBinError(i))
      ErreffTable_3.append(eff3.GetBinError(i))
      yBins.append([etRange[i-1],etRange[i]])

   print "====== L0TIS Efficiency table -  Weights ========="
   dict_temp1 = {'Binsy-B_PT':yBins, 
                'L0TIS-Eff-0':effTable_0,
                'L0TIS-Eff-1':effTable_1,
                'L0TIS-Eff-2':effTable_2,
                'L0TIS-Eff-3':effTable_3,
                'L0TIS-Eff-0-Err':ErreffTable_0,
                'L0TIS-Eff-1-Err':ErreffTable_1,
                'L0TIS-Eff-2-Err':ErreffTable_2,
                'L0TIS-Eff-3-Err':ErreffTable_3}

   dfEff1 = pd.DataFrame(dict_temp1)
   dict_temp2 = {'Binsx-nSPDHits':xBins}
   dfEff2 = pd.DataFrame(dict_temp2)

   dfEff = pd.concat([dfEff1,dfEff2],axis=1)
   print(dfEff)
   print "=========================== END of ===================================="
   print "================ Tag&Probe algorithm for L0TIS ========================\n\n"
   

   return  eff0, eff1, eff2, eff3, dfEff

   #######################################################

def TagAndProbe_HLT(df1, year, Model, Tag, Ada, weight, VERB) :    # Open calibration dataset    

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


   print "================ Tag&Probe algorithm for HLT-{}===================".format(Model)
   print "=================for {}===============".format(Tag)

   if (Model == 'E'):
      if ((year == '11') or (year == '12') or (year == 'RunI')) :
         branchesTag = ['L1_PT','L2_PT','K_PT','Pi_PT',"B_PT","B_Hlt1TrackAllL0Decision_TOS","B_Hlt2Topo2BodyBBDTDecision_TOS","B_Hlt2Topo3BodyBBDTDecision_TOS","B_Hlt2Topo4BodyBBDTDecision_TOS","B_Hlt2TopoE2BodyBBDTDecision_TOS","B_Hlt2TopoE3BodyBBDTDecision_TOS","B_Hlt2TopoE4BodyBBDTDecision_TOS","B_Hlt1Phys_TIS","B_Hlt2Phys_TIS"]
      elif ((year == '15') or (year == '16')or (year == 'RunII')) :
         branchesTag = ['L1_PT','L2_PT','K_PT','Pi_PT',"B_PT","B_Hlt1TrackMVADecision_TOS","B_Hlt1TwoTrackMVADecision_TOS","B_Hlt1ElectronTrackDecision_TOS","B_Hlt1TrackMVALooseDecision_TOS","B_Hlt1TwoTrackMVALooseDecision_TOS","B_Hlt2Topo2BodyDecision_TOS","B_Hlt2Topo3BodyDecision_TOS","B_Hlt2Topo4BodyDecision_TOS","B_Hlt2TopoE2BodyDecision_TOS","B_Hlt2TopoE3BodyDecision_TOS","B_Hlt2TopoE4BodyDecision_TOS","B_Hlt2TopoEE2BodyDecision_TOS","B_Hlt2TopoEE3BodyDecision_TOS","B_Hlt2TopoEE4BodyDecision_TOS","B_Hlt1Phys_TIS","B_Hlt2Phys_TIS"]
   elif(Model == 'M'):

      branchesTag = ['L1_PT','L2_PT','K_PT','Pi_PT',"B_PT","B_Hlt1TrackAllL0Decision_TOS","B_Hlt1TrackMuonDecision_TOS","B_Hlt2Topo2BodyBBDTDecision_TOS","B_Hlt2Topo3BodyBBDTDecision_TOS","B_Hlt2Topo4BodyBBDTDecision_TOS","B_Hlt2TopoMu2BodyBBDTDecision_TOS","B_Hlt2TopoMu3BodyBBDTDecision_TOS","B_Hlt2TopoMu4BodyBBDTDecision_TOS","B_Hlt2DiMuonDetachedDecision_TOS","B_Hlt1Phys_TIS","B_Hlt2Phys_TIS"]

   else:
      print "WARNING!Wrong Model used!"
      exit()
   #Reduce the dataframe to speed up
   df = df1[branchesTag]

   if(not weight):

      df.insert(1, 'weights',1)
      if(VERB):
         print "********* In the following a uniform weight (1) is used for all the calculations *********"
         print(df[["B_PT","weights"]].head())
         print "******************************************************************************************"

   else:
      if(weight not in df.columns):
         df.insert(1, 'weights',1)
         
         print "***************************************==WARNING******************************************************="
         print "********* Was requested a weighted histogram, however no column found with the name 'weights' *********"
         print "************== In the following a uniform weight (1) is used for all the calculations ***************=="
         print(df[["B_PT","weights"]].head())
         print "******************************************************************************************************="
         exit()
      else:

         df = df.rename(columns={weight:'weights'})
         if(VERB):
            print "============ Column of weights successfully loaded =============="
            print(df[["B_PT","weights"]].head())
            print "================================================================="

         
   '''
   '''
         
   print "******************************************************************************=="
   print "****** Number of events in the dataframe used: {}******==".format(df.shape[0])
   print "******************************************************************************=="


   print "TIS sample used: B_HLt1Phys_TIS && B_HLt2Phys_TIS"

   if (Model == "E"):
      arrayTag, arrayTagAndProb = TAP_HLT_E(df, year, VERB)
   elif(Model ==  "M"):
      arrayTag, arrayTagAndProb = TAP_HLT_M(df, year, VERB)
   else:
      print "WARNING!Wrong Model used!"
      exit()


   print "******************************************************************************=="
   print "****** Number of events in the TIS Sample: {}******==".format(len(arrayTag))
   print "******************************************************************************=="
   print "******************************************************************************=="
   print "****** Number of events in the TOS Sample: {}******==".format(len(arrayTagAndProb))
   print "******************************************************************************=="

   # Defining the triggered sample (the one we are interested in), the Tag sample and the intersection of the two

   #

   


   if (VERB):
      print '========================== Appending the minimum ET column  ==============================='
      print "df:\n", df[['PT_min','weights']].head()
      print '=================================='
      print 'The ProbandTag efficiency for the HLT trigger in Data is ', float(len(arrayTagAndProb))/float(len(arrayTag))

      #TO DO: modify the variable used to bin
   if (Ada):
      dfAda = pd.DataFrame(arrayTagAndProb, columns=['PT_min','weights'])

      binPT = AdaptiveBinning1D(dfAda, 'PT_min', 6, VERB)
      print (binPT)
      os.system('mkdir -p Bins')
      os.system('mkdir -p Bins/{}'.format(Tag.split("-")[-2]))
      with open('./Bins/{}/BinBoundaries_TAP_HLT_{}.dat'.format(Tag.split("-")[-2],Tag.split("-")[-2]+"-"+Tag.split("-")[-1]), 'wb') as f:#Change the name!!!
         for a in binPT:
            f.write(str(a)+'\n')
         print "Writing the binning scheme to be adopted in ./Bins/{}/BinBoundaries_TAP_HLT_{}.dat".format(Tag.split("-")[-2],Tag.split("-")[-2]+"-"+Tag.split("-")[-1])
   else:
      with open('./Bins/{}/BinBoundaries_TAP_HLT_{}.dat'.format(Tag.split("-")[-2],Tag.split("-")[-2]+"-"+Tag.split("-")[-1]), 'rb') as f:
        fills = f.read().splitlines()
        binPT = np.asarray(fills)
      print "Writing the binning scheme to be adopted in ./Bins/{}/BinBoundaries_TAP_HLT_{}.dat".format(Tag.split("-")[-2],Tag.split("-")[-2]+"-"+Tag.split("-")[-1])

   #binET = [0., 2500., 3500, 4000., 4500., 5000., 6000.,  8000., 20000.]
   ptRange = np.empty(len(binPT),dtype=np.float32)
   for index, iPT in enumerate(binPT): 
      ptRange[index] = iPT
         
   # Define histograms 
               
   totalHist  = TH1F("totalHist_m_{}".format(Tag) , "totalHist"  , len(ptRange)-1, ptRange)
   passedHist = TH1F("passedHist_m_{}".format(Tag), "passedHist", len(ptRange)-1, ptRange) 
   effHist    = TH1F("effHist_m_{}".format(Tag)   , "HLT efficiency (Tag&Probe) for {}".format(Model), len(ptRange)-1, ptRange)
   totalHist.Sumw2(), passedHist.Sumw2(), effHist.Sumw2() 

   # Loop over the data tree and fill the histos
   for iEntry in xrange(len(arrayTag)): 
      totalHist.Fill(arrayTag[iEntry][0], arrayTag[iEntry][1])

   for iEntry in xrange(len(arrayTagAndProb)): 
      passedHist.Fill(arrayTagAndProb[iEntry][0], arrayTagAndProb[iEntry][1])


   effHist.Divide(passedHist, totalHist, 1, 1, "B")
   #Define the efficiency plots

   # Projections 
   effTable = []
   ErreffTable = []
   ptBins = []

   for i in range(1,len(ptRange)):

      effTable.append(effHist.GetBinContent(i))
      ErreffTable.append(effHist.GetBinError(i))
      ptBins.append([ptRange[i-1],ptRange[i]])

   dict_temp1 = {'Binsy-PT_min':ptBins, 
                 'HLT-Eff':effTable,
                 'HLT-Eff-Err':ErreffTable}
                
   dfEff1 = pd.DataFrame(dict_temp1)

   dict_temp2 = {"Binsx-Year":[[10,13]]}
   dfEff2 = pd.DataFrame(dict_temp2)

   dfEff = pd.concat([dfEff1, dfEff2], axis=1)

   print "======== Efficiency table -  Weights ========="
   print(dfEff)

   print "===== End of HLT TagAndProbe Efficiency ======\n\n"

   return  effHist, dfEff

##################

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
   
