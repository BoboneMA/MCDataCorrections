from ROOT import gROOT, gStyle, TFile, TLorentzVector, TVector3, TRotation, TLorentzRotation, TMath, TH1D, TCanvas, TH2D, TObject, TF1, TH1F, TLegend, kTRUE, gDirectory
import argparse
import pdb



def ComparisonPlots(H, title, yrange, xtitle, ytitle, legend, plot_name):

    '''
    Script that compacts the plotting procedure for two histograms
    '''
    gROOT.SetBatch(kTRUE)
    gROOT.SetStyle("Plain")
    gROOT.ForceStyle()

    gStyle.SetOptStat(1111111)
    #gStyle.SetOptStat(0)
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
    print "Plotting  ==> {}\n".format(plot_name)

    del c
    del H
    del leg0


def GetKeyNames( self, dir = "" ):
    self.GetListOfKeys().Print()
    return [key.GetName() for key in self.GetListOfKeys()]



def Divide_TH1Dlist(hist_Data, hist_MC):
    hist=[]
    for i,histMC in enumerate(hist_MC):
        print "{},{}".format(i,histMC)

        hist_temp = hist_Data[i].Clone()
        hist_temp.Divide(hist_MC[i])
        hist_temp.SetName(hist_Data[i].GetName().replace("_Data",""))
        hist.append(hist_temp)
    return hist

if __name__ == "__main__" :

    parser = argparse.ArgumentParser(description = 'Configuration of the parameters for the SplitAndMerge')

    parser.add_argument("-r", "--run" , dest="r"  , required=False, help="Choose the year for the RI-RII", choices=['RunI','RunII'], default = 'RunI')
    parser.add_argument("-v", "--versionq2" , dest="version"  , required=False, help="", choices= ['All','q2','q2_PVandBDTF'], default = 'All')
    parser.add_argument("--TM", action="store_true", help="TruthMatchedMC")
    parser.add_argument("-low", "--lowq2" , dest="low_val"  , required=False, help="", choices= ['6','7'], default = '6')

    parser.add_argument("--plot", action="store_true", help="Do a test plot")
    parser.add_argument("--VERB", action="store_true", help="VERBOSE")



    args = parser.parse_args()

    # Parameters and configuration                                                                                                                                                 

    r= args.r
    version= args.version
    low_val= args.low_val
    TM = args.TM
    plot = args.plot
    VERB = args.VERB


    if(TM):
        Tag_sample = 'q2{}_lw{}_TM'.format(version, low_val)
    else:
        Tag_sample = 'q2{}_lw{}'.format(version, low_val)
        
    yearlist = ['11','12']
    for iyear in yearlist:

                                                                                                                                                                                 
        #Open the efficiency tables for the correction                                                                                                                                 
        print "Opening the file./EffHisto_Calib_L0_{}_{}_{}-{}.root".format("B02Kst0Jpsi2ee", iyear, "MC", Tag_sample)
        file_ee_root_MC = TFile("./EffHisto_Calib_L0_{}_{}_{}-{}.root".format("B02Kst0Jpsi2ee", iyear, "MC", Tag_sample))
        print "Opening the file./EffHisto_Calib_L0_{}_{}_{}-q2{}_{}.root".format("B02Kst0ee", iyear, "Data", version, low_val)
        file_ee_root_Data = TFile("./EffHisto_Calib_L0_{}_{}_{}-q2{}_lw{}.root".format("B02Kst0ee", iyear, "Data", version, low_val))
        
        file_mm_root_MC = TFile("./EffHisto_Calib_L0_{}_{}_{}-{}.root".format("B02Kst0Jpsi2mm", iyear, "MC", Tag_sample))   
        file_mm_root_Data = TFile("./EffHisto_Calib_L0_{}_{}_{}-q2{}_lw{}.root".format("B02Kst0mm",iyear, "Data", version, low_val))
        TFile.GetKeyNames = GetKeyNames


        hist_MCee = [file_ee_root_MC.Get(key) for key in file_ee_root_MC.GetKeyNames()]
        hist_Dataee =[file_ee_root_Data.Get(key) for key in  file_ee_root_Data.GetKeyNames()]


        hist_MCmm = [file_mm_root_MC.Get(key) for key in file_mm_root_MC.GetKeyNames()]
        hist_Datamm =[file_mm_root_Data.Get(key) for key in  file_mm_root_Data.GetKeyNames()]


        hist_ee = Divide_TH1Dlist(hist_Dataee, hist_MCee)        
        hist_mm = Divide_TH1Dlist(hist_Datamm, hist_MCmm)        

        os.system("mkdir -p Plots")

        print(hist_ee[0])
        ComparisonPlots([hist_ee[0]],
                        " Data/MC ratio for L0E - Inner region " ,
                        [0.,2.0],
                        "e real^{L0CaloTool}_{max(E_{T})}[MeV]",
                        "#epsilon^{Data}/#epsilon^{MC}",
                        ["Ratio Data/MC-{}".format(iyear)],
                        "./Plots/TagandProbe_ee_DataMCratio_In_e_20{}_maxET.pdf".format(iyear))
        
        print(hist_ee[1])
        ComparisonPlots([hist_ee[1]],
                        " Data/MC ratio for L0E - Middle region " ,
                        [0.,2.0],
                        "e real^{L0CaloTool}_{max(E_{T})}[MeV]",
                        "#epsilon^{Data}/#epsilon^{MC}",
                        ["Ratio Data/MC-{}".format(iyear)],
                        "./Plots/TagandProbe_ee_DataMCratio_Mid_e_20{}_maxET.pdf".format(iyear))
        
        print(hist_ee[2])
        ComparisonPlots([hist_ee[2]],
                        " Data/MC ratio for L0E - Outer region " ,
                        [0.,2.0],
                        "e real^{L0CaloTool}_{max(E_{T})}[MeV]",
                        "#epsilon^{Data}/#epsilon^{MC}",
                        ["Ratio Data/MC-{}".format(iyear)],
                        "./Plots/TagandProbe_ee_DataMCratio_Out_e_20{}_maxET.pdf".format(iyear))
        
        
        print(hist_ee[3])
        ComparisonPlots([hist_ee[3]],
                        " Data/MC ratio for L0E - Inner region " ,
                        [0.,2.0],
                        "e real^{L0CaloTool}_{mixedmax(E_{T})}[MeV]",
                        "#epsilon^{Data}/#epsilon^{MC}",
                        ["Ratio Data/MC-{}".format(iyear)],
                        "./Plots/TagandProbe_ee_DataMCratio_In_e_20{}_mixedmaxET.pdf".format(iyear))
        
        print(hist_ee[4])
        ComparisonPlots([hist_ee[4]],
                        
                        " Data/MC ratio for L0E - Middle region " ,
                        [0.,2.0],
                        "e real^{L0CaloTool}_{mixedmax(E_{T})}[MeV]",
                        "#epsilon^{Data}/#epsilon^{MC}",
                        ["Ratio Data/MC-{}".format(iyear)],
                        "./Plots/TagandProbe_ee_DataMCratio_Mid_e_20{}_mixedmaxET.pdf".format(iyear))
        
        print(hist_ee[5])
        ComparisonPlots([hist_ee[5]],
                        
                        " Data/MC ratio for L0E - Outer region " ,
                        [0.,2.0],
                        "e real^{L0CaloTool}_{mixedmax(E_{T})}[MeV]",
                        "#epsilon^{Data}/#epsilon^{MC}",
                        ["Ratio Data/MC-{}".format(iyear)],
                        "./Plots/TagandProbe_ee_DataMCratio_Out_e_20{}_mixedmaxET.pdf".format(iyear))
        
        print(hist_mm[0])
        ComparisonPlots([hist_mm[0]],
                        
                        " Data/MC ratio for L0M - Inner region " ,
                        [0.9,1.1],#                        [0.,2.0],
                        "e real^{L0CaloTool}_{max(E_{T})}[MeV]",
                        "#epsilon^{Data}/#epsilon^{MC}",
                        ["Ratio Data/MC-{}".format(iyear)],
                        "./Plots/TagandProbe_mm_DataMCratio_M_20{}_maxPT.pdf".format(iyear))
        
                
        
        #L0M
        print(hist_mm[1])
        ComparisonPlots([hist_mm[1]],
                        
                        " Data/MC ratio for L0M - Inner region " ,
                        [0.7,1.3],
                        "e real^{L0CaloTool}_{mixedmax(E_{T})}[MeV]",
                        "#epsilon^{Data}/#epsilon^{MC}",
                        ["Ratio Data/MC-{}".format(iyear)],
                        "./Plots/TagandProbe_mm_DataMCratio_M_20{}_mixedmaxPT.pdf".format(iyear))
                        

                        
        print(hist_mm[2],hist_ee[6])
        ComparisonPlots([hist_mm[2],hist_ee[6]],
                        
                        " Data/MC ratio for L0H - Inner-Inner region " ,
                        [0.,2.0],
                        "K^{*0} p_{T}[MeV/c]",
                        "#epsilon^{Data}/#epsilon^{MC}",
                        ["Ratio Data/MC-mm-{}".format(iyear),"Ratio Data/MC-ee-{}".format(iyear)],
                        "./Plots/TagandProbe_DataMCratio_In_KPi_20{}_notL0E.pdf".format(iyear))
        
        print(hist_mm[3],hist_ee[7])
        ComparisonPlots([hist_mm[3],hist_ee[7]],
                        
                        " Data/MC ratio for L0H - Inner-Outer region " ,
                        [0.,2.0],
                        "K^{*0} p_{T}[MeV/c]",
                        "#epsilon^{Data}/#epsilon^{MC}",
                        ["Ratio Data/MC-mm-{}".format(iyear),"Ratio Data/MC-ee-{}".format(iyear)],
                        "./Plots/TagandProbe_ee_DataMCratio_Mid_KPi_20{}_notL0E.pdf".format(iyear))
        
        print(hist_mm[4],hist_ee[8])
        ComparisonPlots([hist_mm[4],hist_ee[8]],
                        
                        " Data/MC ratio for L0H - Outer-Outer region " ,
                        [0.,2.0],
                        "K^{*0} p_{T}[MeV/c]",
                        "#epsilon^{Data}/#epsilon^{MC}",
                        ["Ratio Data/MC-mm-{}".format(iyear),"Ratio Data/MC-ee-{}".format(iyear)],
                        "./Plots/TagandProbe_ee_DataMCratio_Out_KPi_20{}_notL0E.pdf".format(iyear))
        
        print(hist_mm[5],hist_ee[9])
        ComparisonPlots([hist_mm[5],hist_ee[9]],
                        
                        " Data/MC ratio for L0H - Inner-Inner region " ,
                        [0.,2.0],
                        "K^{*0} p_{T}[MeV/c]",
                        "#epsilon^{Data}/#epsilon^{MC}",
                        ["Ratio Data/MC-mm-{}".format(iyear),"Ratio Data/MC-ee-{}".format(iyear)],
                        "./Plots/TagandProbe_ee_DataMCratio_In_KPi_20{}_alsoL0E.pdf".format(iyear))
        
        print(hist_mm[6],hist_ee[10])
        ComparisonPlots([hist_mm[6],hist_ee[10]],
                        
                        " Data/MC ratio for L0H - Inner-Outer region " ,
                        [0.,2.0],
                        "K^{*0} p_{T}[MeV/c]",
                        "#epsilon^{Data}/#epsilon^{MC}",
                        ["Ratio Data/MC-mm-{}".format(iyear),"Ratio Data/MC-ee-{}".format(iyear)],
                        "./Plots/TagandProbe_ee_DataMCratio_Mid_KPi_20{}_alsoL0E.pdf".format(iyear))
        
        print(hist_mm[7],hist_ee[11])
        ComparisonPlots([hist_mm[7],hist_ee[11]],
                        
                        " Data/MC ratio for L0H - Outer-Outer region " ,
                        [0.,2.0],
                        "K^{*0} p_{T}[MeV/c]",
                        "#epsilon^{Data}/#epsilon^{MC}",
                        ["Ratio Data/MC-mm-{}".format(iyear),"Ratio Data/MC-ee-{}".format(iyear)],
                        "./Plots/TagandProbe_ee_DataMCratio_Out_KPi_20{}_alsoL0E.pdf".format(iyear))
        
        
        print(hist_mm[8],hist_ee[12])
        ComparisonPlots([hist_mm[8],hist_ee[12]],
                        "Data/MC ratio for L0TIS - (0,250) region" , 
                        [0.,2.0],
                        "B^{0} p_{T}[MeV/c]",
                        "#epsilon^{Data}/#epsilon^{MC}", 
                        ["Ratio Data/MC-mm-{}".format(iyear),"Ratio Data/MC-ee-{}".format(iyear)],
                        "./Plots/TagandProbe_ee_DataMCratio_TIS_1_notL0EH_20{}.pdf".format(iyear))

        print(hist_mm[9],hist_ee[13])
        ComparisonPlots([hist_mm[9],hist_ee[13]],
                        "Data/MC ratio for L0TIS - (250,350) region" , 
                        [0.,2.0],
                        "B^{0} p_{T}[MeV/c]",
                        "#epsilon^{Data}/#epsilon^{MC}", 
                        ["Ratio Data/MC-mm-{}".format(iyear),"Ratio Data/MC-ee-{}".format(iyear)],
                        "./Plots/TagandProbe_ee_DataMCratio_TIS_2_notL0EH_20{}.pdf".format(iyear))

        print(hist_mm[10],hist_ee[14])
        ComparisonPlots([hist_mm[10],hist_ee[14]],
                        "Data/MC ratio for L0TIS - (350,450) region" , 
                        [0.,2.0],
                        "B^{0} p_{T}[MeV/c]",
                        "#epsilon^{Data}/#epsilon^{MC}", 
                        ["Ratio Data/MC-mm-{}".format(iyear),"Ratio Data/MC-ee-{}".format(iyear)],
                        "./Plots/TagandProbe_ee_DataMCratio_TIS_3_notL0EH_20{}.pdf".format(iyear))

        print(hist_mm[11],hist_ee[15])
        ComparisonPlots([hist_mm[11],hist_ee[15]],
                        "Data/MC ratio for L0TIS - (450,600) region" , 
                        [0.,2.0],
                        "B^{0} p_{T}[MeV/c]",
                        "#epsilon^{Data}/#epsilon^{MC}", 
                        ["Ratio Data/MC-mm-{}".format(iyear),"Ratio Data/MC-ee-{}".format(iyear)],
                        "./Plots/TagandProbe_ee_DataMCratio_TIS_0_notL0EH_20{}.pdf".format(iyear))

        print(hist_mm[12],hist_ee[16])
        ComparisonPlots([hist_mm[12],hist_ee[16]],
                        "Data/MC ratio for L0TIS - (0,250) region" , 
                        [0.,2.0],
                        "B^{0} p_{T}[MeV/c]",
                        "#epsilon^{Data}/#epsilon^{MC}", 
                        ["Ratio Data/MC-mm-{}".format(iyear),"Ratio Data/MC-ee-{}".format(iyear)],
                        "./Plots/TagandProbe_ee_DataMCratio_TIS_1_alsoL0EH_20{}.pdf".format(iyear))
        print(hist_mm[13],hist_ee[17])
        ComparisonPlots([hist_mm[13],hist_ee[17]],
                        "Data/MC ratio for L0TIS - (250,350) region" , 
                        [0.,2.0],
                        "B^{0} p_{T}[MeV/c]",
                        "#epsilon^{Data}/#epsilon^{MC}", 
                        ["Ratio Data/MC-mm-{}".format(iyear),"Ratio Data/MC-ee-{}".format(iyear)],
                        "./Plots/TagandProbe_ee_DataMCratio_TIS_2_alsoL0EH_20{}.pdf".format(iyear))
        print(hist_mm[14],hist_ee[18])
        ComparisonPlots([hist_mm[14],hist_ee[18]],
                        "Data/MC ratio for  L0TIS - (350,450) region" , 
                        [0.,2.0],
                        "B^{0} p_{T}[MeV/c]",
                        "#epsilon^{Data}/#epsilon^{MC}", 
                        ["Ratio Data/MC-mm-{}".format(iyear),"Ratio Data/MC-ee-{}".format(iyear)],
                        "./Plots/TagandProbe_ee_DataMCratio_TIS_3_alsoL0EH_20{}.pdf".format(iyear))
        print(hist_mm[15],hist_ee[19])
        ComparisonPlots([hist_mm[15],hist_ee[19]],
                        "Data/MC ratio for  L0TIS - (450,600) region" , 
                        [0.,2.0],
                        "B^{0} p_{T}[MeV/c]",
                        "#epsilon^{Data}/#epsilon^{MC}", 
                        ["Ratio Data/MC-mm-{}".format(iyear),"Ratio Data/MC-ee-{}".format(iyear)],
                        "./Plots/TagandProbe_ee_DataMCratio_TIS_0_alsoL0EH_20{}.pdf".format(iyear))
        

