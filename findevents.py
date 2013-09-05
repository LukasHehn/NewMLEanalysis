#!/usr/bin/env python

######################################################
#
# Filter and display events from a KData file
# Lukas Hehn, 2013
#
######################################################

import functions
import ROOT

from parameters import BASELINE_CUTS_ERIC
from ROOT import TGraph, gROOT, TCanvas, KDataReader


# Definition of parameters used for skimming event set
DETECTOR_NAME = 'ID3'
KDataFile = 'Data/Run12_ID3_bckg_with_subrecords.root'
OutFileName = False  # 'Data/ID3_eventlist_lowE_corrected.txt'
E_ION_MAX = 10.
E_REC_MAX = 20.
IonFactor = 1.0/1.029
HeatFactor = 1.0/1.017
VOLTAGE = 6.4
E_THRESH = 3.874
FWHM_HEAT = 0.82
FWHM_ION = 0.72
FWHM_REC = functions.fwhm_rec_from_heat(FWHM_HEAT, VOLTAGE, 10)
GAMMA_CUT = 1.96  # gamma cut from ER centroid in standard deviations
WIMP_MASS = 10  # if wimp mass set, show also signal density of wimp signal


CUTS = BASELINE_CUTS_ERIC[DETECTOR_NAME]
TimeList, RecList, IonList = [],  [],  []  # for normal event selection
TimeListCuts, RecListCuts, IonListCuts = [],  [],  []  # for additional cut event selection
eventcounter, cutcounter = 0, 0


infile = KDataReader(KDataFile)

# Get first event and total event number
event = infile.GetEvent()
entries = infile.GetEntries()
print entries, "total number of events in file", KDataFile, "\n"

#stamp = event.GetStamp()
for entry in range(entries):
    #print "entry", entry
    infile.GetEntry(entry)
    bolos = event.GetNumBolos()

    for bolonum in range(bolos):
    #if bolos == 1:
       #bolonum = 0

        #bolo records
        bolo = event.GetBolo(bolonum)

        if IonFactor and HeatFactor:
            EnergyIon = IonFactor*bolo.GetEnergyIonFiducial()
            EnergyHeat = HeatFactor*bolo.GetEnergyHeat(1)
        else:
            EnergyIon = bolo.GetEnergyIonFiducial()
            EnergyHeat = bolo.GetEnergyHeat(1)

        #samba records
        samba = bolo.GetSambaRecord()
    
        EventNumber = samba.GetSambaEventNumber()

        UnixTime = samba.GetNtpDateSec()
        RunName = samba.GetRunName()

        if bolo.GetDetectorName()==DETECTOR_NAME\
            and CUTS['Heat']['Min']<bolo.GetBaselineNoiseHeat(1)<CUTS['Heat']['Max']\
            and CUTS['Coll1']['Min']<bolo.GetBaselineNoiseCollectrode(1)<CUTS['Coll1']['Max']\
            and CUTS['Coll2']['Min']<bolo.GetBaselineNoiseCollectrode(2)<CUTS['Coll2']['Max']\
            and CUTS['Veto1']['Min']<bolo.GetBaselineNoiseVeto(1)<CUTS['Veto1']['Max']\
            and CUTS['Veto2']['Min']<bolo.GetBaselineNoiseVeto(2)<CUTS['Veto2']['Max']\
            and CUTS['Guard1']['Min']<bolo.GetBaselineNoiseGuard(1)<CUTS['Guard1']['Max']\
            and CUTS['Guard2']['Min']<bolo.GetBaselineNoiseGuard(2)<CUTS['Guard2']['Max']\
            and bolo.GetEventFlag()==2\
            and bolo.GetVoltageFlag>=0\
            and bolo.GetIonPulseTimeOffset()<600:
            EnergyRec = functions.recoil_energy_estimator(EnergyHeat, VOLTAGE)
            if 0 < EnergyIon < E_ION_MAX and 0 < EnergyRec < E_REC_MAX:
                Comment = ' '
                if eventcounter%20 == 0:
                    print '-'*80
                    print '{0:10} | {1:12} | {2:12} | {3:10} | {4:10} | {5:10} | {6:10}'.format('FileEntry',  'EventNumber',  'UnixTime',  'EnergyHeat',  'EnergyRec',  'EnergyIon',  'Comment')
                    print '-'*80
                if EnergyIon < 2.:
                    Comment = "Ion < 2keV"
                if bolo.GetChi2Flag() == 1:
                    TimeList.append(UnixTime)
                    IonList.append(EnergyIon)
                    RecList.append(EnergyRec)
                    eventcounter += 1
                else:
                    Comment = "failed additional cut"
                    TimeListCuts.append(UnixTime)
                    IonListCuts.append(EnergyIon)
                    RecListCuts.append(EnergyRec)
                    cutcounter += 1
                print '{0:10} | {1:12} | {2:12} | {3:10.2f} | {4:10.2f} | {5:10.2f} | {6:10}'.format(entry, EventNumber, UnixTime, EnergyHeat, EnergyRec, EnergyIon, Comment)
infile.Close()

print "%i events passed all standard cuts" % eventcounter
print "%i events failed additional cut" % cutcounter

# write events in txt-file
if OutFileName:
    OutFile = open(OutFileName, 'w')
    for i in range(len(TimeList)):
        UnixTime = TimeList[i]
        YearTime = unixtime_to_year(UnixTime)
        EnergyRec = RecList[i]
        EnergyIon = IonList[i]
        OutFile.write("%s %s %s \n" % (YearTime, EnergyRec, EnergyIon))
        OutFile.close()
    print 'Outfile %s written to disc' % OutFileName


# fill events in E_rec/E_ion-TGraph and plot them on canvas
EventGraph = TGraph()
for i in range(eventcounter):
    EnergyRec = RecList[i]
    EnergyIon = IonList[i]
    EventGraph.SetPoint(i, EnergyRec, EnergyIon)
EventGraph.GetXaxis().SetTitle('E_{rec} [keVnr]')
EventGraph.GetYaxis().SetTitle('E_{ion} [keVee]')
EventGraph.SetMarkerStyle(ROOT.kFullDotLarge)
EventGraph.SetMarkerColor(ROOT.kBlack)


# same for additional set of special cut events
EventGraphCuts = TGraph()
for i in range(cutcounter):
    EnergyRec = RecListCuts[i]
    EnergyIon = IonListCuts[i]
    EventGraphCuts.SetPoint(i, EnergyRec, EnergyIon)
EventGraphCuts.GetXaxis().SetTitle('E_{rec} [keVnr]')
EventGraphCuts.GetYaxis().SetTitle('E_{ion} [keVee]')
EventGraphCuts.SetMarkerStyle(ROOT.kFullDotLarge)
EventGraphCuts.SetMarkerColor(ROOT.kMagenta)

# wimp signal
if WIMP_MASS:
    total_efficiency = functions.efficiency_ID3(E_THRESH, FWHM_REC/2.35)
    
    signal = functions.wimp_signal(WIMP_MASS, FWHM_ION, FWHM_REC/2.35)
    signal.Multiply(total_efficiency)
    
    gamma_bckgd = functions.flat_gamma_bckgd(SIGMA_ION.getVal(), SIGMA_REC.getVal())
    gamma_bckgd.Multiply(total_efficiency)
    #gROOT.LoadMacro("/kalinka/home/hehn/PhD/LowMassEric/WimpDistriAdapted.C")
    
    #functions.TRIGGER_EFFICIENCY.SetParameter(0, E_THRESH)
    #functions.TRIGGER_EFFICIENCY.SetParameter(1, FWHM_REC)
    #functions.TRIGGER_EFFICIENCY.SetNpx(2500)
    
    ## read in WIMP spectrum
    #Signal = WimpDistri(str(WIMP_MASS), DETECTOR_NAME, FWHM_REC, FWHM_ION, TRIGGER_EFFICIENCY.GetHistogram(), 0, 0, GAMMA_CUT, VOLTAGE, 1)
    
    #Signal.SetTitle('WIMP Signal %i GeV'%WIMP_MASS)
    #Signal.GetXaxis().SetTitle('E_{recoil} [keV_{nr}]')
    #Signal.GetYaxis().SetTitle('E_{ion} [keV_{ee}]')
    #Signal.GetXaxis().SetRangeUser(3.2, 12.8)
    #Signal.GetYaxis().SetRangeUser(0.8, 5.8)
    #Signal.SetStats(0)
    #Signal.SetContour(20)


# plot everything
c1 = TCanvas('c1',  'Event selection for ID3', 600, 600)
if Signal:
    Signal.Draw('CONT0')
    EventGraph.Draw('SAMEP')
else:
    EventGraph.Draw('AP')
EventGraphCuts.Draw('SAMEP')

# add lines for Electron Recoil and Nuclear Recoil centroids
functions.ER_CENTROID.SetParameter(0, VOLTAGE)
functions.ER_CENTROID.SetLineColor(kBlue)
functions.ER_CENTROID.SetLineWidth(3)
functions.ER_CENTROID.Draw('SAME')
functions.GAMMA_CUT.SetParameter(0, VOLTAGE)
functions.GAMMA_CUT.SetParameter(2, FWHM_ION/2.35)
functions.GAMMA_CUT.SetParameter(3, FWHM_REC/2.35)
functions.GAMMA_CUT.SetParameter(4, GAMMA_CUT)
functions.GAMMA_CUT.SetLineColor(ROOT.kBlue)
functions.GAMMA_CUT.SetLineWidth(3)
functions.GAMMA_CUT.SetLineStyle(7)
functions.GAMMA_CUT.Draw('SAME')
functions.NR_CENTROID.SetLineColor(ROOT.kMagenta)
functions.NR_CENTROID.SetLineWidth(3)
functions.NR_CENTROID.SetLineStyle(7)
functions.NR_CENTROID.Draw('SAME')
