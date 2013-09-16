#!/usr/bin/env python

####################################################################################################
#
# Filter and display events from a KData file
# Lukas Hehn, 2013
#
####################################################################################################

import functions
import parameters
import ROOT

from parameters import BASELINE_CUTS_ERIC
from ROOT import TGraph, gROOT, TCanvas, KDataReader


# Definition of parameters used for skimming event set
DETECTOR_NAME = 'ID401'
#KDataFile = 'Data/Run12_ID3_bckg_with_subrecords.root'
KDataFile = 'Data/Run12_ID3+6+401+404_bckg_with_subrecords.root'
OutFileName = False  # 'Data/ID3_eventlist_lowE_corrected.txt'
WIMP_MASS = 8  # if wimp mass set, show also signal density of wimp signal


# Difinition of energy range and binning
E_ION_MAX = 14.
E_REC_MAX = 25.
BINSIZE = 0.1


# Use energy correction with values defined below
ENERGY_CORRECTION = False
IonFactor = 1.0/1.029
HeatFactor = 1.0/1.017


# Import values for baseline, voltage and threshold from parameter dictionary
VOLTAGE = parameters.AVG_VOLTAGES_ERIC[DETECTOR_NAME]
E_THRESH_EE = parameters.E_THRESHOLD_ERIC[DETECTOR_NAME]
E_THRESH_NR = functions.recoil_energy_estimator(E_THRESH_EE, VOLTAGE)
FWHM_HEAT = parameters.ENERGY_RESOLUTIONS_ERIC[DETECTOR_NAME]['Heat']
FWHM_ION = parameters.ENERGY_RESOLUTIONS_ERIC[DETECTOR_NAME]['Fiducial']
FWHM_REC = functions.fwhm_rec_from_heat(FWHM_HEAT, VOLTAGE, 10)
GAMMA_CUT = 1.96  # gamma cut from ER centroid in standard deviations
print 'Voltage={voltage}V, threshold={thresh:.2f}keVnr, FWHM_ion={ion}keVee, FWHM_heat={heat}keVee, FWHM_rec={rec:.2f}keVnr\n'.format(
    voltage=VOLTAGE, thresh=E_THRESH_NR, ion=FWHM_ION, heat=FWHM_HEAT, rec=FWHM_REC)

# Get dictionary with cuts for all baselines
CUTS = BASELINE_CUTS_ERIC[DETECTOR_NAME]


# Define empty lists for event storage as well as reset event counters
TimeList, RecList, IonList = [], [], []  # for events passing all cuts
TimeListCuts, RecListCuts, IonListCuts = [], [], []  # for events failing additional cut
eventcounter, cutcounter = 0, 0


# Here events are selected from a KData file
infile = KDataReader(KDataFile)

# Get first event and total event number
event = infile.GetEvent()
entries = infile.GetEntries()
print entries, "total number of events in file", KDataFile, "\n"

for entry in range(entries):
    infile.GetEntry(entry)
    bolos = event.GetNumBolos()

    for bolonum in range(bolos):
    #if bolos == 1:
       #bolonum = 0

        # Bolometer records
        bolo = event.GetBolo(bolonum)
        if ENERGY_CORRECTION:
            EnergyIon = IonFactor*bolo.GetEnergyIonFiducial()
            EnergyHeat = HeatFactor*bolo.GetEnergyHeat(1)
        else:
            EnergyIon = bolo.GetEnergyIonFiducial()
            EnergyHeat = bolo.GetEnergyHeat(1)

        # Samba records
        samba = bolo.GetSambaRecord()
        EventNumber = samba.GetSambaEventNumber()
        UnixTime = samba.GetNtpDateSec()
        RunName = samba.GetRunName()

        if bolo.GetDetectorName() == DETECTOR_NAME\
           and CUTS['Heat']['Min'] < bolo.GetBaselineNoiseHeat(1) < CUTS['Heat']['Max']\
           and CUTS['Coll1']['Min'] < bolo.GetBaselineNoiseCollectrode(1) < CUTS['Coll1']['Max']\
           and CUTS['Coll2']['Min'] < bolo.GetBaselineNoiseCollectrode(2) < CUTS['Coll2']['Max']\
           and CUTS['Veto1']['Min'] < bolo.GetBaselineNoiseVeto(1) < CUTS['Veto1']['Max']\
           and CUTS['Veto2']['Min'] < bolo.GetBaselineNoiseVeto(2) < CUTS['Veto2']['Max']\
           and CUTS['Guard1']['Min'] < bolo.GetBaselineNoiseGuard(1) < CUTS['Guard1']['Max']\
           and CUTS['Guard2']['Min'] < bolo.GetBaselineNoiseGuard(2) < CUTS['Guard2']['Max']\
           and bolo.GetEventFlag() == 2\
           and bolo.GetVoltageFlag >= 0\
           and bolo.GetIonPulseTimeOffset() < 600:

            EnergyRec = functions.recoil_energy_estimator(EnergyHeat, VOLTAGE)  # estimate recoil energy

            if 0 < EnergyIon < E_ION_MAX and 0 < EnergyRec < E_REC_MAX:

                if eventcounter%20 == 0:
                    print '-'*100
                    print '{0:10} | {1:12} | {2:12} | {3:10} | {4:10} | {5:10} | {6:10}'.format('FileEntry', 'EventNumber', 'UnixTime', 'EnergyHeat', 'EnergyRec', 'EnergyIon', 'Comment')
                    print '-'*100

                if EnergyIon < 2.:
                    Comment = 'Ion < 2keV'
                else:
                    Comment = ' '

                if bolo.GetChi2Flag():
                    TimeList.append(UnixTime)
                    IonList.append(EnergyIon)
                    RecList.append(EnergyRec)
                    eventcounter += 1
                else:
                    Comment = "failed GetChi2Flag"
                    TimeListCuts.append(UnixTime)
                    IonListCuts.append(EnergyIon)
                    RecListCuts.append(EnergyRec)
                    cutcounter += 1

                print '{0:10} | {1:12} | {2:12} | {3:10.2f} | {4:10.2f} | {5:10.2f} | {6:10}'.format(entry, EventNumber, UnixTime, EnergyHeat, EnergyRec, EnergyIon, Comment)
infile.Close()

print '{events} events passed all standard cuts'.format(events=eventcounter)
print '{events} events failed additional cut'.format(events=cutcounter)


# Write events in txt-file
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


# Fill events in E_rec/E_ion-TGraph and plot them on canvas
EventGraph = TGraph()
for i in range(eventcounter):
    EnergyRec = RecList[i]
    EnergyIon = IonList[i]
    EventGraph.SetPoint(i, EnergyRec, EnergyIon)
EventGraph.GetXaxis().SetTitle('E_{rec} (keVnr)')
EventGraph.GetYaxis().SetTitle('E_{ion} (keVee)')
EventGraph.SetMarkerStyle(ROOT.kFullDotLarge)
EventGraph.SetMarkerColor(ROOT.kBlack)


# Same for additional set of special cut events
EventGraphCuts = TGraph()
for i in range(cutcounter):
    EnergyRec = RecListCuts[i]
    EnergyIon = IonListCuts[i]
    EventGraphCuts.SetPoint(i, EnergyRec, EnergyIon)
EventGraphCuts.GetXaxis().SetTitle('E_{rec} (keVnr)')
EventGraphCuts.GetYaxis().SetTitle('E_{ion} (keVee)')
EventGraphCuts.SetMarkerStyle(ROOT.kFullDotLarge)
EventGraphCuts.SetMarkerColor(ROOT.kMagenta)


# Get efficiency and flat gamma background
total_efficiency = functions.simple_efficiency(DETECTOR_NAME, E_THRESH_NR, FWHM_REC/2.35,
                                               binsize=BINSIZE, rec_max=E_REC_MAX, ion_max=E_ION_MAX
                                               )

gamma_bckgd = functions.flat_gamma_bckgd(FWHM_ION/2.35, FWHM_REC/2.35, VOLTAGE,
                                         binsize=BINSIZE, rec_max=E_REC_MAX, ion_max=E_ION_MAX
                                         )
gamma_bckgd.Multiply(total_efficiency)
gamma_bckgd.Scale(1./gamma_bckgd.GetMaximum())  # maximum value is 1 afterwards
gamma_bckgd.SetStats(0)


# Also get wimp signal PDF
if WIMP_MASS:
    signal = functions.wimp_signal(WIMP_MASS, FWHM_ION/2.35, FWHM_REC/2.35,
                                   binsize=BINSIZE, rec_max=E_REC_MAX, ion_max=E_ION_MAX
                                   )
    signal.Multiply(total_efficiency)
    signal.Scale(1./signal.GetMaximum())  # maximum value is 1 afterwards
    signal.SetStats(0)
    print 'Created 2D WIMP signal for %i GeV mass'%WIMP_MASS

    bckgd_plus_sig = gamma_bckgd.Clone('bckgd_plus_sig')
    bckgd_plus_sig.Add(signal)
    bckgd_plus_sig.SetStats(0)
    bckgd_plus_sig.SetContour(30)
    bckgd_plus_sig.SetTitle('{detector} Flat Gamma Bckgd + {mass}GeV WIMP signal'.format(detector=DETECTOR_NAME, mass=WIMP_MASS))
    print "Combined flat gamma background and WIMP signal both normalized to maximum of 1."


# Plot PDF, data and ER and NR centroids
c1 = TCanvas('c1', 'Event selection for {detector}'.format(detector=DETECTOR_NAME), 600, 600)
if signal:
    bckgd_plus_sig.Draw('CONT0Z')
else:
    gamma_bckgd.Draw('CONT0Z')
EventGraph.Draw('SAMEP')
#EventGraphCuts.Draw('SAMEP')

functions.ER_CENTROID.SetParameter(0, VOLTAGE)
functions.ER_CENTROID.SetLineColor(ROOT.kBlue)
functions.ER_CENTROID.SetLineWidth(3)
functions.ER_CENTROID.Draw('SAME')
functions.NR_CENTROID.SetLineColor(ROOT.kBlack)
functions.NR_CENTROID.SetLineWidth(3)
functions.NR_CENTROID.SetLineStyle(1)
functions.NR_CENTROID.Draw('SAME')

# Plot Gamma cut line
if GAMMA_CUT:
    functions.GAMMA_CUT.FixParameter(0, VOLTAGE)
    functions.GAMMA_CUT.FixParameter(2, FWHM_ION/2.35)
    functions.GAMMA_CUT.FixParameter(3, FWHM_REC/2.35)
    functions.GAMMA_CUT.FixParameter(4, GAMMA_CUT)
    functions.GAMMA_CUT.SetLineColor(ROOT.kBlue)
    functions.GAMMA_CUT.SetLineWidth(3)
    functions.GAMMA_CUT.SetLineStyle(7)
    functions.GAMMA_CUT.Draw('SAME')
