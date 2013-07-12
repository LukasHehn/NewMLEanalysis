#!/usr/bin/env python
from DetectorClass import *
from ROOT import *
import Functions


# define parameters used for skimming event set
DetectorName = 'ID3'
KDataFile = 'Data/Run12_ID3_bckg_with_subrecords.root'
OutFileName = False#'Data/ID3_Recalibration_Test.txt'
EnergyIonMax = 14
EnergyRecMax = 25
IonFactor = 1.0
HeatFactor = 1.0
Voltage = 6.4
E_Thresh = 3.874
FWHM_heat = 0.82
FWHM_ion = 0.72

# Additional WIMP signal parameters
GammaCut = 1.96
wimp_mass = 10 #if wimp mass set, show also signal density of wimp signal


CUTS = EricsLowEnergyCuts[DetectorName]
TimeList,RecList,IonList = [],[],[] #for normal event selection
TimeListCuts,RecListCuts,IonListCuts = [],[],[] #for additional cut event selection
eventcounter, cutcounter = 0, 0


infile = KDataReader(KDataFile)

# get first event and total number of events
event = infile.GetEvent()
entries = infile.GetEntries()
print entries,"total number of events in file",KDataFile

#stamp = event.GetStamp()
for entry in range(entries):
  #print "entry",entry
  infile.GetEntry(entry)
  bolos = event.GetNumBolos()

  for bolonum in range(bolos):
  #if bolos == 1:
    #bolonum = 0

    #bolo records
    bolo = event.GetBolo(bolonum)

    EnergyIon = IonFactor*bolo.GetEnergyIonFiducial()
    EnergyHeat = HeatFactor*bolo.GetEnergyHeat(1)

    #samba records
    samba = bolo.GetSambaRecord()
    
    EventNumber = samba.GetSambaEventNumber()

    UnixTime = samba.GetNtpDateSec()
    RunName = samba.GetRunName()

    if bolo.GetDetectorName()==DetectorName\
       and CUTS['Heat']['Min']<bolo.GetBaselineNoiseHeat(1)<CUTS['Heat']['Max']\
       and CUTS['Coll1']['Min']<bolo.GetBaselineNoiseCollectrode(1)<CUTS['Coll1']['Max']\
       and CUTS['Coll2']['Min']<bolo.GetBaselineNoiseCollectrode(2)<CUTS['Coll2']['Max']\
       and CUTS['Veto1']['Min']<bolo.GetBaselineNoiseVeto(1)<CUTS['Veto1']['Max']\
       and CUTS['Veto2']['Min']<bolo.GetBaselineNoiseVeto(2)<CUTS['Veto2']['Max']\
       and CUTS['Guard1']['Min']<bolo.GetBaselineNoiseGuard(1)<CUTS['Guard1']['Max']\
       and CUTS['Guard2']['Min']<bolo.GetBaselineNoiseGuard(2)<CUTS['Guard2']['Max']\
       and bolo.GetEventFlag()==2\
       and bolo.GetVoltageFlag>=0\
       and bolo.GetIonPulseTimeOffset()<600: #       and bolo.TestCutsBit(6)==True\
        EnergyRec = Functions.GetEnergyRecoilFromEstimator(EnergyHeat, Voltage)
        if 0 < EnergyIon < EnergyIonMax and 0 < EnergyRec < EnergyRecMax:
          print 'Entry: {0:6};Event: {1:6}; time: {2:8}; E_heat: {3:4.1f}; E_rec: {4:4.1f}; E_ion: {5:4.1f}'.format(entry,EventNumber, UnixTime, EnergyHeat, EnergyRec, EnergyIon)
          if EnergyIon < 2:
            print "low E_ion"
          if bolo.TestCutsBit(6) == True:
            TimeList.append(UnixTime)
            IonList.append(EnergyIon)
            RecList.append(EnergyRec)
            eventcounter += 1
          else:
            print "Add. cut failed"
            TimeListCuts.append(UnixTime)
            IonListCuts.append(EnergyIon)
            RecListCuts.append(EnergyRec)
            cutcounter += 1
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
  EventGraph.SetPoint(i,EnergyRec,EnergyIon)
EventGraph.GetXaxis().SetTitle('E_{rec} [keVnr]')
EventGraph.GetYaxis().SetTitle('E_{ion} [keVee]')
EventGraph.SetMarkerStyle(kFullDotLarge)
EventGraph.SetMarkerColor(kBlack)


# same for additional set of special cut events
EventGraphCuts = TGraph()
for i in range(cutcounter):
  EnergyRec = RecListCuts[i]
  EnergyIon = IonListCuts[i]
  EventGraphCuts.SetPoint(i,EnergyRec,EnergyIon)
EventGraphCuts.GetXaxis().SetTitle('E_{rec} [keVnr]')
EventGraphCuts.GetYaxis().SetTitle('E_{ion} [keVee]')
EventGraphCuts.SetMarkerStyle(kFullDotLarge)
EventGraphCuts.SetMarkerColor(kMagenta)

# wimp signal
if wimp_mass:
  gROOT.LoadMacro("/kalinka/home/hehn/PhD/LowMassEric/WimpDistriAdapted.C")
  
  FWHM_rec = Functions.RecoilResolutionFromHeat(FWHM_heat,Voltage,10)
    
  TriggerEfficiency.SetParameter(0, E_Thresh)
  TriggerEfficiency.SetParameter(1, 1.6)
  TriggerEfficiency.SetNpx(2500)

  # read in WIMP spectrum
  Signal = WimpDistri(str(wimp_mass), DetectorName, FWHM_rec, FWHM_ion, TriggerEfficiency.GetHistogram(), 0, 0, GammaCut, Voltage, 1)
  
  Signal.SetTitle('WIMP Signal %i GeV'%wimp_mass)
  Signal.GetXaxis().SetTitle('E_{recoil} [keV_{nr}]')
  Signal.GetYaxis().SetTitle('E_{ion} [keV_{ee}]')
  Signal.GetXaxis().SetRangeUser(3.2,12.8)
  Signal.GetYaxis().SetRangeUser(0.8,5.8)
  Signal.SetStats(0)
  Signal.SetContour(20)


# plot everything
c1 = TCanvas('c1','Event selection for ID3',600,600)
if Signal:
  Signal.Draw('CONT0')
  EventGraph.Draw('SAMEP')
else:
  EventGraph.Draw('AP')
EventGraphCuts.Draw('SAMEP')

# add lines for Electron Recoil and Nuclear Recoil centroids
Functions.ER_centroid.SetParameter(0,Voltage)
Functions.ER_centroid.Draw('SAME')
Functions.ER_centroid.SetLineColor(kBlue)
Functions.ER_centroid.SetLineWidth(3)
Functions.NR_centroid.Draw('SAME')
Functions.NR_centroid.SetLineColor(kBlue)
Functions.NR_centroid.SetLineWidth(3)
Functions.NR_centroid.SetLineStyle(7)
