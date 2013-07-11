#!/usr/bin/env python
from DetectorClass import *
from ROOT import *
import Functions


# define parameters used for skimming event set
DetectorName = 'ID3'
KDataFile = 'Data/Run12_ID3_bckg_with_subrecords.root'
OutFileName = 'Data/ID3_Recalibration_Test.txt'
EnergyIonMax = 14
EnergyRecMax = 25
IonFactor = 1.0
HeatFactor = 1.0
Voltage = 6.4

CUTS = EricsLowEnergyCuts['ID3']
TimeList,RecList,IonList = [],[],[]
eventcounter = 0

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

  #for bolonum in range(bolos):
  if bolos == 1:
    bolonum = 0

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
       and bolo.TestCutsBit(6)==True\
       and bolo.GetIonPulseTimeOffset()<600:
        EnergyRec = Functions.GetEnergyRecoilFromEstimator(EnergyHeat, Voltage)
        if 0 < EnergyIon < EnergyIonMax and 0 < EnergyRec < EnergyRecMax:
          TimeList.append(UnixTime)
          IonList.append(EnergyIon)
          RecList.append(EnergyRec)
          print 'Event: {0:6}, time: {1:8}, E_heat: {2:4.1f}, E_rec: {3:4.1f}, E_ion: {4:4.1f}'.format(EventNumber, UnixTime, EnergyHeat, EnergyRec, EnergyIon)
          eventcounter += 1

infile.Close()

print "KData: %i events passed all cuts" % eventcounter

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

c1 = TCanvas('c1','Event selection for ID3',800,600)
EventGraph.SetMarkerStyle(kPlus)
EventGraph.Draw('AP')

# add lines for Electron Recoil and Nuclear Recoil centroids
Functions.ER_centroid.SetParameter(0,6.4)
Functions.ER_centroid.Draw('SAME')
Functions.NR_centroid.Draw('SAME')
