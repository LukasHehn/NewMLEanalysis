#!/usr/bin/env python
from DetectorClass import *
from ROOT import *
from Parameters import *

ID = Detector('ID3')

KDataFile = 'Run12_ID3_bckg_with_subrecords.root'

def FindKDataEvents(ID, KDataFile):
  eventcounter = 0
  infile = KDataReader(KDataFile)
  event = infile.GetEvent()
  entries = infile.GetEntries()
  stamp = event.GetStamp()
  for entry in range(entries):
    infile.GetEntry(entry)
    bolos = event.GetNumBolos()

    #for bolonum in range(bolos):
    if bolos == 1:
      bolonum = 0

      #bolo records
      bolo = event.GetBolo(bolonum)

      EnergyIon = bolo.GetEnergyIonFiducial()

      #samba records
      samba = bolo.GetSambaRecord()

      UnixTime = samba.GetNtpDateSec()
      RunName = samba.GetRunName()

      if bolo.GetDetectorName()==ID.GetName()\
         and bolo.GetEventFlag()==2\
         and bolo.TestCutsBit(6)==True\
         and bolo.GetVoltageFlag>=0\
         and bolo.GetIonPulseTimeOffset()<600:
        EnergyRec = ID.GetNuclearRecoilEnergy(bolo.GetEnergyHeat(1), bolo.GetVoltageFlag())
        if ID.IsGoodEvent(UnixTime, EnergyRec, EnergyIon) == True:
          ID.AddEvent(UnixTime, EnergyRec, EnergyIon)
          print 'time: {0:8}, E_rec: {1:4.1f}, E_ion: {2:4.1f}'.format(UnixTime, EnergyRec, EnergyIon)
          eventcounter += 1

  infile.Close()

  return "KData: %i events read" % eventcounter

  # save events to txt-file
  #ID.WriteEventsToFile('ID3-eventlist_new.txt')
