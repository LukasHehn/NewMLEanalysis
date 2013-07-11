#!/usr/bin/env python
import datetime as dt
import calendar
#from ROOT import TF1, exp, log, TH1F, TH2F
from ROOT import *
from Parameters import *
from couchdbkit import Server, Database
import couchdbkit
import numpy as np


# global
global inpath; inpath = 'Run12PeriodInformation/'


# date and time binning functions
March1st = dt.datetime(2009, 3, 1, 0, 0, 0)
OffsetJules = dt.timedelta(hours=-1)
global dayzero; DayZero = March1st + OffsetJules


def unixtime_to_year(unixtime):
  date = dt.datetime.utcfromtimestamp(unixtime)
  delta = date - DayZero
  return (delta.days/365.24)+(delta.seconds/(365.24*86400))


def unixtime_to_day(unixtime):
  date = dt.datetime.utcfromtimestamp(unixtime)
  delta = date - DayZero
  return delta.days


def unixtime_to_dayhour(unixtime):
  date = dt.datetime.utcfromtimestamp(unixtime)
  delta = date - DayZero
  return delta.days + delta.seconds/86400.


def unixtime_to_hour(unixtime):
  date = dt.datetime.utcfromtimestamp(unixtime)
  delta = date - DayZero
  return float(delta.seconds)/3600


def unixtime_to_date(unixtime):
  date = dt.datetime.utcfromtimestamp(unixtime)
  return date


def day_to_date(days):
  delta = dt.timedelta(days)
  date = DayZero + delta
  return date


def date_to_day(date):
  delta = date - DayZero
  return delta.days


def day_to_unixtime(day):
  date = day_to_date(day)
  unixtime = calendar.timegm(date.timetuple())
  return unixtime

def date_to_unixtime(date):
  unixtime = calendar.timegm(date.timetuple())
  return unixtime


# root functions
# trigger efficiency from DAQ trigger threshold and resolution on heat channel
TriggerEfficiency = TF1('trigger_efficiency','0.5*(1+ROOT::Math::erf(((x-[0])/([1]*sqrt(2)))))', 0, 20)
TriggerEfficiency.SetNpx(1000)
TriggerEfficiency.SetParName(0, 'Threshold')
TriggerEfficiency.SetParName(1, 'Resolution')
TriggerEfficiency.SetTitle('Trigger Efficiency;E_{Heat} (keV);Efficiency')


# measured fiducial efficiency (division of neutron histogram with/without fiducial cut)
FiducialEfficiency = TF1('fiducial_efficiency','[2]*(1-exp([0]*(x-[1])))', 0, 30)
FiducialEfficiency.SetNpx(1000)
FiducialEfficiency.SetParName(0,'Slope')
FiducialEfficiency.SetParName(1,'Cut off')
FiducialEfficiency.SetParName(2,'Maximum')
FiducialEfficiency.SetTitle('Fiducial Efficiency;E_{ion} (keV_{ee});Efficiency')


# Lindhard quenching relation (nuclear recoil)
LindhardQuenching = TF1('lindhard_quenching','[0]*(x^[1])', 0, 30)
LindhardQuenching.SetNpx(1000)
LindhardQuenching.SetParName(0,'a')
LindhardQuenching.SetParName(1,'b')
LindhardQuenching.FixParameter(0, 0.16)
LindhardQuenching.FixParameter(1, 0.18)
LindhardQuenching.SetTitle('Lindhard Quenching for Nuclear Recoils;E_{Recoil} [keV];Q(E_{Rec})')


# Recoil energy estimator for nuclear recoils
RecoilEstimator = TF1('recoil_energy_estimator','(x/(1+[0]/[1]))*(1+[0]/[1]*0.16*x^0.18)', 0, 30)
RecoilEstimator.SetNpx(1000)
RecoilEstimator.SetParName(0,'Voltage')
RecoilEstimator.SetParName(1,'Creation Potential')
RecoilEstimator.FixParameter(1,3.0)
RecoilEstimator.SetTitle('E_{Rec} estimator from E_{Heat};E_{Rec} [keV_{nr}];E_{Heat} [keV_{ee}]')


# centroids of ER and NR band
def ER_centroid_function(x, par):
      xx =x[0]
      f = xx*(1+(0.16*xx**(0.18+0j))*(par[0]/3))/(1+par[0]/3)
      return f.real
ER_centroid_real = TF1("ER_centroid_real", ER_centroid_function, -10, 30, 1)
ER_centroid_real.SetNpx(1000)
ER_centroid_real.SetParameter(0, 6.4)
ER_centroid_real.SetTitle('Electron Recoil Centroid (Real Part Only);E_{Rec} [keV_{nr}];E_{Heat} [keV_{ee}]')

ER_centroid = TF1('ER_centroid','x*(1+(0.16*x^0.18)*([0]/3))/(1+[0]/3)',0,30)
ER_centroid.SetNpx(1000)
ER_centroid.SetParName(0,'voltage')
ER_centroid_real.SetTitle('Electron Recoil Centroid;E_{Rec} [keV_{nr}];E_{Heat} [keV_{ee}]')

NR_centroid = TF1('NR_centroid','0.16*x^1.18',0,30)
NR_centroid.SetNpx(1000)
NR_centroid.SetTitle('Nuclear Recoil Centroid;E_{Rec} [keV_{nr}];E_{Heat} [keV_{ee}]')


def RecoilResolutionFromHeat(FWHM_heat,voltage,Erec):
  Q = LindhardQuenching.Eval(Erec)
  FWHM_rec = FWHM_heat * ((1+voltage/3.0)/(1+1.18*Q*voltage/3.0))
  return FWHM_rec


def GetEnergyRecoilFromEstimator(energy_ee,voltage):
  function = RecoilEstimator
  function.SetParameter(0,voltage)
  function.FixParameter(1,3.0) #electron-hole creation potential
  energy_recoil = function.GetX(energy_ee)
  return energy_recoil


def GetVoltageFlagList(Detector, Runtype):
  RunName, VoltageFlag = [],[]
  infile = open(inpath+Detector+'_voltage_'+Runtype+'.txt','r')
  for line in infile:
    runname, voltage = line.split()
    RunName.append(str(runname))
    VoltageFlag.append(int(voltage))
  infile.close()
  return [RunName, VoltageFlag]


def GetVoltageFlag(InList, RunName):
  try:
    flag = InList[1][InList[0].index(RunName)]
    return flag
  except ValueError:
    print "error finding run",RunName
    return 0


def GetSambaNumber(Detector):
  if Detector in ['ID5','ID403','FID401']:
    return 'S1'
  elif Detector in ['ID3','ID401','FID402']:
    return 'S2'
  elif Detector in ['ID4','ID6','ID404']:
    return 'S3'
  elif Detector in ['ID2','ID402','ID405']:
    return 'S4'


def GetERAConstants(Detector):
  s = couchdbkit.Server('https://edelweissuser:edwdbw1mp@edelweiss.cloudant.com')
  db = s['analysis']
  vr = db.view('constants/run12', include_docs=True)
  constantDoc = {}
  for row in vr:
    constantDoc[row['key']] = row['doc']
  return constantDoc[Detector]


def ReadInWimpSpectrumEric(mass_of_wimp):
  basedir="/kalinka/home/hehn/PhD/LowMassEric/"
  filename=basedir+"wimpsignal_M"+str(mass_of_wimp)+".txt"

  Energies,Rates = [], []

  #print "reading WIMP spectrum from:",filename

  infile = open(filename,'r')
  for line in infile:
    energy, rate = line.split()
    Energies.append(float(energy))
    Rates.append(float(rate))
  infile.close()

  Energybins = np.array(Energies, dtype=np.float)

  Hist = TH1F('wimp_spectrum_%sGeV'%mass_of_wimp,'Spectrum Eric %sGeV;E_{recoil} [keV];Rate [evts/kg/day/0.02keV]'%mass_of_wimp,Energybins.size-1, Energybins.flatten('C'))

  for i in range(len(Energies)):
    Hist.Fill(Energies[i], Rates[i])

  return Hist


def WimpSignal2DEric(mass_of_wimp,sigma_ion,sigma_rec,spectrum):
  denom_i=1./(2*pow(sigma_ion,2))
  denom_r=1./(2*pow(sigma_rec,2))
  hist = TH2F('wimp_signal_%sGeV'%mass_of_wimp,'WIMP signal %sGeV;E_{rec} (keVnr);E_{ion} (keVee);Rate (cts/kg*day)'%mass_of_wimp,Energy['rec']['bins'],Energy['rec']['min'],Energy['rec']['max'],Energy['ion']['bins'],Energy['ion']['min'],Energy['ion']['max'])
  for recbin in range(1,hist.GetNbinsX()+1):
    Erec = hist.GetXaxis().GetBinCenter(recbin)
    for ionbin in range(1,hist.GetNbinsY()+1):
      Eion = hist.GetYaxis().GetBinCenter(ionbin)
      summe=0
      for specbin in range(1,spectrum.GetNbinsX()+1):
	ErecSpec = spectrum.GetXaxis().GetBinCenter(specbin)
	Q = LindhardQuenching.Eval(ErecSpec)
	wimprate = spectrum.GetBinContent(specbin)
	tutu = denom_r*pow((Erec-ErecSpec),2)
	kernel = TMath.exp(-tutu-denom_i*pow((Eion-Q*ErecSpec),2))
	summe += (kernel*wimprate)
      hist.SetBinContent(recbin,ionbin,0.02002*summe/(2*3.141592*sigma_rec*sigma_ion))
  return hist


def FlatGammaBckgd2DEric(sigma_ion,sigma_rec):
  denom_i=1./(2*sigma_ion**2)
  denom_r=1./(2*sigma_rec**2)
  spectrum = TH1F('flat_gamma_spectrum','flat gamma spectrum',Energy['rec']['bins'],Energy['rec']['min'],Energy['rec']['max'])
  hist = TH2F('flat_gamma_bckgd','Flat Gamma Bckgd;E_{rec} (keVnr);E_{ion} (keVee);Rate (cts/kg*day)',Energy['rec']['bins'],Energy['rec']['min'],Energy['rec']['max'],Energy['ion']['bins'],Energy['ion']['min'],Energy['ion']['max'])
  for recbin in range(1,hist.GetNbinsX()+1):
    Erec = hist.GetXaxis().GetBinCenter(recbin)
    for ionbin in range(1,hist.GetNbinsY()+1):
      Eion = hist.GetYaxis().GetBinCenter(ionbin)
      summe=0
      for specbin in range(1,spectrum.GetNbinsX()+1):
	ErecSpec = spectrum.GetXaxis().GetBinCenter(specbin)
	EionSpec = ER_centroid_real.Eval(ErecSpec)
	tutu = denom_r*(Erec-ErecSpec)**2
	kernel = TMath.exp(-tutu-denom_i*(Eion-EionSpec)**2)
	summe += kernel
      hist.SetBinContent(recbin,ionbin,summe)
  return hist


def Simple2DEfficiencyID3():
  hist = TH2F('efficiency','efficiency;E_{rec} (keVnr);E_{ion} (keVee);Efficiency',Energy['rec']['bins'],Energy['rec']['min'],Energy['rec']['max'],Energy['ion']['bins'],Energy['ion']['min'],Energy['ion']['max'])

  TriggerEfficiency.SetParameter(0, 3.874)
  TriggerEfficiency.SetParameter(1, 0.82)
  FiducialEfficiency.SetParameter(0, -1.876)
  FiducialEfficiency.SetParameter(1, 1.247)
  FiducialEfficiency.SetParameter(2, 0.947)

  for recbin in range(1,hist.GetNbinsX()+1):
    Erec = hist.GetXaxis().GetBinCenter(recbin)
    eff_trigger = TriggerEfficiency.Eval(Erec)
    for ionbin in range(1,hist.GetNbinsY()+1):
      Eion = hist.GetYaxis().GetBinCenter(ionbin)
      eff_fiducial = FiducialEfficiency.Eval(Eion)

      if eff_trigger >= 0 and eff_fiducial >= 0:
	eff_total = eff_trigger * eff_fiducial
      else:
	eff_total = 0
      hist.SetBinContent(recbin,ionbin,eff_total)
      hist.SetBinError(recbin,ionbin,0)
  return hist


def GetGammaCutEfficiency(self,voltage):
  GammaCutEfficiency = TH2F('gamma_cut_efficiency','Gamma Cut Efficiency;E_{Rec} (keV_{nr});E_{Ion} (keV_{nr});Efficiency',200,0,20,100,0,10)
  hist = GammaCutEfficiency
  ER_centroid.SetParameter(0,voltage)
  for xbin in range(1,hist.GetNbinsX()+1):
    Erec = hist.GetXaxis().GetBinCenter(xbin)
    ioncut = ER_centroid.Eval(Erec)-offset
    for ybin in range(1,hist.GetNbinsY()+1):
      Eion = hist.GetYaxis().GetBinCenter(ybin)
      if Eion < ioncut:
	hist.SetBinContent(xbin, ybin, 1.0)
      else:
	hist.SetBinContent(xbin, ybin, 0.0)
  return True


def TGraphFromDataSet(dataset):
  entries = dataset.numEntries()
  tgraph = TGraph()
  for event in range(entries):
    time = dataset.get(event).getRealValue('time')
    rec = dataset.get(event).getRealValue('rec')
    ion = dataset.get(event).getRealValue('ion')
    #print event,time,rec,ion
    tgraph.SetPoint(event,rec,ion)
  return tgraph