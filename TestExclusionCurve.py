#!/usr/bin/env python
from ROOT import *
from Parameters import *
from DetectorClass import *
from Functions import *


def Perform(mass_of_wimp):
  Results[mass_of_wimp] = {}
  result = Results[mass_of_wimp]

  Container[mass_of_wimp] = {}
  this = Container[mass_of_wimp]

  # spectrum & signal histogram
  this['spectrum_list'], this['spectrum_hist'], this['integral'] = ReadInWimpSpectrum(mass_of_wimp)
  this['wimp_hist'] = WimpSignal2DEric(mass_of_wimp,sigma_ion,sigma_rec,this['spectrum_hist'])

  result['N_wimp'] = this['integral']*exposure

  # scale to expected WIMP event count
  this['wimp_hist'].Scale(result['N_wimp']/this['wimp_hist'].Integral('WIDTH'))
  this['signal_hist'] = this['wimp_hist'].Clone('signal_hist_%sGeV'%mass_of_wimp)

  # correct for detector efficiency
  this['signal_hist'].Multiply(total_efficiency)

  # correction for  additional gamma cut
  #this['signal_hist'].Multiply(gamma_cut_efficiency)

  # remaining number of wimp events after efficiency correction
  result['N_signal'] = this['signal_hist'].Integral('WIDTH')

  # calculate cross section limit for no observed events and 90% poisson error
  result['cross section limit'] = (2.35/result['N_signal'])*pow(10,-6)
  return True


DetectorName = 'ID3'
ID = Detector(DetectorName)
total_efficiency = ID.GetProjectedEnergyEfficiency()
gamma_cut_efficiency = ID.GetGammaCutEfficiency()
exposure = ID.GetExposure()


# list of WIMP masses
WIMP_Masses = [8,10,15,30]
#WIMP_Masses = [10]

# dictionaries
Container = {}
Results = {}

# definition of observables
ion = RooRealVar('ion','E_{ion}',Energy['ion']['min'],Energy['ion']['max'],'keV_{ee}')
rec = RooRealVar('rec','E_{rec}',Energy['rec']['min'],Energy['rec']['max'],'keV_{nr}')
time = RooRealVar('time','time',ID.GetTimeBinning()[0],ID.GetTimeBinning()[-1],'years')


# detector specific parameters
voltage = ID.GetAverageValue('voltage')

FWHM_heat = ID.GetAverageValue('heat')
sigma_heat = FWHM_heat/2.35

sigma_rec = RecoilResolutionFromHeatBaseline(sigma_heat,voltage,10)

FWHM_ion = ID.GetAverageValue('fiducialmean')
sigma_ion = FWHM_ion/2.35

print sigma_ion, sigma_rec

# dataset
realdata = RooDataSet.read('ID3-eventlist_test.txt',RooArgList(time,rec,ion))
realdata_scatter = realdata.createHistogram(rec,ion)
events = int(realdata.numEntries())

# loop over several wimp masses
for mass_of_wimp in WIMP_Masses:
  Perform(mass_of_wimp)
  r = Results[mass_of_wimp]
  print mass_of_wimp,r['N_wimp'],r['N_signal'],r['N_signal']/r['N_wimp'],r['cross section limit']
