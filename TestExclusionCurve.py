#!/usr/bin/env python
from ROOT import *
from Parameters import *
from DetectorClass import *
from Functions import *
import pyWIMP.DMModels.wimp_model as wimp_model
from pyWIMP.DMModels.flat_model import FlatModel
from pyWIMP.DMModels.base_model import BaseVariables


def WIMPsignal(mass_of_wimp,sigma_ion,sigma_rec,spectrum):
  denom_i=1./(2*pow(sigma_ion,2))
  denom_r=1./(2*pow(sigma_rec,2))
  signal_hist = TH2F('signal_hist_%sGeV'%mass_of_wimp,'WIMP Density for %s GeV;E_{Rec} (keV_{nr});E_{Ion} (keV_{nr});Density'%mass_of_wimp,Energy['rec']['bins'],Energy['rec']['min'],Energy['rec']['max'],Energy['ion']['bins'],Energy['ion']['min'],Energy['ion']['max'])
  for recbin in range(1,signal_hist.GetNbinsX()+1):
    Erec = signal_hist.GetXaxis().GetBinCenter(recbin)
    for ionbin in range(1,signal_hist.GetNbinsY()+1):
      Eion = signal_hist.GetYaxis().GetBinCenter(ionbin)
      summe=0
      for recbinspec in range(1,spectrum.GetNbinsX()+1):
	ErecSpec = spectrum.GetXaxis().GetBinCenter(recbinspec)
	Q = LindhardQuenching.Eval(ErecSpec)
	wimprate = spectrum.GetBinContent(recbinspec)
	kernel=TMath.Exp(-denom_r*pow((Erec-ErecSpec),2)-denom_i*pow((Eion-Q*ErecSpec),2))
	summe+=(kernel*wimprate)
      signal_hist.SetBinContent(recbin,ionbin,0.1*summe/(2*3.141592*sigma_rec*sigma_ion))
  return signal_hist


def Perform(mass_of_wimp):
  Results[mass_of_wimp] = {}
  result = Results[mass_of_wimp]

  Container[mass_of_wimp] = {}
  this = Container[mass_of_wimp]

  this['pywimp_pdf'] = wm.get_WIMP_model_with_escape_vel(mass_of_wimp)

  result['f_norm'] = wm.get_normalization().getVal()
  result['N_wimp'] = this['pywimp_pdf'].expectedEvents(RooArgSet(time,rec)) # number of expected events per nucleon cross section [pb] for detector mass in time and energy range

  this['pywimp_hist'] = this['pywimp_pdf'].createHistogram('pywimp_hist_%sGeV'%mass_of_wimp,rec,RooFit.Binning(Energy['rec']['bins']))
  # scale to differential rate in nucleon cross section
  this['pywimp_hist'].Scale(result['N_wimp']/this['pywimp_hist'].Integral())

  # signal histogram
  this['wimp_hist'] = WIMPsignal(mass_of_wimp,sigma_ion,sigma_rec,this['pywimp_hist'])

  # scale to expected WIMP event count
  this['wimp_hist'].Scale(result['N_wimp']/this['wimp_hist'].Integral())
  this['signal_hist'] = this['wimp_hist'].Clone('signal_hist_%sGeV'%mass_of_wimp)

  # correct for detector efficiency
  this['signal_hist'].Multiply(total_efficiency)

  # correction for  additional gamma cut
  this['signal_hist'].Multiply(gamma_cut_efficiency)

  # remaining number of wimp events after efficiency correction
  result['N_signal'] = this['signal_hist'].Integral()

  # calculate cross section limit for no observed events and 90% poisson error
  result['cross section limit'] = result['f_norm']*(2.35/result['N_signal'])
  return True


DetectorName = 'ID3'
ID = Detector(DetectorName)
total_efficiency = ID.GetProjectedEfficiency('total')
gamma_cut_efficiency = ID.GetGammaCutEfficiency()


# list of WIMP masses
WIMP_Masses = [7,8,10,12,15,20,25,30]

# dictionaries
Container = {}
Results = {}

# pyWIMP library initialization
basevars = BaseVariables(ID.GetTimeBinning()[0],ID.GetTimeBinning()[-1],Energy['rec']['min'],Energy['rec']['max'],True,time_offset=-2./365.25) #offset to start at 1st Mar 2009
wm = wimp_model.WIMPModel(basevars,WIMP_Masses[0],ID.GetMass(),constant_quenching=True,nucl_recoil=True)

# definition of observables
ion = RooRealVar('ion','E_{ion}',Energy['ion']['min'],Energy['ion']['max'],'keV_{ee}')
time = basevars.get_time()
rec = basevars.get_energy()
rec.SetNameTitle('rec','E_{rec}')
rec.setUnit('keV_{nr}')

# detector specific parameters
voltage = ID.GetWeightedAverage('voltage')

FWHM_heat = ID.GetWeightedAverage('heat')
sigma_heat = FWHM_heat/2.35

sigma_rec = RecoilResolutionFromHeatBaseline(sigma_heat,voltage,10)

FWHM_ion = ID.GetWeightedAverage('fiducialmean')
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
  print mass_of_wimp,r['N_wimp'],r['N_signal'],r['f_norm'],r['cross section limit']
