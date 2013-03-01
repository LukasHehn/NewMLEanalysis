#!/usr/bin/env python
from ROOT import *
from Parameters import *
from DetectorClass import *
import pyWIMP.DMModels.wimp_model as wimp_model
from pyWIMP.DMModels.flat_model import FlatModel
from pyWIMP.DMModels.base_model import BaseVariables


def Perform(mass_of_wimp):
  Results[mass_of_wimp] = {}
  result = Results[mass_of_wimp]

  Container[mass_of_wimp] = {}
  this = Container[mass_of_wimp]

  # position of bands
  this['NR_centroid'] = RooFormulaVar('NR_centroid_%sGeV'%mass_of_wimp,'0.16*@0^1.18',RooArgList(rec))

  # gaussians in ionization energy with shifting mean in recoil energy
  this['ion_gauss_NR'] = RooGaussian('ion_gauss_NR_%sGeV'%mass_of_wimp,'gauss with shifted mean',ion,this['NR_centroid'],fiducial_sigma)

  this['pywimp_pdf'] = wm.get_WIMP_model_with_escape_vel(mass_of_wimp)
  this['wimp_spectrum'] = this['pywimp_pdf']
  result['f_norm'] = wm.get_normalization().getVal()
  result['N_wimp'] = this['pywimp_pdf'].expectedEvents(RooArgSet(time,rec)) # number of expected events per nucleus cross section [pb] for detector mass in time and energy range

  # signal
  this['final_wimp_pdf'] = RooProdPdf('final_wimp_pdf_%sGeV'%mass_of_wimp,'final_wimp_pdf',this['wimp_spectrum'],this['ion_gauss_NR'])

  # efficiency correction
  this['wimp_hist'] = this['final_wimp_pdf'].createHistogram('wimp_hist_%sGeV'%mass_of_wimp,rec,RooFit.Binning(Energy['rec']['bins']),RooFit.YVar(ion,RooFit.Binning(Energy['ion']['bins'])))
  this['wimp_hist'].Scale(result['N_wimp']/this['wimp_hist'].Integral('WIDTH'))

  this['signal_hist'] = this['wimp_hist'].Clone('signal_hist_%sGeV'%mass_of_wimp);this['signal_hist'].Multiply(total_efficiency)

  #additional gamma cut correction
  this['signal_hist'].Multiply(gamma_cut_efficiency)

  result['N_signal'] = this['signal_hist'].Integral('WIDTH')


  # efficiency corrected histograms back to pdf
  this['signal_datahist'] = RooDataHist('signal_datahist_%sGeV'%mass_of_wimp,'signal_datahist',RooArgList(rec,ion),this['signal_hist'])
  this['signal_pdf'] = RooHistPdf('signal_pdf_%sGeV'%mass_of_wimp,'signal_pdf',RooArgSet(rec,ion),this['signal_datahist'])

  result['cross section limit'] = result['f_norm'] * (2.35/result['N_signal'])
  return True


DetectorName = 'ID3'
ID = Detector(DetectorName)
total_efficiency = ID.GetProjectedEfficiency('total')
gamma_cut_efficiency = ID.GetGammaCutEfficiency()


# list of WIMP masses
WIMP_Masses = [5,6,7,8,9,10,12,15,20,25,30]

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
voltage = RooRealVar('voltage','applied voltage',ID.GetWeightedAverage('voltage'))
voltage.setConstant()

heat_baseline_nr = GetEnergyRecoilFromEstimator(ID.GetWeightedAverage('heat'),voltage.getVal())
heat_FWHM = RooRealVar('heat_FWHM','FWHM heat channel baseline',heat_baseline_nr)
heat_sigma = RooFormulaVar('heat_sigma','heat channel baseline resolution','@0/sqrt(8*log(2))',RooArgList(heat_FWHM))

fiducial_FWHM = RooRealVar('fiducial_FWHM','FWHM fiducial electrode baseline',ID.GetWeightedAverage('fiducialmean'))
fiducial_sigma = RooFormulaVar('fiducial_sigma','fiducial electrode baseline resolution','@0/sqrt(8*log(2))', RooArgList(fiducial_FWHM))

# dataset
realdata = RooDataSet.read('ID3-eventlist_test.txt',RooArgList(time,rec,ion))
realdata_scatter = realdata.createHistogram(rec,ion)
events = int(realdata.numEntries())

# loop over several wimp masses
for mass_of_wimp in WIMP_Masses:
  Perform(mass_of_wimp)
