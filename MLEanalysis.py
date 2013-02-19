#!/usr/bin/env python
from ROOT import *
from Parameters import *
from DetectorClass import *
import pyWIMP.DMModels.wimp_model as wimp_model
from pyWIMP.DMModels.flat_model import FlatModel
from pyWIMP.DMModels.base_model import BaseVariables


DetectorName = 'ID3'
ID = Detector(DetectorName)
total_efficiency = ID.GetProjectedEfficiency('total')
mass_of_wimp = 10

# pyWIMP models
basevars = BaseVariables(ID.GetTimeBinning()[0],ID.GetTimeBinning()[-1],Energy['rec']['min'],Energy['rec']['max'],True,time_offset=-2./365.25) #offset to start at 1st Mar 2009
wm = wimp_model.WIMPModel(basevars,mass_of_wimp,ID.GetMass(),constant_quenching=True,nucl_recoil=True)
pywimp_pdf = wm.get_WIMP_model_with_escape_vel(mass_of_wimp)


# observables
ion = RooRealVar('ion','E_{ion}',Energy['ion']['min'],Energy['ion']['max'],'keV_{ee}')
time = basevars.get_time()
rec = basevars.get_energy()
rec.SetNameTitle('rec','E_{rec}')
rec.setUnit('keV_{nr}')

f_norm = wm.get_normalization().getVal()
n_wimp = pywimp_pdf.expectedEvents(RooArgSet(time,rec)) # number of expected events per nucleus cross section [pb] for detector mass in time and energy range


# detector parameters
voltage = RooRealVar('voltage','applied voltage',ID.GetWeightedAverage('voltage'))
voltage.setConstant()

heat_baseline_nr = GetEnergyRecoilFromEstimator(ID.GetWeightedAverage('heat'),voltage.getVal())
heat_FWHM = RooRealVar('heat_FWHM','FWHM heat channel baseline',heat_baseline_nr)
heat_sigma = RooFormulaVar('heat_sigma','heat channel baseline resolution','@0/sqrt(8*log(2))',RooArgList(heat_FWHM))

fiducial_FWHM = RooRealVar('fiducial_FWHM','FWHM fiducial electrode baseline',ID.GetWeightedAverage('fiducialmean'))
fiducial_sigma = RooFormulaVar('fiducial_sigma','fiducial electrode baseline resolution','@0/sqrt(8*log(2))', RooArgList(fiducial_FWHM))


# fit parameters
ratio = RooRealVar('ratio','ratio',0.5,0,1)

# position of bands
# for root plotting
ER_band = TF1('ER_centroid_func','x*(1+(0.16*x^0.18)*([0]/3))/(1+[0]/3)',Energy['rec']['min'],Energy['rec']['max'])
ER_band.SetParName(0,'voltage')
ER_band.SetParameter(0,voltage.getVal())
NR_band = TF1('NR_centroid_func','0.16*x^1.18',Energy['rec']['min'],Energy['rec']['max'])

# for RooFit
ER_centroid = RooFormulaVar('ER_centroid','@0*(1+(0.16*@0^0.18)*(@1/3))/(1+@1/3)',RooArgList(rec,voltage))
NR_centroid = RooFormulaVar('NR_centroid','0.16*@0^1.18',RooArgList(rec))


# gaussians in ionization energy with shifting mean in recoil energy
ion_gauss_ER = RooGaussian('ion_gauss_ER','Gaussian in E_ion with shifted mean in E_rec for ER band',ion,ER_centroid,fiducial_sigma)
ion_gauss_NR = RooGaussian('ion_gauss_NR','Gaussian in E_ion with shifted mean in E_rec for NR band',ion,NR_centroid,fiducial_sigma)

gamma_linear = RooPolynomial('gamma_spectrum','gamma_spectrum',ion)
gamma_spectrum = gamma_linear

wimp_spectrum = pywimp_pdf

## alternative expontial wimp spectrum
#wimp_slope = RooRealVar('wimp_slope','wimp_slope',-0.2)
#wimp_slope.setConstant()
#wimp_exponential = RooExponential('wimp_exponential','wimp_exponential',rec,wimp_slope)
#wimp_spectrum = wimp_exponential


# heat convolution
heat_mean = RooRealVar('heat_mean','heat_mean',0)
heat_smearing = RooGaussian('heat_smearing','heat_smearing',rec,heat_mean,heat_sigma)
rec.setBins(10000,'fft')

# main pdfs
gamma_pdf = RooProdPdf('gamma_pdf','gamma_pdf',ion_gauss_ER,gamma_spectrum)
wimp_pdf = RooProdPdf('wimp_pdf','wimp_pdf',ion_gauss_NR,wimp_spectrum)

gamma_pdf_smeared = RooFFTConvPdf('gamma_pdf_smeared','gamma pdf * heat gaussian',rec,gamma_pdf,heat_smearing)
wimp_pdf_smeared = RooFFTConvPdf('wimp_pdf_smeared','wimp pdf * heat gaussian',rec,wimp_pdf,heat_smearing)

final_gamma_pdf = gamma_pdf_smeared
#final_gamma_pdf = gamma_pdf

#final_wimp_pdf = wimp_pdf_smeared
final_wimp_pdf = wimp_pdf

gamma_hist = final_gamma_pdf.createHistogram('gamma_hist',rec,RooFit.Binning(Energy['rec']['bins']),RooFit.YVar(ion,RooFit.Binning(Energy['ion']['bins'])))
wimp_hist = final_wimp_pdf.createHistogram('wimp_hist',rec,RooFit.Binning(Energy['rec']['bins']),RooFit.YVar(ion,RooFit.Binning(Energy['ion']['bins'])))

bckgd_hist = gamma_hist.Clone('bckgd_hist');bckgd_hist.Multiply(total_efficiency)
signal_hist = wimp_hist.Clone('signal_hist');signal_hist.Multiply(total_efficiency)


# efficiency corrected histograms back to pdf
bckgd_datahist = RooDataHist('bckgd_datahist','bckgd_datahist',RooArgList(rec,ion),bckgd_hist)
bckgd_pdf = RooHistPdf('bckgd_pdf','bckgd_pdf',RooArgSet(rec,ion),bckgd_datahist)
signal_datahist = RooDataHist('signal_datahist','signal_datahist',RooArgList(rec,ion),signal_hist)
signal_pdf = RooHistPdf('signal_pdf','signal_pdf',RooArgSet(rec,ion),signal_datahist)


# dataset
data = RooDataSet.read('ID3-eventlist.txt',RooArgList(time,rec,ion))
scatter = data.createHistogram(rec,ion)
valid_events = int(data.numEntries())


# definition of final PDF model to fit the data to
#final_pdf = RooAddPdf('final_pdf','r*signal+(1-r)*bckgd',RooArgList(signal_pdf,bckgd_pdf),RooArgList(ratio))
final_pdf = bckgd_pdf
#final_pdf = signal_pdf


## two different possibillities to fit data
#final_pdf.fitTo(data)

#nll = RooNLLVar('nll','nll',final_pdf,data, RooFit.PrintEvalErrors(2))
#minuit = RooMinuit(nll)
#fitresults = minuit.fit('hvr')


# toy data
toydata = final_pdf.generate(RooArgSet(rec,ion),valid_events)
toydata_scatter = toydata.createHistogram(rec,ion)


final_hist = final_pdf.createHistogram('final_hist',rec,RooFit.Binning(Energy['rec']['bins']), RooFit.YVar(ion,RooFit.Binning(Energy['ion']['bins'])))
final_hist.SetTitle('\gamma-bckgd PDF only')

# correct choice of binning
recbins = 1*(Energy['rec']['max']-Energy['rec']['min'])
ionbins = 2*(Energy['ion']['max']-Energy['ion']['min'])

recframe = rec.frame()
data.plotOn(recframe,RooFit.Binning(recbins))
final_pdf.plotOn(recframe,RooFit.LineColor(kBlue))
#signal_pdf.plotOn(recframe,RooFit.LineColor(kRed))
#bckgd_pdf.plotOn(recframe,RooFit.LineColor(kGreen))
toydata.plotOn(recframe,RooFit.Binning(recbins),RooFit.MarkerColor(15),RooFit.LineColor(15))
data.plotOn(recframe,RooFit.Binning(recbins))


ionframe = ion.frame()
data.plotOn(ionframe, RooFit.Binning(ionbins))
final_pdf.plotOn(ionframe,RooFit.LineColor(kBlue))
#signal_pdf.plotOn(ionframe,RooFit.LineColor(kRed))
#bckgd_pdf.plotOn(ionframe,RooFit.LineColor(kGreen))
toydata.plotOn(ionframe,RooFit.Binning(ionbins),RooFit.MarkerColor(15),RooFit.LineColor(15))
data.plotOn(ionframe,RooFit.Binning(ionbins))


c = TCanvas()
c.Divide(2,2)
c.cd(1)
final_hist.Draw('COLZ')
scatter.Draw('SAMES')
scatter.SetMarkerColor(kBlack)
toydata_scatter.Draw('SAMES')
toydata_scatter.SetMarkerColor(kWhite)
ER_centroid_func.Draw('SAMES')
ER_centroid_func.SetLineColor(kBlack)
ER_centroid_func.SetLineWidth(2)
NR_centroid_func.Draw('SAMES')
NR_centroid_func.SetLineColor(kBlack)
NR_centroid_func.SetLineWidth(2)
c.cd(2)
ionframe.Draw()
c.cd(3)
recframe.Draw()
c.cd(4)
final_hist.Draw('LEGO2')