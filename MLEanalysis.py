#!/usr/bin/env python
from ROOT import *
from Parameters import *
from DetectorClass import *
import pyWIMP.DMModels.wimp_model as wimp_model
from pyWIMP.DMModels.flat_model import FlatModel
from pyWIMP.DMModels.base_model import BaseVariables


def PerformMLE(mass_of_wimp, MC_Simulations):
  Container[mass_of_wimp] = {}
  c = Container[mass_of_wimp]

  # fit parameter
  c['ratio'] = RooRealVar('ratio_%GeV'%mass_of_wimp,'ratio_%GeV'%mass_of_wimp,0.5,0,1)

  # position of bands
  c['ER_centroid'] = RooFormulaVar('ER_centroid_%GeV'%mass_of_wimp,'@0*(1+(0.16*@0^0.18)*(@1/3))/(1+@1/3)',RooArgList(rec,voltage))
  c['NR_centroid'] = RooFormulaVar('NR_centroid_%GeV'%mass_of_wimp,'0.16*@0^1.18',RooArgList(rec))

  # gaussians in ionization energy with shifting mean in recoil energy
  c['ion_gauss_ER'] = RooGaussian('ion_gauss_ER_%GeV'%mass_of_wimp,'Gaussian in E_ion with shifted mean in E_rec for ER band',ion,c['ER_centroid'],fiducial_sigma)
  c['ion_gauss_NR'] = RooGaussian('ion_gauss_NR_%GeV'%mass_of_wimp,'Gaussian in E_ion with shifted mean in E_rec for NR band',ion,c['NR_centroid'],fiducial_sigma)

  c['gamma_linear'] = RooPolynomial('gamma_linear_%GeV'%mass_of_wimp,'linear gamma spectrum',ion)
  c['gamma_spectrum'] = c['gamma_linear']

  c['pywimp_pdf'] = wm.get_WIMP_model_with_escape_vel(mass_of_wimp)
  c['wimp_spectrum'] = c['pywimp_pdf']
  c['f_norm'] = wm.get_normalization().getVal()
  c['N_wimp'] = c['pywimp_pdf'].expectedEvents(RooArgSet(time,rec)) # number of expected events per nucleus cross section [pb] for detector mass in time and energy range

  # heat convolution
  c['heat_mean'] = RooRealVar('heat_mean_%GeV'%mass_of_wimp,'heat_mean',0)
  c['heat_smearing'] = RooGaussian('heat_smearing_%GeV'%mass_of_wimp,'heat_smearing',rec,c['heat_mean'],heat_sigma)
  rec.setBins(1000,'fft')

  # background
  c['gamma_pdf'] = RooProdPdf('gamma_pdf_%GeV'%mass_of_wimp,'gamma_pdf',c['ion_gauss_ER'],c['gamma_spectrum'])
  c['gamma_pdf_smeared'] = RooFFTConvPdf('gamma_pdf_smeared_%GeV'%mass_of_wimp,'gamma pdf * heat gaussian',rec,c['gamma_pdf'],c['heat_smearing'],2)
  c['final_gamma_pdf'] = c['gamma_pdf_smeared']

  # signal
  c['wimp_spectrum_smeared'] = RooFFTConvPdf('wimp_spectrum_smeared_%GeV'%mass_of_wimp,'wimp spectrum * heat gaussian',rec,c['wimp_spectrum'],c['heat_smearing'],2)
  c['wimp_spectrum_smeared'].setBufferFraction(0.7)
  c['final_wimp_pdf'] = RooProdPdf('final_wimp_pdf_%GeV'%mass_of_wimp,'final_wimp_pdf',c['ion_gauss_NR'],c['wimp_spectrum_smeared'])


  # efficiency correction
  c['gamma_hist'] = c['final_gamma_pdf'].createHistogram('gamma_hist_%GeV'%mass_of_wimp,rec,RooFit.Binning(Energy['rec']['bins']),RooFit.YVar(ion,RooFit.Binning(Energy['ion']['bins'])))

  c['wimp_hist'] = c['final_wimp_pdf'].createHistogram('wimp_hist_%GeV'%mass_of_wimp,rec,RooFit.Binning(Energy['rec']['bins']),RooFit.YVar(ion,RooFit.Binning(Energy['ion']['bins'])))
  c['wimp_hist'].Scale(c['N_wimp']/c['wimp_hist'].Integral('WIDTH'))

  c['bckgd_hist'] = c['gamma_hist'].Clone('bckgd_hist_%GeV'%mass_of_wimp);c['bckgd_hist'].Multiply(total_efficiency)
  c['signal_hist'] = c['wimp_hist'].Clone('signal_hist_%GeV'%mass_of_wimp);c['signal_hist'].Multiply(total_efficiency)

  c['N_signal'] = c['signal_hist'].Integral('WIDTH')


  # efficiency corrected histograms back to pdf
  c['bckgd_datahist'] = RooDataHist('bckgd_datahist_%GeV'%mass_of_wimp,'bckgd_datahist',RooArgList(rec,ion),c['bckgd_hist'])
  c['bckgd_pdf'] = RooHistPdf('bckgd_pdf_%GeV'%mass_of_wimp,'bckgd_pdf',RooArgSet(rec,ion),c['bckgd_datahist'])
  c['signal_datahist'] = RooDataHist('signal_datahist_%GeV'%mass_of_wimp,'signal_datahist',RooArgList(rec,ion),c['signal_hist'])
  c['signal_pdf'] = RooHistPdf('signal_pdf_%GeV'%mass_of_wimp,'signal_pdf',RooArgSet(rec,ion),c['signal_datahist'])

  # definition of final PDF model to fit the data to
  c['final_pdf'] = RooAddPdf('final_pdf_%GeV'%mass_of_wimp,'r*signal+(1-r)*bckgd',RooArgList(c['signal_pdf'],c['bckgd_pdf']),RooArgList(c['ratio']))


  # maximum likelihood fit to data
  c['nll'] = RooNLLVar('nll_%GeV'%mass_of_wimp,'nll',c['final_pdf'],data,RooFit.PrintEvalErrors(2))
  c['minuit'] = RooMinuit(c['nll'])
  c['fitresults'] = c['minuit'].fit('hvr')

  c['final_hist'] = c['final_pdf'].createHistogram('final_hist_%GeV'%mass_of_wimp,rec,RooFit.Binning(Energy['rec']['bins']), RooFit.YVar(ion,RooFit.Binning(Energy['ion']['bins'])))
  c['final_hist'].SetTitle('fitted signal + bckgd PDF for % GeV'%mass_of_wimp)

  if MC_Simulations != 0:
    # testing the fit parameter distribution for toy data
    c['MC_study'] = RooMCStudy(c['final_pdf'],RooArgSet(rec,ion))
    c['MC_study'].generateAndFit(MC_Simulations,events,kTRUE)

  Results[mass_of_wimp] = {}
  result = Results[mass_of_wimp]
  result['ratio'] = c['ratio'].getVal()
  result['ratio error high'] = c['ratio'].getErrorHi()
  result['ratio error 90'] = result['ratio'] + 1.64*result['ratio error high']
  result['event limit'] = result['ratio error 90'] * events
  result['cross section limit'] = c['f_norm'] * (result['event limit']/c['N_signal'])
  return True


def PlotFitResults(mass_of_wimp):
  ## toy data generation
  #toydata = final_pdf.generate(RooArgSet(rec,ion),events)
  #toydata_scatter = toydata.createHistogram(rec,ion)

  # correct choice of binning
  recbins = 1*(Energy['rec']['max']-Energy['rec']['min'])
  ionbins = 2*(Energy['ion']['max']-Energy['ion']['min'])

  # RooFit frames
  recframe = rec.frame()
  data.plotOn(recframe,RooFit.Binning(recbins))
  Container[mass_of_wimp]['final_pdf'].plotOn(recframe,RooFit.LineColor(kBlue))
  Container[mass_of_wimp]['signal_pdf'].plotOn(recframe,RooFit.LineColor(kRed))
  Container[mass_of_wimp]['bckgd_pdf'].plotOn(recframe,RooFit.LineColor(kGreen))
  #toydata.plotOn(recframe,RooFit.Binning(recbins),RooFit.MarkerColor(15),RooFit.LineColor(17))
  data.plotOn(recframe,RooFit.Binning(recbins))

  ionframe = ion.frame()
  data.plotOn(ionframe, RooFit.Binning(ionbins))
  Container[mass_of_wimp]['final_pdf'].plotOn(ionframe,RooFit.LineColor(kBlue))
  Container[mass_of_wimp]['signal_pdf'].plotOn(ionframe,RooFit.LineColor(kRed))
  Container[mass_of_wimp]['bckgd_pdf'].plotOn(ionframe,RooFit.LineColor(kGreen))
  #toydata.plotOn(ionframe,RooFit.Binning(ionbins),RooFit.MarkerColor(15),RooFit.LineColor(17))
  data.plotOn(ionframe,RooFit.Binning(ionbins))

  paramframe = Container[mass_of_wimp]['MC_study'].plotParam(Container[mass_of_wimp]['ratio'],RooFit.FrameRange(0.,0.2),RooFit.Binning(200))

  nllframe = Container[mass_of_wimp]['MC_study'].plotNLL()

  c0 = TCanvas('c0','MLE fit results and MC toy fit for %s GeV' % mass_of_wimp,800,600)
  c0.Divide(2,2)
  c0.cd(1)
  Container[mass_of_wimp]['final_hist'].Draw('COLZ')
  scatter.Draw('SAMES')
  scatter.SetMarkerColor(kBlack)
  #toydata_scatter.Draw('SAMES')
  #toydata_scatter.SetMarkerColor(kWhite)
  ER_centroid.SetParameter(0,voltage.getVal())
  ER_centroid.Draw('SAMES')
  ER_centroid.SetLineColor(kBlack)
  ER_centroid.SetLineWidth(2)
  NR_centroid.Draw('SAMES')
  NR_centroid.SetLineColor(kBlack)
  NR_centroid.SetLineWidth(2)
  c0.cd(2)
  ionframe.Draw()
  c0.cd(3)
  recframe.Draw()
  c0_4.Divide(1,2)
  # root line
  line = TLine()
  line.SetLineWidth(2)
  line.SetLineColor(kRed)
  c0_4.cd(1)
  paramframe.Draw()
  c0.Update()
  # ratio marker line for real data
  line.SetLineStyle(0)
  line.DrawLine(Results[mass_of_wimp]['ratio'],gPad.GetUymin(),Results[mass_of_wimp]['ratio'],gPad.GetUymax())
  line.SetLineStyle(7)
  line.DrawLine(Results[mass_of_wimp]['ratio']+Results[mass_of_wimp]['ratio error high'],gPad.GetUymin(),Results[mass_of_wimp]['ratio']+Results[mass_of_wimp]['ratio error high'],gPad.GetUymax())
  c0_4.cd(2)
  nllframe.Draw()
  c0.Update()
  # nll marker line for real data
  line.SetLineStyle(0)
  line.DrawLine(Container[mass_of_wimp]['nll'].getVal(),gPad.GetUymin(),Container[mass_of_wimp]['nll'].getVal(),gPad.GetUymax())
  return c0


DetectorName = 'ID3'
ID = Detector(DetectorName)
total_efficiency = ID.GetProjectedEfficiency('total')

# list of WIMP masses
WIMP_Masses = [8, 10]

# dictionaries
Container = {}
Results = {}

# pyWIMP model
basevars = BaseVariables(ID.GetTimeBinning()[0],ID.GetTimeBinning()[-1],Energy['rec']['min'],Energy['rec']['max'],True,time_offset=-2./365.25) #offset to start at 1st Mar 2009
wm = wimp_model.WIMPModel(basevars,WIMP_Masses[0],ID.GetMass(),constant_quenching=True,nucl_recoil=True)

# observable definition
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
data = RooDataSet.read('ID3-eventlist_test.txt',RooArgList(time,rec,ion))
scatter = data.createHistogram(rec,ion)
events = int(data.numEntries())


# loop over several wimp masses
print "fit results:"
print '{0:8} | {1:8} | {2:8} | {3:8} | {4:8} | {5:8} | {6:8} | {7:8}'.format('m_chi','r','r_sigma','r_90','N_max','N_obs','N_chi','xs')
for mass_of_wimp in WIMP_Masses:
  PerformMLE(mass_of_wimp)
  for mass_of_wimp in WIMP_Masses:
    print '{0:8} | {1:8.2g} | {2:8.2g} | {3:8.2g} | {4:8.2g} | {5:8.2g} | {6:8.2g} | {7:8.1e}'.format(mass_of_wimp,R[i],R_MAX[i],R_90[i],N_OBS[i],N_CHI[i],N_SIG[i],XS[i])
    i += 1
