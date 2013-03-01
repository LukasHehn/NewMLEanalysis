#!/usr/bin/env python
from ROOT import *
from Parameters import *
from DetectorClass import *
import pyWIMP.DMModels.wimp_model as wimp_model
from pyWIMP.DMModels.flat_model import FlatModel
from pyWIMP.DMModels.base_model import BaseVariables

def GetGammaCutEfficiency():
  hist = TH2F('gamma_cut_efficiency', 'Gamma Cut Efficiency;E_{Rec} (keV_{nr});Efficiency', \
                                      ID.fEnergyRecBinning.size-1, ID.fEnergyRecBinning.flatten('C'), ID.fEnergyIonBinning.size-1, ID.fEnergyIonBinning.flatten('C'))
  ER_centroid.SetParameter(0,6.4)
  for xbin in range(1,hist.GetNbinsX()+1):
    Erec = hist.GetXaxis().GetBinCenter(xbin)
    ioncut = ER_centroid.Eval(Erec)-0.9
    for ybin in range(1,hist.GetNbinsY()+1):
      Eion = hist.GetYaxis().GetBinCenter(ybin)
      if Eion < ioncut:
        hist.SetBinContent(xbin, ybin, 1)
      else:
        hist.SetBinContent(xbin, ybin, 0)
  return hist


def PerformMLEfit(mass_of_wimp,MC_Simulations):
  Results[mass_of_wimp] = {}
  result = Results[mass_of_wimp]

  Container[mass_of_wimp] = {}
  this = Container[mass_of_wimp]

  # fit parameter
  this['ratio'] = RooRealVar('ratio_%sGeV'%mass_of_wimp,'ratio_%sGeV'%mass_of_wimp,ratio_range['start'],ratio_range['min'],ratio_range['max'])

  # position of bands
  this['ER_centroid'] = RooFormulaVar('ER_centroid_%sGeV'%mass_of_wimp,'@0*(1+(0.16*@0^0.18)*(@1/3))/(1+@1/3)',RooArgList(rec,voltage))
  this['NR_centroid'] = RooFormulaVar('NR_centroid_%sGeV'%mass_of_wimp,'0.16*@0^1.18',RooArgList(rec))

  # gaussians in ionization energy with shifting mean in recoil energy
  this['ion_gauss_ER'] = RooGaussian('ion_gauss_ER_%sGeV'%mass_of_wimp,'gauss with shifted mean',ion,this['ER_centroid'],fiducial_sigma)
  this['ion_gauss_NR'] = RooGaussian('ion_gauss_NR_%sGeV'%mass_of_wimp,'gauss with shifted mean',ion,this['NR_centroid'],fiducial_sigma)

  this['gamma_linear'] = RooPolynomial('gamma_linear_%sGeV'%mass_of_wimp,'linear gamma spectrum',ion)
  this['gamma_spectrum'] = this['gamma_linear']

  this['pywimp_pdf'] = wm.get_WIMP_model_with_escape_vel(mass_of_wimp)
  this['wimp_spectrum'] = this['pywimp_pdf']
  result['f_norm'] = wm.get_normalization().getVal()
  result['N_wimp'] = this['pywimp_pdf'].expectedEvents(RooArgSet(time,rec)) # number of expected events per nucleus cross section [pb] for detector mass in time and energy range

  # heat convolution
  this['heat_mean'] = RooRealVar('heat_mean_%sGeV'%mass_of_wimp,'heat_mean',0)
  this['heat_smearing'] = RooGaussian('heat_smearing_%sGeV'%mass_of_wimp,'heat_smearing',rec,this['heat_mean'],heat_sigma)
  rec.setBins(10000,'fft')


  # background
  this['gamma_pdf'] = RooProdPdf('gamma_pdf_%sGeV'%mass_of_wimp,'gamma_pdf',this['gamma_spectrum'],this['ion_gauss_ER'])
  this['gamma_pdf_smeared'] = RooFFTConvPdf('gamma_pdf_smeared_%sGeV'%mass_of_wimp,'gamma pdf * heat gaussian',rec,this['gamma_pdf'],this['heat_smearing'])
  this['final_gamma_pdf'] = this['gamma_pdf_smeared']

  # signal
  #this['wimp_spectrum_smeared'] = RooFFTConvPdf('wimp_spectrum_smeared_%sGeV'%mass_of_wimp,'wimp spectrum * heat gaussian',rec,this['wimp_spectrum'],this['heat_smearing'],2)
  #this['wimp_spectrum_smeared'].setBufferFraction(1.0)
  #this['final_wimp_pdf'] = RooProdPdf('final_wimp_pdf_%sGeV'%mass_of_wimp,'final_wimp_pdf',this['wimp_spectrum_smeared'],this['ion_gauss_NR'])
  this['final_wimp_pdf'] = RooProdPdf('final_wimp_pdf_%sGeV'%mass_of_wimp,'final_wimp_pdf',this['wimp_spectrum'],this['ion_gauss_NR'])


  # efficiency correction
  this['gamma_hist'] = this['final_gamma_pdf'].createHistogram('gamma_hist_%sGeV'%mass_of_wimp,rec,RooFit.Binning(Energy['rec']['bins']),RooFit.YVar(ion,RooFit.Binning(Energy['ion']['bins'])))

  this['wimp_hist'] = this['final_wimp_pdf'].createHistogram('wimp_hist_%sGeV'%mass_of_wimp,rec,RooFit.Binning(Energy['rec']['bins']),RooFit.YVar(ion,RooFit.Binning(Energy['ion']['bins'])))
  this['wimp_hist'].Scale(result['N_wimp']/this['wimp_hist'].Integral('WIDTH'))

  this['bckgd_hist'] = this['gamma_hist'].Clone('bckgd_hist_%sGeV'%mass_of_wimp);this['bckgd_hist'].Multiply(total_efficiency)
  this['signal_hist'] = this['wimp_hist'].Clone('signal_hist_%sGeV'%mass_of_wimp);this['signal_hist'].Multiply(total_efficiency)

  #additional gamma cut correction
  this['signal_hist'].Multiply(gamma_cut_efficiency)

  result['N_signal'] = this['signal_hist'].Integral('WIDTH')


  # efficiency corrected histograms back to pdf
  this['bckgd_datahist'] = RooDataHist('bckgd_datahist_%sGeV'%mass_of_wimp,'bckgd_datahist',RooArgList(rec,ion),this['bckgd_hist'])
  this['bckgd_pdf'] = RooHistPdf('bckgd_pdf_%sGeV'%mass_of_wimp,'bckgd_pdf',RooArgSet(rec,ion),this['bckgd_datahist'])
  this['signal_datahist'] = RooDataHist('signal_datahist_%sGeV'%mass_of_wimp,'signal_datahist',RooArgList(rec,ion),this['signal_hist'])
  this['signal_pdf'] = RooHistPdf('signal_pdf_%sGeV'%mass_of_wimp,'signal_pdf',RooArgSet(rec,ion),this['signal_datahist'])

  # definition of final PDF model to fit the realdata_scatter to
  this['final_pdf'] = RooAddPdf('final_pdf_%sGeV'%mass_of_wimp,'r*signal+(1-r)*bckgd',RooArgList(this['signal_pdf'],this['bckgd_pdf']),RooArgList(this['ratio']))


  # maximum likelihood fit to data
  this['nll'] = RooNLLVar('nll_%sGeV'%mass_of_wimp,'nll',this['final_pdf'],realdata,RooFit.PrintEvalErrors(2))
  this['minuit'] = RooMinuit(this['nll'])
  this['fitresults'] = this['minuit'].fit('hvr')

  this['final_hist'] = this['final_pdf'].createHistogram('final_hist_%sGeV'%mass_of_wimp,rec,RooFit.Binning(Energy['rec']['bins']), RooFit.YVar(ion,RooFit.Binning(Energy['ion']['bins'])))
  this['final_hist'].SetTitle('fitted signal + bckgd PDF for % GeV'%mass_of_wimp)

  if MC_Simulations != 0:
    # testing the fit parameter distribution for toy data
    this['MC_study'] = RooMCStudy(this['final_pdf'],RooArgSet(rec,ion))
    this['MC_study'].generateAndFit(MC_Simulations,events,kTRUE)


  result['ratio'] = this['ratio'].getVal()
  result['ratio error high'] = this['ratio'].getErrorHi()
  result['upper ratio 90'] = result['ratio'] + 1.64*result['ratio error high']
  result['event limit'] = result['upper ratio 90'] * events
  #result['cross section limit'] = result['f_norm'] * (result['event limit']/result['N_signal'])
  result['cross section limit'] = result['f_norm'] * (2.35/result['N_signal'])
  return True


def PlotFitResults(mass_of_wimp,toydata):
  this = Container[mass_of_wimp]

  if toydata == True:
    # toy data generation
    this['toyevents'] = this['final_pdf'].generate(RooArgSet(rec,ion),events)
    this['toyevents_scatter'] = this['toyevents'].createHistogram(rec,ion)

  # correct choice of binning
  recbins = 1*(Energy['rec']['max']-Energy['rec']['min'])
  ionbins = 2*(Energy['ion']['max']-Energy['ion']['min'])

  # RooFit frames
  recframe = rec.frame()
  realdata.plotOn(recframe,RooFit.Binning(recbins))
  this['final_pdf'].plotOn(recframe,RooFit.LineColor(kBlue))
  this['signal_pdf'].plotOn(recframe,RooFit.LineColor(kRed))
  this['bckgd_pdf'].plotOn(recframe,RooFit.LineColor(kGreen))
  if toydata == True:
    this['toyevents'].plotOn(recframe,RooFit.Binning(recbins),RooFit.MarkerColor(15),RooFit.LineColor(17))
  realdata.plotOn(recframe,RooFit.Binning(recbins))

  ionframe = ion.frame()
  realdata.plotOn(ionframe, RooFit.Binning(ionbins))
  this['final_pdf'].plotOn(ionframe,RooFit.LineColor(kBlue))
  this['signal_pdf'].plotOn(ionframe,RooFit.LineColor(kRed))
  this['bckgd_pdf'].plotOn(ionframe,RooFit.LineColor(kGreen))
  if toydata == True:
    this['toyevents'].plotOn(ionframe,RooFit.Binning(ionbins),RooFit.MarkerColor(15),RooFit.LineColor(17))
  realdata.plotOn(ionframe,RooFit.Binning(ionbins))

  ratiobins = (ratio_range['max']-ratio_range['min'])*10
  #paramframe = this['MC_study'].plotParam(this['ratio'],RooFit.FrameRange(ratio_range['min'],ratio_range['max']),RooFit.Binning(150))

  #nllframe = this['MC_study'].plotNLL()

  c0 = TCanvas('c0','MLE fit results and MC toy fit for %s GeV' % mass_of_wimp,800,600)
  c0.Divide(2,2)
  c0.cd(1)
  this['final_hist'].Draw('COLZ')
  realdata_scatter.Draw('SAMES')
  realdata_scatter.SetMarkerColor(kBlack)
  if toydata == True:
    this['toydata_scatter'].Draw('SAMES')
    this['toydata_scatter'].SetMarkerColor(kWhite)
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
  #sub = c0.cd(4)
  #sub.Divide(1,2)
  ## root line
  #line = TLine()
  #line.SetLineWidth(2)
  #line.SetLineColor(kRed)
  #sub.cd(1)
  #paramframe.Draw()
  #c0.Update()
  ## ratio marker line for real data
  #line.SetLineStyle(0)
  #line.DrawLine(Results[mass_of_wimp]['ratio'],gPad.GetUymin(),Results[mass_of_wimp]['ratio'],gPad.GetUymax())
  #line.SetLineStyle(7)
  #line.DrawLine(Results[mass_of_wimp]['ratio']+Results[mass_of_wimp]['ratio error high'],gPad.GetUymin(),Results[mass_of_wimp]['ratio']+Results[mass_of_wimp]['ratio error high'],gPad.GetUymax())
  #sub.cd(2)
  #nllframe.Draw()
  #c0.Update()
  ## nll marker line for real data
  #line.SetLineStyle(0)
  #line.DrawLine(this['nll'].getVal(),gPad.GetUymin(),this['nll'].getVal(),gPad.GetUymax())
  return c0


def PlotExclusionCurve():
  Exclusion['exclusion curve Eric'] = TGraph(len(EricsResults[DetectorName]))
  Exclusion['exclusion curve MLE'] = TGraph(len(Results))
  i = 0
  for mass_of_wimp in EricsResults[DetectorName]:
    Exclusion['exclusion curve Eric'].SetPoint(i, mass_of_wimp, EricsResults[DetectorName][mass_of_wimp])
    i += 1
  j = 0
  for mass_of_wimp in Results:
    Exclusion['exclusion curve MLE'].SetPoint(j, mass_of_wimp, Results[mass_of_wimp]['cross section limit'])
    j += 1
  c0 = TCanvas('canvas','Exclusion Plot',800,600)
  #gPad.SetLogy()
  Exclusion['exclusion curve Eric'].SetTitle('Exclusion Plot %s GeV' % DetectorName)
  Exclusion['exclusion curve Eric'].Draw('ALP')
  Exclusion['exclusion curve Eric'].GetYaxis().SetTitle('\sigma_{\chi - nucleon} [pb]')
  Exclusion['exclusion curve Eric'].GetXaxis().SetTitle('m_{\chi } [GeV]')
  Exclusion['exclusion curve Eric'].SetLineColor(kBlue)
  Exclusion['exclusion curve Eric'].SetLineWidth(3)
  Exclusion['exclusion curve MLE'].Draw('LP')
  Exclusion['exclusion curve MLE'].SetLineColor(kRed)
  Exclusion['exclusion curve MLE'].SetLineWidth(3)
  return c0


DetectorName = 'ID3'
ID = Detector(DetectorName)
total_efficiency = ID.GetProjectedEfficiency('total')
gamma_cut_efficiency = GetGammaCutEfficiency()

ratio_range = {'min' : 0, 'max' : 1.0, 'start' : 0.0}

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


MC_Simulations = 0

# loop over several wimp masses
print 'fit results for {0} with {1} kg of mass'.format(DetectorName,ID.GetMass())
print "-----------------------------------------------------------------"
print '{0:8} | {1:8} | {2:8} | {3:8} | {4:8} | {5:8} | {6:8}'.format('m_chi','r','r_sigma','N_90','N_pywimp','N_signal','xs')
for mass_of_wimp in WIMP_Masses:
  PerformMLEfit(mass_of_wimp,MC_Simulations)
  this = Container[mass_of_wimp]
  result = Results[mass_of_wimp]
  print '{0:8} | {1:8.2g} | {2:8.2g} | {3:8.2g} | {4:8.2g} | {5:8.2g} | {6:8.1e}'.format(mass_of_wimp,result['ratio'],result['ratio error high'],result['event limit'],result['N_wimp'],result['N_signal'],result['cross section limit'])
