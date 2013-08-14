#!/usr/bin/env python
from ROOT import *
from Functions import *
from DetectorClass import *


# define maximal energies
EnergyIonMax = 10.
EnergyRecMax = 20.


#switches and input parameters to control script
wimp_mass = 8 #set wimp mass or switch signal of entirely (with False)
MC_sets = int(1e3) #set number of MC simulations: 0 means none at all
SavePlots = True #flag decides whether plots are saved

DataFile = '/kalinka/home/hehn/PhD/LowMassEric/ID3_eventlist.txt'


# definition of observables
ion = RooRealVar('ion','E_{ion}',0.,EnergyIonMax,'keV_{ee}')
rec = RooRealVar('rec','E_{rec}',0.,EnergyRecMax,'keV_{nr}')


# detector specific parameters
voltage = RooRealVar('voltage','applied voltage',6.4)

E_thresh = 3.874
FWHM_heat = 0.82 #Eric@10keV
FWHM_rec = RecoilResolutionFromHeat(FWHM_heat,voltage.getVal(),10)
sigma_rec = RooRealVar('sigma_rec','recoil energy resolution',FWHM_rec/2.35)

FWHM_ion = 0.72 #Eric
sigma_ion = RooRealVar('sigma_ion','ionization energy resolution',FWHM_ion/2.35)


# detector efficiency
total_efficiency = Simple2DEfficiencyID3(E_thresh,FWHM_rec/2.35)
total_efficiency_datahist = RooDataHist('total_efficiency_datahist','total_efficiency_datahist',RooArgList(rec,ion),total_efficiency)
total_efficiency_pdf = RooHistPdf('total_efficiency_pdf','total_efficiency_pdf',RooArgSet(rec,ion),total_efficiency_datahist)


# dataset
realdata = RooDataSet.read(DataFile,RooArgList(rec,ion))
realdata_scatter = realdata.createHistogram(rec,ion,Energy['rec']['bins'],Energy['ion']['bins'])
realdata_graph = TGraphFromDataSet(realdata)
events = int(realdata.numEntries())


# -----------------------------------------------------------------------------------------
# definition of gamma peaks
ion_scaling = RooRealVar('ion_scaling','scaling factor ionization energy',1.0,0.95,1.05)
rec_scaling = RooRealVar('rec_scaling','scaling factor recoil energy',1.0,0.95,1.05)
ion_scaling.setConstant(kTRUE)
rec_scaling.setConstant(kTRUE)
ER_centroid.SetParameter(0,6.4) #ER_centroid used for calculation of peak position in Erec

V49_ion_energy = RooRealVar('V49_ion_energy','V49 peak ion energy',4.97)
V49_ion_pos = RooFormulaVar('V49_ion_pos','@0*@1',RooArgList(V49_ion_energy,ion_scaling))
V49_ion_pdf = RooGaussian('V49_ion_pdf','V49 peak gauss pdf in ion',ion,V49_ion_pos,sigma_ion)
V49_rec_energy = RooRealVar('V49_rec_energy','recoil energy V49 peak',ER_centroid.GetX(V49_ion_energy.getVal()))
V49_rec_pos = RooFormulaVar('V49_rec_pos','@0*@1',RooArgList(V49_rec_energy,rec_scaling))
V49_rec_pdf = RooGaussian('V49_rec_pdf','V49 peak pdf in recoil energy',rec,V49_rec_pos,sigma_rec)
V49_pdf = RooProdPdf('V49_pdf','V49 peak pdf',V49_ion_pdf,V49_rec_pdf)
V49_pdf_eff = RooProdPdf('V49_pdf_eff','eff corr V49 peak pdf',V49_pdf,total_efficiency_pdf)
N_V49 = RooRealVar('N_V49','evts of V49 peak (4.97keV)',16.,0.,30.)


Cr51_ion_energy = RooRealVar('Cr51_ion_energy','Cr51 peak ion energy',5.46)
Cr51_ion_pos = RooFormulaVar('Cr51_ion_pos','@0*@1',RooArgList(Cr51_ion_energy,ion_scaling))
Cr51_ion = RooGaussian('Cr51_ion_pdf','Cr51 peak gauss pdf with shifted mean',ion,Cr51_ion_pos,sigma_ion)
Cr51_rec_energy = RooRealVar('Cr51_rec_energy','v_rec_energy',ER_centroid.GetX(Cr51_ion_energy.getVal()))
Cr51_rec_pos = RooFormulaVar('Cr51_rec_pos','@0*@1',RooArgList(Cr51_rec_energy,rec_scaling))
Cr51_rec = RooGaussian('Cr51_rec_pdf','Cr51_rec_pdf with shifted mean',rec,Cr51_rec_pos,sigma_rec)
Cr51_pdf = RooProdPdf('Cr51_pdf','Cr51 peak pdf',Cr51_ion,Cr51_rec)
N_Cr51 = RooRealVar('N_Cr51','evts of 51Cr peak (5.46keV)',11.,0.,30.)


Mn54_ion_energy = RooRealVar('Mn54_ion_energy','Mn54_ion_energy',5.99)
Mn54_ion_pos = RooFormulaVar('Mn54_ion_pos','@0*@1',RooArgList(Mn54_ion_energy,ion_scaling))
Mn54_ion = RooGaussian('Mn54_ion_pdf','Mn54_ion_pdf with shifted mean',ion,Mn54_ion_pos,sigma_ion)
Mn54_rec_energy = RooRealVar('Mn54_rec_energy','Mn54_rec_energy',ER_centroid.GetX(Mn54_ion_energy.getVal()))
Mn54_rec_pos = RooFormulaVar('Mn54_rec_pos','@0*@1',RooArgList(Mn54_rec_energy,rec_scaling))
Mn54_rec = RooGaussian('Mn54_rec_pdf','Mn54_rec_pdf with shifted mean',rec,Mn54_rec_pos,sigma_rec)
Mn54_pdf = RooProdPdf('Mn54_pdf','Mn54 peak pdf',Mn54_ion,Mn54_rec)
N_Mn54 = RooRealVar('N_Mn54','evts of 54Mn peak (5.99keV)',4.,0.,10.)


Fe55_ion_energy = RooRealVar('Fe55_ion_energy','Fe55_ion_energy',6.54)
Fe55_ion_pos = RooFormulaVar('Fe55_ion_pos','@0*@1',RooArgList(Fe55_ion_energy,ion_scaling))
Fe55_ion = RooGaussian('Fe55_ion_pdf','Fe55_ion_pdf with shifted mean',ion,Fe55_ion_pos,sigma_ion)
Fe55_rec_energy = RooRealVar('Fe55_rec_energy','Fe55_rec_energy',ER_centroid.GetX(Fe55_ion_energy.getVal()))
Fe55_rec_pos = RooFormulaVar('Fe55_rec_pos','@0*@1',RooArgList(Fe55_rec_energy,rec_scaling))
Fe55_rec = RooGaussian('Fe55_rec_pdf','Fe55_rec_pdf with shifted mean',rec,Fe55_rec_pos,sigma_rec)
Fe55_pdf = RooProdPdf('Fe55_pdf','Fe55 peak pdf',Fe55_ion,Fe55_rec)
N_Fe55 = RooRealVar('N_Fe55','evts of 55Fe peak (6.54keV)',31.,0.,60.)


Co57_ion_energy = RooRealVar('Co57_ion_energy','Co57_ion_energy',7.11)
Co57_ion_pos = RooFormulaVar('Co57_ion_pos','@0*@1',RooArgList(Co57_ion_energy,ion_scaling))
Co57_ion = RooGaussian('Co57_ion_pdf','Co57_ion_pdf with shifted mean',ion,Co57_ion_pos,sigma_ion)
Co57_rec_energy = RooRealVar('Co57_rec_energy','Co57_rec_energy',ER_centroid.GetX(Co57_ion_energy.getVal()))
Co57_rec_pos = RooFormulaVar('Co57_rec_pos','@0*@1',RooArgList(Co57_rec_energy,rec_scaling))
Co57_rec = RooGaussian('Co57_rec_pdf','Co57_rec_pdf with shifted mean',rec,Co57_rec_pos,sigma_rec)
Co57_pdf = RooProdPdf('Co57_pdf','Co57 peak pdf',Co57_ion,Co57_rec)
N_Co57 = RooRealVar('N_Co57','evts of 57Co peak (7.11keV)',0.,0.,20.)


Zn65_ion_energy = RooRealVar('Zn65_ion_energy','Zn65_ion_energy',8.98)
Zn65_ion_pos = RooFormulaVar('Zn65_ion_pos','@0*@1',RooArgList(Zn65_ion_energy,ion_scaling))
Zn65_ion = RooGaussian('Zn65_ion_pdf','Zn65_ion_pdf with shifted mean',ion,Zn65_ion_pos,sigma_ion)
Zn65_rec_energy = RooRealVar('Zn65_rec_energy','Zn65_rec_energy',ER_centroid.GetX(Zn65_ion_energy.getVal()))
Zn65_rec_pos = RooFormulaVar('Zn65_rec_pos','@0*@1',RooArgList(Zn65_rec_energy,rec_scaling))
Zn65_rec = RooGaussian('Zn65_rec_pdf','Zn65_rec_pdf with shifted mean',rec,Zn65_rec_pos,sigma_rec)
Zn65_pdf = RooProdPdf('Zn65_pdf','Zn65 peak pdf',Zn65_ion,Zn65_rec)
N_Zn65 = RooRealVar('N_Zn65','evts of 65Zn peak (8.98keV)',110.,70.,130.)


Ga68_ion_energy = RooRealVar('Ga68_ion_energy','Ga68_ion_energy',9.66)
Ga68_ion_pos = RooFormulaVar('Ga68_ion_pos','@0*@1',RooArgList(Ga68_ion_energy,ion_scaling))
Ga68_ion = RooGaussian('Ga68_ion_pdf','Ga68_ion_pdf with shifted mean',ion,Ga68_ion_pos,sigma_ion)
Ga68_rec_energy = RooRealVar('Ga68_rec_energy','Ga68_rec_energy',ER_centroid.GetX(Ga68_ion_energy.getVal()))
Ga68_rec_pos = RooFormulaVar('Ga68_rec_pos','@0*@1',RooArgList(Ga68_rec_energy,rec_scaling))
Ga68_rec = RooGaussian('Ga68_rec_pdf','Ga68_rec_pdf with shifted mean',rec,Ga68_rec_pos,sigma_rec)
Ga68_pdf = RooProdPdf('Ga68_pdf','Ga68 peak pdf',Ga68_ion,Ga68_rec)
N_Ga68 = RooRealVar('N_Ga68','evts of 68Ga peak (9.66keV)',32.,0.,60.)


Ge68_ion_energy = RooRealVar('Ge68_ion_energy','Ge68_ion_energy',10.37)
Ge68_ion_pos = RooFormulaVar('Ge68_ion_pos','@0*@1',RooArgList(Ge68_ion_energy,ion_scaling))
Ge68_ion = RooGaussian('Ge68_ion_pdf','Ge68_ion_pdf with shifted mean',ion,Ge68_ion_pos,sigma_ion)
Ge68_rec_energy = RooRealVar('Ge68_rec_energy','Ge68_rec_energy',ER_centroid.GetX(Ge68_ion_energy.getVal()))
Ge68_rec_pos = RooFormulaVar('Ge68_rec_pos','@0*@1',RooArgList(Ge68_rec_energy,rec_scaling))
Ge68_rec = RooGaussian('Ge68_rec_pdf','Ge68 peak pdf in rec',rec,Ge68_rec_pos,sigma_rec)
Ge68_pdf = RooProdPdf('Ge68_pdf','Ge68 peak pdf',Ge68_ion,Ge68_rec)
N_Ge68 = RooRealVar('N_Ge68','evts of 68Ge peak (10.37keV)',0.,0.,100.)
# -----------------------------------------------------------------------------------------

# wimp signal
if wimp_mass:

  # read in WIMP spectrum
  signal_hist = WimpSignal2DEric(wimp_mass,sigma_ion.getVal(),sigma_rec.getVal())
  signal_hist.Multiply(total_efficiency)
  signal_datahist = RooDataHist('signal_datahist','signal_datahist',RooArgList(rec,ion),signal_hist)
  signal_pdf = RooHistPdf('signal_pdf','signal_pdf',RooArgSet(rec,ion),signal_datahist)
  N_signal = RooRealVar('N_signal','WIMP signal events',0.,0.,10.)


# flat gamma background component
flat_gamma_bckgd_hist = FlatGammaBckgd2DEric(sigma_ion.getVal(),sigma_rec.getVal())
flat_gamma_bckgd_hist.Multiply(total_efficiency)
flat_gamma_bckgd_datahist = RooDataHist('flat_gamma_bckgd_datahist','flat_gamma_bckgd_datahist',RooArgList(rec,ion),flat_gamma_bckgd_hist)
flat_gamma_bckgd_pdf = RooHistPdf('flat_gamma_bckgd_pdf','flat_gamma_bckgd_pdf',RooArgSet(rec,ion),flat_gamma_bckgd_datahist)
N_flat = RooRealVar('N_flat','bckgd events',70.,40.,100.)


if wimp_mass:
  # definition of pdf to be fitted
  final_pdf = RooAddPdf('final_pdf','final_pdf',RooArgList(signal_pdf,flat_gamma_bckgd_pdf, V49_pdf, Cr51_pdf, Mn54_pdf, Fe55_pdf, Co57_pdf, Zn65_pdf, Ga68_pdf),RooArgList(N_signal,N_flat, N_V49, N_Cr51, N_Mn54, N_Fe55, N_Co57, N_Zn65, N_Ga68))
else:
  final_pdf = RooAddPdf('final_pdf','final_pdf',RooArgList(flat_gamma_bckgd_pdf, V49_pdf, Cr51_pdf, Mn54_pdf, Fe55_pdf, Co57_pdf, Zn65_pdf, Ga68_pdf, Ge68_pdf),RooArgList(N_flat, N_V49, N_Cr51, N_Mn54, N_Fe55, N_Co57, N_Zn65, N_Ga68, N_Ge68))


# manual mode
nll = RooNLLVar('nll','nll',final_pdf,realdata,RooFit.NumCPU(2),RooFit.PrintEvalErrors(2),RooFit.Extended(kTRUE),RooFit.Verbose(kFALSE))
minuit = RooMinuit(nll)
minuit.migrad() #find minimum
minuit.hesse() #symmetric errors
minuit.minos() #asymmetric errors
FitResult = minuit.save('realfit','fit to real data')
ndf = FitResult.floatParsFinal().getSize()

FitResult.Print('v')


#calculate and print cross section limit for this wimp mass:
if wimp_mass:
  signal_parameter = FitResult.floatParsFinal().find('N_signal')
  wimp_events = signal_parameter.getVal()
  error_low = signal_parameter.getAsymErrorLo()
  error_high = signal_parameter.getAsymErrorHi()
  wimp_events_limit = wimp_events + 1.28 * error_high
  livetime = 197. #wimp search livetime in days
  mass = 0.160 #detector mass in kg
  rate = signal_hist.Integral('WIDTH')
  signal_events = rate * livetime * mass * 1e6
  sigma_limit = wimp_events_limit / signal_events

  print '{0:9} | {1:10} | {2:8} | {3:8} | {4:10} | {5:10}'.format('wimp_mass','Rate','N_signal','ErrorLow','ErrorHigh','XS-limit [pb]')
  print "-----------------------------------------------------------------------------------"
  print '{0:9} | {1:10} | {2:8} | {3:8} | {4:10} | {5:10}'.format(wimp_mass,rate,wimp_events,error_low,error_high,sigma_limit)



# histogram
final_hist = final_pdf.createHistogram('final_hist',rec,RooFit.Binning(int((rec.getMax()-rec.getMin())*10)),RooFit.YVar(ion,RooFit.Binning(int((ion.getMax()-ion.getMin())*10))))


# RooFit frames
recbins = int((rec.getMax()-rec.getMin())*2)
ionbins = int((ion.getMax()-ion.getMin())*5)


ionframe = ion.frame()
ionframe.SetTitle('Projection in E_{ion}')
realdata.plotOn(ionframe, RooFit.Name('data'), RooFit.Binning(ionbins), RooFit.MarkerSize(1.0))
final_pdf.plotOn(ionframe, RooFit.Components("flat_gamma_bckgd_pdf"), RooFit.LineColor(kGreen), RooFit.LineWidth(2))
final_pdf.plotOn(ionframe, RooFit.Name('model'), RooFit.LineColor(kBlue), RooFit.LineWidth(2))
final_pdf.plotOn(ionframe, RooFit.Components("V49_ion_pdf"), RooFit.LineColor(kRed), RooFit.LineWidth(2), RooFit.LineStyle(kDashed))
final_pdf.plotOn(ionframe, RooFit.Components("Cr51_ion_pdf"), RooFit.LineColor(kRed), RooFit.LineWidth(2), RooFit.LineStyle(kDashed))
final_pdf.plotOn(ionframe, RooFit.Components("Mn54_ion_pdf"), RooFit.LineColor(kRed), RooFit.LineWidth(2), RooFit.LineStyle(kDashed))
final_pdf.plotOn(ionframe, RooFit.Components("Fe55_ion_pdf"), RooFit.LineColor(kRed), RooFit.LineWidth(2), RooFit.LineStyle(kDashed))
final_pdf.plotOn(ionframe, RooFit.Components("Co57_ion_pdf"), RooFit.LineColor(kRed), RooFit.LineWidth(2), RooFit.LineStyle(kDashed))
final_pdf.plotOn(ionframe, RooFit.Components("Zn65_ion_pdf"), RooFit.LineColor(kRed), RooFit.LineWidth(2), RooFit.LineStyle(kDashed))
final_pdf.plotOn(ionframe, RooFit.Components("Ge68_ion_pdf"), RooFit.LineColor(kRed), RooFit.LineWidth(2), RooFit.LineStyle(kDashed))
final_pdf.plotOn(ionframe, RooFit.Components("Ga68_ion_pdf"), RooFit.LineColor(kRed), RooFit.LineWidth(2), RooFit.LineStyle(kDashed))
final_pdf.plotOn(ionframe, RooFit.Components("signal_pdf"), RooFit.LineColor(kMagenta), RooFit.LineWidth(3))
final_pdf.paramOn(ionframe,RooFit.Format('NEU',RooFit.AutoPrecision(2)),RooFit.Layout(0.1,0.55,0.9),RooFit.ShowConstants(kFALSE))
red_chi2_ion = ionframe.chiSquare("model", "data", ndf)
#print "reduced chi2 for ion:",red_chi2_ion


recframe = rec.frame()
recframe.SetTitle('Projection in E_{rec}')
realdata.plotOn(recframe, RooFit.Name("data"), RooFit.Binning(recbins), RooFit.MarkerColor(kBlack), RooFit.MarkerSize(1.0))
final_pdf.plotOn(recframe, RooFit.Components("flat_gamma_bckgd_pdf"), RooFit.LineColor(kGreen), RooFit.LineWidth(2))
final_pdf.plotOn(recframe, RooFit.Name('model'), RooFit.LineColor(kBlue), RooFit.LineWidth(2))
final_pdf.plotOn(recframe, RooFit.Components("V49_rec_pdf"), RooFit.LineColor(kRed), RooFit.LineWidth(2), RooFit.LineStyle(kDashed))
final_pdf.plotOn(recframe, RooFit.Components("Cr51_rec_pdf"), RooFit.LineColor(kRed), RooFit.LineWidth(2), RooFit.LineStyle(kDashed))
final_pdf.plotOn(recframe, RooFit.Components("Mn54_rec_pdf"), RooFit.LineColor(kRed), RooFit.LineWidth(2), RooFit.LineStyle(kDashed))
final_pdf.plotOn(recframe, RooFit.Components("Fe55_rec_pdf"), RooFit.LineColor(kRed), RooFit.LineWidth(2), RooFit.LineStyle(kDashed))
final_pdf.plotOn(recframe, RooFit.Components("Co57_rec_pdf"), RooFit.LineColor(kRed), RooFit.LineWidth(2), RooFit.LineStyle(kDashed))
final_pdf.plotOn(recframe, RooFit.Components("Zn65_rec_pdf"), RooFit.LineColor(kRed), RooFit.LineWidth(2), RooFit.LineStyle(kDashed))
final_pdf.plotOn(recframe, RooFit.Components("Ge68_rec_pdf"), RooFit.LineColor(kRed), RooFit.LineWidth(2), RooFit.LineStyle(kDashed))
final_pdf.plotOn(recframe, RooFit.Components("Ga68_rec_pdf"), RooFit.LineColor(kRed), RooFit.LineWidth(2), RooFit.LineStyle(kDashed))
final_pdf.plotOn(recframe, RooFit.Components("signal_pdf"), RooFit.LineColor(kMagenta), RooFit.LineWidth(3))
red_chi2_rec = recframe.chiSquare("model", "data", ndf)
#print "reduced chi2 for rec:",red_chi2_rec


# Plotting of fit results (PDFs and projected spectra)
c1 = TCanvas('c1','Fit Results',1000,750)
c1.Divide(2,2)
# final pdf and data
c1.cd(1)
c1.cd(1).SetLogz()
final_hist.Draw('COL')
final_hist.SetStats(0)
final_hist.SetTitle('ID3 WIMP search data + best fit PDF')
realdata_graph.SetMarkerColor(kMagenta)
realdata_graph.SetMarkerStyle(kPlus)
realdata_graph.Draw('SAMESP')
ER_centroid.SetLineColor(kBlack)
ER_centroid.SetLineWidth(1)
ER_centroid.DrawCopy('SAME')
NR_centroid.SetLineColor(kBlack)
NR_centroid.SetLineWidth(1)
NR_centroid.DrawCopy('SAME')
# wimp signal and data only
if wimp_mass:
  c1.cd(3)
  signal_hist.Draw('CONT0')
  signal_hist.SetStats(0)
  signal_hist.SetContour(30)
  realdata_graph.Draw('SAMESP')
  ER_centroid.DrawCopy('SAME')
  NR_centroid.DrawCopy('SAME')
# data and fit in ionization energy
c1.cd(2)
ionframe.Draw()
# data and fit in recoil energy
c1.cd(4)
recframe.Draw()


# Monte Carlo statistics output
if MC_sets:
  MC_study = RooMCStudy(final_pdf,RooArgSet(rec,ion),RooFit.Verbose(kFALSE),RooFit.FitOptions(RooFit.PrintLevel(-1),RooFit.NumCPU(2)),RooFit.Extended(kTRUE))
  MC_study.generateAndFit(MC_sets,events,kTRUE)#keep generated data

  FitParams = FitResult.floatParsFinal()
  NumFitParams = FitParams.getSize()

  ParamNLLFrameList = []
  ParamDistriFrameList = []
  ParamPullFrameList = []
  ParamLineList = []

  c2 = TCanvas('c2','Statistical interpretation (MC study)')
  c2.Divide(1,NumFitParams,0.001,0.001)

  nll_line = TLine()
  nll_line.SetLineWidth(2)
  nll_line.SetLineStyle(7)
  nll_line.SetLineColor(kBlue)

  zeroline = TLine()
  zeroline.SetLineWidth(2)
  zeroline.SetLineStyle(7)
  zeroline.SetLineColor(kBlack)

  for i in range(NumFitParams):
    parameter = FitParams[i]
    paramname = parameter.GetName()

    paramline = TLine()
    paramline.SetLineWidth(2)
    paramline.SetLineStyle(1)
    paramline.SetLineColor(kBlue)
    ParamLineList.append(paramline)

    pad = c2.cd(i+1)
    pad.Divide(3,1,0.001,0.001)

    if paramname == 'N_signal':
      pad.SetFillColor(kYellow-9)
      pad.SetFillStyle(4100)

    paramnllframe = parameter.frame()
    paramnllframe.SetTitle('NLL fit '+paramname)
    nll.plotOn(paramnllframe, RooFit.Precision(1e-5),RooFit.ShiftToZero())
    paramnllframe.SetMaximum(20.)
    paramnllframe.SetMinimum(0.)
    ParamNLLFrameList.append(paramnllframe)
    pad.cd(1)
    paramnllframe.Draw()
    gPad.Update()
    paramline.DrawLine(parameter.getVal(),gPad.GetUymin(),parameter.getVal(),gPad.GetUymax())

    paramdistriframe = MC_study.plotParam(parameter)
    paramdistriframe.SetTitle('MC distri '+paramname)
    ParamDistriFrameList.append(paramdistriframe)
    pad.cd(2)
    paramdistriframe.Draw()
    paramdistriframe.getHist().Fit('gaus','QEM')
    gPad.Update()
    paramline.DrawLine(parameter.getVal(),gPad.GetUymin(),parameter.getVal(),gPad.GetUymax())

    parampullframe = MC_study.plotPull(parameter,RooFit.FitGauss(kTRUE))
    parampullframe.SetTitle('MC pull distri '+paramname)
    ParamPullFrameList.append(parampullframe)
    pad.cd(3)
    parampullframe.Draw()
    gPad.Update()
    zeroline.DrawLine(0,gPad.GetUymin(),0,gPad.GetUymax())

  c2.SetCanvasSize(1200,2400)

  # extra canvas for signal parameter only
  if wimp_mass:
    parameter = FitResult.floatParsFinal().find('N_signal')
    paramname = parameter.GetName()

    nllvalue = nll.getVal()

    c3 = TCanvas('c3','Fit Statistics WIMP Signal',1000,750)
    c3.Divide(2,2)
    c3.cd(1)
    paramnllframe = parameter.frame()
    paramnllframe.SetTitle('NLL fit '+paramname)
    nll.plotOn(paramnllframe, RooFit.Precision(1e-5), RooFit.ShiftToZero())
    paramnllframe.SetMaximum(20.)
    paramnllframe.SetMinimum(0.)
    paramnllframe.Draw()
    gPad.Update()
    paramline.DrawLine(parameter.getVal(),gPad.GetUymin(),parameter.getVal(),gPad.GetUymax())

    c3.cd(2)
    paramdistriframe = MC_study.plotParam(parameter)
    paramdistriframe.SetTitle('MC distri '+paramname)
    paramdistriframe.Draw()
    paramdistriframe.getHist().Fit('gaus','QEM')
    gPad.Update()
    paramline.DrawLine(parameter.getVal(),gPad.GetUymin(),parameter.getVal(),gPad.GetUymax())

    c3.cd(3)
    parampullframe = (MC_study.plotPull(parameter,RooFit.FitGauss(kTRUE)))
    parampullframe.SetTitle('MC pull distri '+paramname)
    parampullframe.Draw()
    gPad.Update()
    zeroline.DrawLine(0,gPad.GetUymin(),0,gPad.GetUymax())

    c3.cd(4)
    MCnllframe = MC_study.plotNLL()
    MCnllframe.Draw()
    MCnllframe.getHist().Fit('gaus','QEM')
    gPad.Update()
    nll_line.SetLineStyle(1)
    nll_line.DrawLine(nllvalue,gPad.GetUymin(),nllvalue,gPad.GetUymax())
    paramline.DrawLine(parameter.getVal(),gPad.GetUymin(),parameter.getVal(),gPad.GetUymax())


if SavePlots:
  c1.SaveAs('%iGeV_PDFs.png'%wimp_mass)
  if MC_sets:
    c2.SaveAs('%iGeV_param-stats.png'%wimp_mass)
    c3.SaveAs('%iGeV_signal-stats.png'%wimp_mass)


N_sig_list = []
for i in range(MC_sets):
  value = MC_study.fitParams(i).find('N_signal').getVal()
  print i, value 
  N_sig_list.append(value)
N_sig_array = np.array(N_sig_list, float)
N_sig_array.sort()
c = N_sig_array[0.9*MC_sets]
print N_sig_array[0.9*MC_sets]