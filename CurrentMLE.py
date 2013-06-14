#!/usr/bin/env python
from ROOT import *
from Functions import *
from DetectorClass import *


# load Eric's WIMP signal file
gROOT.LoadMacro("/kalinka/home/hehn/PhD/LowMassEric/WimpDistri.C")


#switches and input parameters to control script
wimp_mass = False #set wimp mass or switch signal of entirely
MC_sims = False #set number of MC simulations: 0 means none at all
cutset = False #use event set with 3 outlying events cut


# definition of observables
ion = RooRealVar('ion','E_{ion}',Energy['ion']['min'],Energy['ion']['max'],'keV_{ee}')
rec = RooRealVar('rec','E_{rec}',Energy['rec']['min'],Energy['rec']['max'],'keV_{nr}')
time = RooRealVar('time','time',0.0,1.2,'years')


# detector efficiency
total_efficiency = Simple2DEfficiencyID3()

efficiency_datahist = RooDataHist('efficiency_datahist','efficiency_datahist',RooArgList(rec,ion),total_efficiency)
efficiency_pdf = RooHistPdf('efficiency_pdf','efficiency_pdf',RooArgSet(rec,ion),efficiency_datahist)# detector specific parameters
voltage = RooRealVar('voltage','applied voltage',6.4)

#FWHM_heat = 0.71 #me@baseline
FWHM_heat = 0.82 #Eric@10keV
FWHM_rec = RecoilResolutionFromHeat(FWHM_heat,voltage.getVal(),10)
sigma_rec = RooRealVar('sigma_rec','recoil energy resolution',FWHM_rec/2.35)

#FWHM_ion = 0.69 #me
FWHM_ion = 0.72 #Eric
sigma_ion = RooRealVar('sigma_ion','ionization energy resolution',FWHM_ion/2.35)


# dataset
if cutset: realdata = RooDataSet.read('ID3-eventlist_30keV_cut3.txt',RooArgList(time,rec,ion)) #without 2 events near NR band
else: realdata = RooDataSet.read('ID3-eventlist_30keV.txt',RooArgList(time,rec,ion))
realdata_scatter = realdata.createHistogram(rec,ion,Energy['rec']['bins'],Energy['ion']['bins'])
realdata_graph = TGraphFromDataSet(realdata)
events = int(realdata.numEntries())


# -----------------------------------------------------------------------------------------
# definition of gamma peaks
energy_correction_ion = RooRealVar('energy_correction_ion','scaling factor ionization energy',0.5,1.5)
energy_correction_rec = RooRealVar('energy_correction_rec','scaling factor recoil energy',0.5,1.5)
ER_centroid.SetParameter(0,6.4) #ER_centroid used for calculation of peak position in Erec

V49_ion_energy = RooRealVar('V49_ion_energy','V49 peak ion energy',4.97)
V49_ion_pos = RooFormulaVar('V49_ion_pos','@0*@1',RooArgList(V49_ion_energy,energy_correction_ion))
V49_ion_pdf = RooGaussian('V49_ion_pdf','V49 peak gauss pdf in ion',ion,V49_ion_pos,sigma_ion)
V49_rec_energy = RooRealVar('V49_rec_energy','recoil energy V49 peak',ER_centroid.GetX(V49_ion_energy.getVal()))
V49_rec_pos = RooFormulaVar('V49_rec_pos','@0*@1',RooArgList(V49_rec_energy,energy_correction_rec))
V49_rec_pdf = RooGaussian('V49_rec_pdf','V49 peak pdf in recoil energy',rec,V49_rec_pos,sigma_rec)
V49_pdf = RooProdPdf('V49_pdf','V49 peak pdf',V49_ion_pdf,V49_rec_pdf)
V49_coeff = RooRealVar('V49_coeff','fraction of V49 peak (4.97keV)',0.5,0.0,1.0)


Cr51_ion_energy = RooRealVar('Cr51_ion_energy','Cr51 peak ion energy',5.46)
Cr51_ion_pos = RooFormulaVar('Cr51_ion_pos','@0*@1',RooArgList(Cr51_ion_energy,energy_correction_ion))
Cr51_ion = RooGaussian('Cr51_ion_pdf','Cr51 peak gauss pdf with shifted mean',ion,Cr51_ion_pos,sigma_ion)
Cr51_rec_energy = RooRealVar('Cr51_rec_energy','v_rec_energy',ER_centroid.GetX(Cr51_ion_energy.getVal()))
Cr51_rec_pos = RooFormulaVar('Cr51_rec_pos','@0*@1',RooArgList(Cr51_rec_energy,energy_correction_rec))
Cr51_rec = RooGaussian('Cr51_rec_pdf','Cr51_rec_pdf with shifted mean',rec,Cr51_rec_pos,sigma_rec)
Cr51_pdf = RooProdPdf('Cr51_pdf','Cr51 peak pdf',Cr51_ion,Cr51_rec)
Cr51_coeff = RooRealVar('Cr51_coeff','fraction of 51Cr peak (5.46keV)',0.5,0.0,1.0)


Mn54_ion_energy = RooRealVar('Mn54_ion_energy','Mn54_ion_energy',5.99)
Mn54_ion_pos = RooFormulaVar('Mn54_ion_pos','@0*@1',RooArgList(Mn54_ion_energy,energy_correction_ion))
Mn54_ion = RooGaussian('Mn54_ion_pdf','Mn54_ion_pdf with shifted mean',ion,Mn54_ion_pos,sigma_ion)
Mn54_rec_energy = RooRealVar('Mn54_rec_energy','Mn54_rec_energy',ER_centroid.GetX(Mn54_ion_energy.getVal()))
Mn54_rec_pos = RooFormulaVar('Mn54_rec_pos','@0*@1',RooArgList(Mn54_rec_energy,energy_correction_rec))
Mn54_rec = RooGaussian('Mn54_rec_pdf','Mn54_rec_pdf with shifted mean',rec,Mn54_rec_pos,sigma_rec)
Mn54_pdf = RooProdPdf('Mn54_pdf','Mn54 peak pdf',Mn54_ion,Mn54_rec)
Mn54_coeff = RooRealVar('Mn54_coeff','fraction of 54Mn peak (5.99keV)',0.5,0.0,1.0)


Fe55_ion_energy = RooRealVar('Fe55_ion_energy','Fe55_ion_energy',6.54)
Fe55_ion_pos = RooFormulaVar('Fe55_ion_pos','@0*@1',RooArgList(Fe55_ion_energy,energy_correction_ion))
Fe55_ion = RooGaussian('Fe55_ion_pdf','Fe55_ion_pdf with shifted mean',ion,Fe55_ion_pos,sigma_ion)
Fe55_rec_energy = RooRealVar('Fe55_rec_energy','Fe55_rec_energy',ER_centroid.GetX(Fe55_ion_energy.getVal()))
Fe55_rec_pos = RooFormulaVar('Fe55_rec_pos','@0*@1',RooArgList(Fe55_rec_energy,energy_correction_rec))
Fe55_rec = RooGaussian('Fe55_rec_pdf','Fe55_rec_pdf with shifted mean',rec,Fe55_rec_pos,sigma_rec)
Fe55_pdf = RooProdPdf('Fe55_pdf','Fe55 peak pdf',Fe55_ion,Fe55_rec)
Fe55_coeff = RooRealVar('Fe55_coeff','fraction of 55Fe peak (6.54keV)',0.5,0.0,1.0)


Co57_ion_energy = RooRealVar('Co57_ion_energy','Co57_ion_energy',7.11)
Co57_ion_pos = RooFormulaVar('Co57_ion_pos','@0*@1',RooArgList(Co57_ion_energy,energy_correction_ion))
Co57_ion = RooGaussian('Co57_ion_pdf','Co57_ion_pdf with shifted mean',ion,Co57_ion_pos,sigma_ion)
Co57_rec_energy = RooRealVar('Co57_rec_energy','Co57_rec_energy',ER_centroid.GetX(Co57_ion_energy.getVal()))
Co57_rec_pos = RooFormulaVar('Co57_rec_pos','@0*@1',RooArgList(Co57_rec_energy,energy_correction_rec))
Co57_rec = RooGaussian('Co57_rec_pdf','Co57_rec_pdf with shifted mean',rec,Co57_rec_pos,sigma_rec)
Co57_pdf = RooProdPdf('Co57_pdf','Co57 peak pdf',Co57_ion,Co57_rec)
Co57_coeff = RooRealVar('Co57_coeff','fraction of 57Co peak (7.11keV)',0.5,0.0,1.0)


Zn65_ion_energy = RooRealVar('Zn65_ion_energy','Zn65_ion_energy',8.98)
Zn65_ion_pos = RooFormulaVar('Zn65_ion_pos','@0*@1',RooArgList(Zn65_ion_energy,energy_correction_ion))
Zn65_ion = RooGaussian('Zn65_ion_pdf','Zn65_ion_pdf with shifted mean',ion,Zn65_ion_pos,sigma_ion)
Zn65_rec_energy = RooRealVar('Zn65_rec_energy','Zn65_rec_energy',ER_centroid.GetX(Zn65_ion_energy.getVal()))
Zn65_rec_pos = RooFormulaVar('Zn65_rec_pos','@0*@1',RooArgList(Zn65_rec_energy,energy_correction_rec))
Zn65_rec = RooGaussian('Zn65_rec_pdf','Zn65_rec_pdf with shifted mean',rec,Zn65_rec_pos,sigma_rec)
Zn65_pdf = RooProdPdf('Zn65_pdf','Zn65 peak pdf',Zn65_ion,Zn65_rec)
Zn65_coeff = RooRealVar('Zn65_coeff','fraction of 65Zn peak (8.98keV)',0.5,0.0,1.0)


Ga68_ion_energy = RooRealVar('Ga68_ion_energy','Ga68_ion_energy',9.66)
Ga68_ion_pos = RooFormulaVar('Ga68_ion_pos','@0*@1',RooArgList(Ga68_ion_energy,energy_correction_ion))
Ga68_ion = RooGaussian('Ga68_ion_pdf','Ga68_ion_pdf with shifted mean',ion,Ga68_ion_pos,sigma_ion)
Ga68_rec_energy = RooRealVar('Ga68_rec_energy','Ga68_rec_energy',ER_centroid.GetX(Ga68_ion_energy.getVal()))
Ga68_rec_pos = RooFormulaVar('Ga68_rec_pos','@0*@1',RooArgList(Ga68_rec_energy,energy_correction_rec))
Ga68_rec = RooGaussian('Ga68_rec_pdf','Ga68_rec_pdf with shifted mean',rec,Ga68_rec_pos,sigma_rec)
Ga68_pdf = RooProdPdf('Ga68_pdf','Ga68 peak pdf',Ga68_ion,Ga68_rec)
Ga68_coeff = RooRealVar('Ga68_coeff','fraction of 68Ga peak (9.66keV)',0.5,0.0,1.0)


Ge68_ion_energy = RooRealVar('Ge68_ion_energy','Ge68_ion_energy',10.37)
Ge68_ion_pos = RooFormulaVar('Ge68_ion_pos','@0*@1',RooArgList(Ge68_ion_energy,energy_correction_ion))
Ge68_ion = RooGaussian('Ge68_ion_pdf','Ge68_ion_pdf with shifted mean',ion,Ge68_ion_pos,sigma_ion)
Ge68_rec_energy = RooRealVar('Ge68_rec_energy','Ge68_rec_energy',ER_centroid.GetX(Ge68_ion_energy.getVal()))
Ge68_rec_pos = RooFormulaVar('Ge68_rec_pos','@0*@1',RooArgList(Ge68_rec_energy,energy_correction_rec))
Ge68_rec = RooGaussian('Ge68_rec_pdf','Ge68 peak pdf in rec',rec,Ge68_rec_pos,sigma_rec)
Ge68_pdf = RooProdPdf('Ge68_pdf','Ge68 peak pdf',Ge68_ion,Ge68_rec)
Ge68_coeff = RooRealVar('Ge68_coeff','fraction of 68Ge peak (10.37keV)',0.5,0.0,1.0)
# -----------------------------------------------------------------------------------------


# wimp signal
if wimp_mass:
  TriggerEfficiency.SetParameter(0, 3.874)
  TriggerEfficiency.SetParameter(1, FWHM_rec)

  TriggerEfficiency.SetNpx(1000)
  efficiency = TriggerEfficiency.GetHistogram()

  signal_hist = WimpDistri(str(wimp_mass), 'ID3', FWHM_rec, FWHM_ion, efficiency, 0, 0, 0, 6.4, 1)
  signal_hist.SetTitle('WIMP signal %sGeV'%wimp_mass)
  signal_datahist = RooDataHist('signal_datahist','signal_datahist',RooArgList(rec,ion),signal_hist)
  signal_pdf = RooHistPdf('signal_pdf','signal_pdf',RooArgSet(rec,ion),signal_datahist)


# flat gamma background component
flat_gamma_bckgd_hist = FlatGammaBckgd2DEric(sigma_ion.getVal(),sigma_rec.getVal())
flat_gamma_bckgd_hist.Multiply(total_efficiency)
flat_gamma_bckgd_datahist = RooDataHist('flat_gamma_bckgd_datahist','flat_gamma_bckgd_datahist',RooArgList(rec,ion),flat_gamma_bckgd_hist)
flat_gamma_bckgd_pdf = RooHistPdf('flat_gamma_bckgd_pdf','flat_gamma_bckgd_pdf',RooArgSet(rec,ion),flat_gamma_bckgd_datahist)


# normal gamma bckgd model with peaks
gamma_bckgd_pdf = RooAddPdf('combined_bckgd_pdf','combined_bckgd_pdf',RooArgList(V49_pdf, Cr51_pdf, Mn54_pdf, Fe55_pdf, Co57_pdf, Zn65_pdf, Ga68_pdf, Ge68_pdf, flat_gamma_bckgd_pdf),RooArgList(V49_coeff, Cr51_coeff, Mn54_coeff, Fe55_coeff, Co57_coeff, Zn65_coeff, Ga68_coeff, Ge68_coeff),kTRUE)


## initial background fit only
#gamma_bckgd_pdf.fitTo(realdata)
#for param in [V49_coeff, Cr51_coeff, Mn54_coeff, Fe55_coeff, Co57_coeff, Zn65_coeff, Ga68_coeff, Ge68_coeff, energy_correction_ion,energy_correction_rec]:
  #param.setConstant(kTRUE)


if wimp_mass:
  # combine signal and background
  signal_ratio = RooRealVar('signal_ratio','signal_ratio',0.5,-0.1,1.0)
  final_pdf = RooAddPdf('final_pdf','final_pdf',signal_pdf,gamma_bckgd_pdf,signal_ratio)
  params = RooArgSet(signal_ratio,energy_correction_rec,energy_correction_ion)
else:
  final_pdf = gamma_bckgd_pdf
  params = RooArgSet(energy_correction_rec,energy_correction_ion)

# manual mode
nll = RooNLLVar('nll','nll',final_pdf,realdata,RooFit.PrintEvalErrors(2))
minuit = RooMinuit(nll)
FitResults = minuit.fit('hvr')


ndf = FitResults.floatParsFinal().getSize()


# histogram
final_hist = final_pdf.createHistogram('final_hist',rec,RooFit.Binning(int((rec.getMax()-rec.getMin())*10)),RooFit.YVar(ion,RooFit.Binning(int((ion.getMax()-ion.getMin())*10))))


recbins = int((rec.getMax()-rec.getMin())*5)
ionbins = int((ion.getMax()-ion.getMin())*10)



# RooFit frames
recframe = rec.frame()
realdata.plotOn(recframe, RooFit.Name("data"), RooFit.Binning(recbins), RooFit.MarkerColor(kBlack), RooFit.MarkerSize(1.0))
final_pdf.plotOn(recframe, RooFit.Name('model'), RooFit.LineColor(kBlue), RooFit.LineWidth(2))
final_pdf.plotOn(recframe, RooFit.Components("V49_rec_pdf"), RooFit.LineColor(kRed), RooFit.LineWidth(2), RooFit.LineStyle(kDashed))
final_pdf.plotOn(recframe, RooFit.Components("Cr51_rec_pdf"), RooFit.LineColor(kRed), RooFit.LineWidth(2), RooFit.LineStyle(kDashed))
final_pdf.plotOn(recframe, RooFit.Components("Mn54_rec_pdf"), RooFit.LineColor(kRed), RooFit.LineWidth(2), RooFit.LineStyle(kDashed))
final_pdf.plotOn(recframe, RooFit.Components("Fe55_rec_pdf"), RooFit.LineColor(kRed), RooFit.LineWidth(2), RooFit.LineStyle(kDashed))
final_pdf.plotOn(recframe, RooFit.Components("Co57_rec_pdf"), RooFit.LineColor(kRed), RooFit.LineWidth(2), RooFit.LineStyle(kDashed))
final_pdf.plotOn(recframe, RooFit.Components("Zn65_rec_pdf"), RooFit.LineColor(kRed), RooFit.LineWidth(2), RooFit.LineStyle(kDashed))
final_pdf.plotOn(recframe, RooFit.Components("Ge68_rec_pdf"), RooFit.LineColor(kRed), RooFit.LineWidth(2), RooFit.LineStyle(kDashed))
final_pdf.plotOn(recframe, RooFit.Components("Ga68_rec_pdf"), RooFit.LineColor(kRed), RooFit.LineWidth(2), RooFit.LineStyle(kDashed))
final_pdf.paramOn(recframe,RooFit.Parameters(params),RooFit.Format('NEU',RooFit.AutoPrecision(1)),RooFit.Layout(0.1,0.5,0.9),RooFit.ShowConstants(kTRUE))
red_chi2_rec = recframe.chiSquare("model", "data", ndf)
print "reduced chi2 for rec:",red_chi2_rec


ionframe = ion.frame()
realdata.plotOn(ionframe, RooFit.Name('data'), RooFit.Binning(ionbins), RooFit.MarkerSize(1.0))
final_pdf.plotOn(ionframe, RooFit.Name('model'), RooFit.LineColor(kBlue), RooFit.LineWidth(2))
final_pdf.plotOn(ionframe, RooFit.Components("V49_ion_pdf"), RooFit.LineColor(kRed), RooFit.LineWidth(2), RooFit.LineStyle(kDashed))
final_pdf.plotOn(ionframe, RooFit.Components("Cr51_ion_pdf"), RooFit.LineColor(kRed), RooFit.LineWidth(2), RooFit.LineStyle(kDashed))
final_pdf.plotOn(ionframe, RooFit.Components("Mn54_ion_pdf"), RooFit.LineColor(kRed), RooFit.LineWidth(2), RooFit.LineStyle(kDashed))
final_pdf.plotOn(ionframe, RooFit.Components("Fe55_ion_pdf"), RooFit.LineColor(kRed), RooFit.LineWidth(2), RooFit.LineStyle(kDashed))
final_pdf.plotOn(ionframe, RooFit.Components("Co57_ion_pdf"), RooFit.LineColor(kRed), RooFit.LineWidth(2), RooFit.LineStyle(kDashed))
final_pdf.plotOn(ionframe, RooFit.Components("Zn65_ion_pdf"), RooFit.LineColor(kRed), RooFit.LineWidth(2), RooFit.LineStyle(kDashed))
final_pdf.plotOn(ionframe, RooFit.Components("Ge68_ion_pdf"), RooFit.LineColor(kRed), RooFit.LineWidth(2), RooFit.LineStyle(kDashed))
final_pdf.plotOn(ionframe, RooFit.Components("Ga68_ion_pdf"), RooFit.LineColor(kRed), RooFit.LineWidth(2), RooFit.LineStyle(kDashed))
final_pdf.paramOn(ionframe,RooFit.Parameters(params),RooFit.Format('NEU',RooFit.AutoPrecision(1)),RooFit.Layout(0.1,0.5,0.9),RooFit.ShowConstants(kTRUE))
red_chi2_ion = ionframe.chiSquare("model", "data", ndf)
red_chi2_ion_label = TPaveLabel()
red_chi2_ion_label.SetLabel('test')
ionframe.addObject(red_chi2_ion_label)
print "reduced chi2 for ion:",red_chi2_ion


FitResults.Print('v')


#print "wimp mass:",wimp_mass
#print "signal ratio:",signal_ratio.getVal()
#print "upper error:",signal_ratio.getErrorHi()
#print "90% C.L. upper events:",(signal_ratio.getVal()+1.64*signal_ratio.getErrorHi())*events


# plotting
c1 = TCanvas('c1','fit results',1000,750)
c1.Divide(2,2)
# final pdf and data
c1.cd(1)
c1.cd(1).SetLogz()
final_hist.Draw('COLZ')
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
  ratioframe = signal_ratio.frame()
  nll.plotOn(ratioframe)
  c1.cd(2)
  signal_hist.Draw('COLZ')
  signal_hist.SetStats(0)
  final_hist.SetTitle('ID3 WIMP search data + signal PDF only')
  realdata_graph.Draw('SAMESP')
  ER_centroid.DrawCopy('SAME')
  NR_centroid.DrawCopy('SAME')
# data and fit in ionization energy
c1.cd(3)
ionframe.Draw()
ionframe.SetTitle('Projection in E_{ion}')
ionframe.GetXaxis().SetRangeUser(1,13.5)
# data and fit in recoil energy
c1.cd(4)
recframe.Draw()
recframe.SetTitle('Projection in E_{rec}')
recframe.GetXaxis().SetRangeUser(3,25)


# Monte Carlo statistics output
# Monte Carlo
if MC_sims:
  MC_study = RooMCStudy(final_pdf,RooArgSet(rec,ion))
  MC_study.generateAndFit(MC_sims,events,kTRUE)
  nllframe = MC_study.plotNLL()
  paramframe = MC_study.plotParam(signal_ratio)#,RooFit.FrameRange(ratio_range['min'],ratio_range['max']),RooFit.Binning(150))

  c2 = TCanvas('c2','Monte Carlo statitics',1000,750)
  c2.Divide(2,2)
  # MC NLL value distribution
  c2.cd(1)
  nllframe.Draw()
  nllframe.SetTitle('MC distribution of NLL-values')
  nll_line = TLine()
  nll_line.SetLineWidth(2)
  nll_line.SetLineColor(kRed)
  nll_line.DrawLine(nll.getVal(),gPad.GetUymin(),nll.getVal(),gPad.GetUymax())
  gPad.Update()
  # NLL function of fit to real data set
  c2.cd(2)
  ratioframe.Draw()
  ratioframe.SetTitle('NLL')
  gPad.Update()
  nll_line.DrawLine(gPad.GetUxmin(),nll.getVal(),gPad.GetUxmax(),nll.getVal())
  ratio_line = TLine()
  ratio_line.SetLineWidth(2)
  ratio_line.SetLineColor(kRed)
  ratio_line.DrawLine(signal_ratio.getVal(),gPad.GetUymin(),signal_ratio.getVal(),gPad.GetUymax())
  ratio_line.SetLineStyle(7)
  ratio_line.DrawLine(signal_ratio.getVal()+signal_ratio.getErrorHi(),gPad.GetUymin(),signal_ratio.getVal()+signal_ratio.getErrorHi(),gPad.GetUymax())
  # MC ratio distribution
  c2.cd(3)
  paramframe.Draw()
  ratio_line.SetLineWidth(2)
  ratio_line.SetLineColor(kRed)
  ratio_line.DrawLine(signal_ratio.getVal(),gPad.GetUymin(),signal_ratio.getVal(),gPad.GetUymax())
  ratio_line.SetLineStyle(7)
  ratio_line.DrawLine(signal_ratio.getVal()+signal_ratio.getErrorHi(),gPad.GetUymin(),signal_ratio.getVal()+signal_ratio.getErrorHi(),gPad.GetUymax())
