#!/usr/bin/env python
from ROOT import *
from Functions import *
from DetectorClass import *


# load Eric's WIMP signal file
gROOT.LoadMacro("/kalinka/home/hehn/PhD/LowMassEric/WimpDistri.C")


# definition of observables
ion = RooRealVar('ion','E_{ion}',Energy['ion']['min'],Energy['ion']['max'],'keV_{ee}')
rec = RooRealVar('rec','E_{rec}',Energy['rec']['min'],Energy['rec']['max'],'keV_{nr}')
time = RooRealVar('time','time',0.0,1.2,'years')


# detector efficiency
total_efficiency = Simple2DEfficiencyID3()

efficiency_datahist = RooDataHist('efficiency_datahist','efficiency_datahist',RooArgList(rec,ion),total_efficiency)
efficiency_pdf = RooHistPdf('efficiency_pdf','efficiency_pdf',RooArgSet(rec,ion),efficiency_datahist)


# detector specific parameters
voltage = RooRealVar('voltage','applied voltage',6.4)

#FWHM_heat = 0.71 #me@baseline
FWHM_heat = 0.82 #Eric@10keV
FWHM_rec = RecoilResolutionFromHeat(FWHM_heat,voltage.getVal(),10)
sigma_rec = RooRealVar('sigma_rec','recoil energy resolution',FWHM_rec/2.35)

#FWHM_ion = 0.69 #me
FWHM_ion = 0.72 #Eric
sigma_ion = RooRealVar('sigma_ion','ionization energy resolution',FWHM_ion/2.35)


# dataset
#realdata = RooDataSet.read('ID3-eventlist_30keV_cut3.txt',RooArgList(time,rec,ion)) #without 2 events near NR band
realdata = RooDataSet.read('ID3-eventlist_30keV.txt',RooArgList(time,rec,ion))
realdata_scatter = realdata.createHistogram(rec,ion,Energy['rec']['bins'],Energy['ion']['bins'])
realdata_graph = TGraphFromDataSet(realdata)
events = int(realdata.numEntries())


## smearing in recoil
#mean_rec = RooRealVar('mean_rec','mean_rec mean',0)
#recoil_smearing = RooGaussian('recoil_smearing','recoil_smearing',rec,mean_rec,sigma_rec)
#rec.setBins(10000,'fft')

# position of bands
#ER_centroid = RooFormulaVar('ER_centroid','@0*(1+(0.16*@0^0.18)*(@1/3))/(1+@1/3)',RooArgList(rec,voltage))
#ER_centroid = RooFormulaVar('ER_centroid','0.5*@0',RooArgList(rec))
#ER_centroid = RooFit.bindFunction(ER_centroid_real,rec,RooArgList(voltage)) #use only real defined function for Electron Recoil centroid
#NR_centroid = RooFormulaVar('NR_centroid','0.16*@0^1.18',RooArgList(rec))

# gaussians in ionization energy with shifting mean in recoil energy
#ion_gauss_ER = RooGaussian('ion_gauss_ER','gauss with shifted mean',ion,ER_centroid,sigma_ion)


# -----------------------------------------------------------------------------------------
# definition of gamma peaks
energy_correction_ion = RooRealVar('energy_correction_ion','energy_correction_ion',0.5,1.5)
energy_correction_rec = RooRealVar('energy_correction_rec','energy_correction_rec',0.5,1.5)
ER_centroid.SetParameter(0,6.4) #ER_centroid used for calculation of peak position in Erec

V49_ion_energy = RooRealVar('V49_ion_energy','V49_ion_energy',4.97)
V49_ion_pos = RooFormulaVar('V49_ion_pos','@0*@1',RooArgList(V49_ion_energy,energy_correction_ion))
V49_ion = RooGaussian('V49_ion_pdf','V49_ion_pdf with shifted mean',ion,V49_ion_pos,sigma_ion)
V49_rec_energy = RooRealVar('V49_rec_energy','V49_rec_energy',ER_centroid.GetX(V49_ion_energy.getVal()))
V49_rec_pos = RooFormulaVar('V49_rec_pos','@0*@1',RooArgList(V49_rec_energy,energy_correction_rec))
V49_rec = RooGaussian('V49_rec_pdf','V49_rec_pdf with shifted mean',rec,V49_rec_pos,sigma_rec)
V49_pdf = RooProdPdf('V49_pdf','5keV peak',V49_ion,V49_rec)
V49_coeff = RooRealVar('V49_coeff','scaling factor of 5keV peak',0.5,0.0,1.0)
#n_V49 = RooRealVar('n_V49','events for 49V peak',0.,0.,1000.0)

Cr51_ion_energy = RooRealVar('Cr51_ion_energy','Cr51_ion_energy',5.46)
Cr51_ion_pos = RooFormulaVar('Cr51_ion_pos','@0*@1',RooArgList(Cr51_ion_energy,energy_correction_ion))
Cr51_ion = RooGaussian('Cr51_ion_pdf','Cr51_ion_pdf with shifted mean',ion,Cr51_ion_pos,sigma_ion)
Cr51_rec_energy = RooRealVar('Cr51_rec_energy','v_rec_energy',ER_centroid.GetX(Cr51_ion_energy.getVal()))
Cr51_rec_pos = RooFormulaVar('Cr51_rec_pos','@0*@1',RooArgList(Cr51_rec_energy,energy_correction_rec))
Cr51_rec = RooGaussian('Cr51_rec_pdf','Cr51_rec_pdf with shifted mean',rec,Cr51_rec_pos,sigma_rec)
Cr51_pdf = RooProdPdf('Cr51_pdf','5keV peak',Cr51_ion,Cr51_rec)
Cr51_coeff = RooRealVar('Cr51_coeff','scaling factor of 5keV peak',0.5,0.0,1.0)
#n_Cr51 = RooRealVar('n_Cr51','events for 51Cr peak',0.,0.,1000.0)

Mn54_ion_energy = RooRealVar('Mn54_ion_energy','Mn54_ion_energy',5.99)
Mn54_ion_pos = RooFormulaVar('Mn54_ion_pos','@0*@1',RooArgList(Mn54_ion_energy,energy_correction_ion))
Mn54_ion = RooGaussian('Mn54_ion_pdf','Mn54_ion_pdf with shifted mean',ion,Mn54_ion_pos,sigma_ion)
Mn54_rec_energy = RooRealVar('Mn54_rec_energy','Mn54_rec_energy',ER_centroid.GetX(Mn54_ion_energy.getVal()))
Mn54_rec_pos = RooFormulaVar('Mn54_rec_pos','@0*@1',RooArgList(Mn54_rec_energy,energy_correction_rec))
Mn54_rec = RooGaussian('Mn54_rec_pdf','Mn54_rec_pdf with shifted mean',rec,Mn54_rec_pos,sigma_rec)
Mn54_pdf = RooProdPdf('Mn54_pdf','5keV peak',Mn54_ion,Mn54_rec)
Mn54_coeff = RooRealVar('Mn54_coeff','scaling factor of 5keV peak',0.5,0.0,1.0)
#n_Mn54 = RooRealVar('n_Mn54','events for 54Mn peak',0.,0.,1000.0)

Fe55_ion_energy = RooRealVar('Fe55_ion_energy','Fe55_ion_energy',6.54)
Fe55_ion_pos = RooFormulaVar('Fe55_ion_pos','@0*@1',RooArgList(Fe55_ion_energy,energy_correction_ion))
Fe55_ion = RooGaussian('Fe55_ion_pdf','Fe55_ion_pdf with shifted mean',ion,Fe55_ion_pos,sigma_ion)
Fe55_rec_energy = RooRealVar('Fe55_rec_energy','Fe55_rec_energy',ER_centroid.GetX(Fe55_ion_energy.getVal()))
Fe55_rec_pos = RooFormulaVar('Fe55_rec_pos','@0*@1',RooArgList(Fe55_rec_energy,energy_correction_rec))
Fe55_rec = RooGaussian('Fe55_rec_pdf','Fe55_rec_pdf with shifted mean',rec,Fe55_rec_pos,sigma_rec)
Fe55_pdf = RooProdPdf('Fe55_pdf','6keV peak',Fe55_ion,Fe55_rec)
Fe55_coeff = RooRealVar('Fe55_coeff','scaling factor of 6keV peak',0.5,0.0,1.0)
#n_Fe55 = RooRealVar('n_Fe55','events for 55Fe peak',0.,0.,1000.0)

Co57_ion_energy = RooRealVar('Co57_ion_energy','Co57_ion_energy',7.11)
Co57_ion_pos = RooFormulaVar('Co57_ion_pos','@0*@1',RooArgList(Co57_ion_energy,energy_correction_ion))
Co57_ion = RooGaussian('Co57_ion_pdf','Co57_ion_pdf with shifted mean',ion,Co57_ion_pos,sigma_ion)
Co57_rec_energy = RooRealVar('Co57_rec_energy','Co57_rec_energy',ER_centroid.GetX(Co57_ion_energy.getVal()))
Co57_rec_pos = RooFormulaVar('Co57_rec_pos','@0*@1',RooArgList(Co57_rec_energy,energy_correction_rec))
Co57_rec = RooGaussian('Co57_rec_pdf','Co57_rec_pdf with shifted mean',rec,Co57_rec_pos,sigma_rec)
Co57_pdf = RooProdPdf('Co57_pdf','7keV peak',Co57_ion,Co57_rec)
Co57_coeff = RooRealVar('Co57_coeff','scaling factor of 7keV peak',0.5,0.0,1.0)
#n_Co57 = RooRealVar('n_Co57','events for 56Co peak',0.,0.,1000.0)

Zn65_ion_energy = RooRealVar('Zn65_ion_energy','Zn65_ion_energy',8.98)
Zn65_ion_pos = RooFormulaVar('Zn65_ion_pos','@0*@1',RooArgList(Zn65_ion_energy,energy_correction_ion))
Zn65_ion = RooGaussian('Zn65_ion_pdf','Zn65_ion_pdf with shifted mean',ion,Zn65_ion_pos,sigma_ion)
Zn65_rec_energy = RooRealVar('Zn65_rec_energy','Zn65_rec_energy',ER_centroid.GetX(Zn65_ion_energy.getVal()))
Zn65_rec_pos = RooFormulaVar('Zn65_rec_pos','@0*@1',RooArgList(Zn65_rec_energy,energy_correction_rec))
Zn65_rec = RooGaussian('Zn65_rec_pdf','Zn65_rec_pdf with shifted mean',rec,Zn65_rec_pos,sigma_rec)
Zn65_pdf = RooProdPdf('Zn65_pdf','9keV peak',Zn65_ion,Zn65_rec)
Zn65_coeff = RooRealVar('Zn65_coeff','scaling factor of 9keV peak',0.5,0.0,1.0)
#n_Zn65 = RooRealVar('n_Zn65','events for 65Zn peak',0.,0.,1000.0)

Ge68_ion_energy = RooRealVar('Ge68_ion_energy','Ge68_ion_energy',10.37)
Ge68_ion_pos = RooFormulaVar('Ge68_ion_pos','@0*@1',RooArgList(Ge68_ion_energy,energy_correction_ion))
Ge68_ion = RooGaussian('Ge68_ion_pdf','Ge68_ion_pdf with shifted mean',ion,Ge68_ion_pos,sigma_ion)
Ge68_rec_energy = RooRealVar('Ge68_rec_energy','Ge68_rec_energy',ER_centroid.GetX(Ge68_ion_energy.getVal()))
Ge68_rec_pos = RooFormulaVar('Ge68_rec_pos','@0*@1',RooArgList(Ge68_rec_energy,energy_correction_rec))
Ge68_rec = RooGaussian('Ge68_rec_pdf','Ge68_rec_pdf with shifted mean',rec,Ge68_rec_pos,sigma_rec)
Ge68_pdf = RooProdPdf('Ge68_pdf','10keV peak',Ge68_ion,Ge68_rec)
#Ge68_coeff = RooRealVar('Ge68_coeff','scaling factor of 10keV peak',0.5,0.0,1.0)
#n_Ge68 = RooRealVar('n_Ge68','events for 68Ge peak',0.,0.,1000.0)

Ga68_ion_energy = RooRealVar('Ga68_ion_energy','Ga68_ion_energy',9.66)
Ga68_ion_pos = RooFormulaVar('Ga68_ion_pos','@0*@1',RooArgList(Ga68_ion_energy,energy_correction_ion))
Ga68_ion = RooGaussian('Ga68_ion_pdf','Ga68_ion_pdf with shifted mean',ion,Ga68_ion_pos,sigma_ion)
Ga68_rec_energy = RooRealVar('Ga68_rec_energy','Ga68_rec_energy',ER_centroid.GetX(Ga68_ion_energy.getVal()))
Ga68_rec_pos = RooFormulaVar('Ga68_rec_pos','@0*@1',RooArgList(Ga68_rec_energy,energy_correction_rec))
Ga68_rec = RooGaussian('Ga68_rec_pdf','Ga68_rec_pdf with shifted mean',rec,Ga68_rec_pos,sigma_rec)
Ga68_pdf = RooProdPdf('Ga68_pdf','96keV peak',Ga68_ion,Ga68_rec)
Ga68_coeff = RooRealVar('Ga68_coeff','scaling factor of 96keV peak',0.1,0.0,1.0)#0.1,0.0,1.0)
#Ga68_coeff.setConstant(kTRUE)
#n_Ga68_small = RooFormulaVar('n_Ga68_small','0.1*@0',RooArgList(n_Ga68))


GeGa68_pdf = RooAddPdf('GeGa68_pdf', '68Ge and fixed 68Ga peak',Ga68_pdf,Ge68_pdf,Ga68_coeff)
GeGa68_coeff = RooRealVar('GeGa68_coeff','GeGa68_coeff',0.0,1.0)

#n_flat = RooRealVar('n_flat','events flat spectrum',0.,0.,1000.0)
# -----------------------------------------------------------------------------------------


# wimp signal
TriggerEfficiency.SetParameter(0, 3.874)
TriggerEfficiency.SetParameter(1, FWHM_rec)

TriggerEfficiency.SetNpx(1000)
efficiency = TriggerEfficiency.GetHistogram()

signal_hist = WimpDistri('10', 'ID3', FWHM_rec, FWHM_ion, efficiency, 0, 0, 0, 6.4, 1)
signal_datahist = RooDataHist('signal_datahist','signal_datahist',RooArgList(rec,ion),signal_hist)
signal_pdf = RooHistPdf('signal_pdf','signal_pdf',RooArgSet(rec,ion),signal_datahist)


# gamma background
flat_gamma_bckgd_hist = FlatGammaBckgd2DEric(sigma_ion.getVal(),sigma_rec.getVal())
flat_gamma_bckgd_hist.Multiply(total_efficiency)
flat_gamma_bckgd_datahist = RooDataHist('flat_gamma_bckgd_datahist','flat_gamma_bckgd_datahist',RooArgList(rec,ion),flat_gamma_bckgd_hist)
flat_gamma_bckgd_pdf = RooHistPdf('flat_gamma_bckgd_pdf','flat_gamma_bckgd_pdf',RooArgSet(rec,ion),flat_gamma_bckgd_datahist)

# normal gamma bckgd model
gamma_bckgd_pdf = RooAddPdf('combined_bckgd_pdf','combined_bckgd_pdf',RooArgList(V49_pdf, Cr51_pdf, Mn54_pdf, Fe55_pdf, Co57_pdf, Zn65_pdf, GeGa68_pdf, flat_gamma_bckgd_pdf),RooArgList(V49_coeff, Cr51_coeff, Mn54_coeff, Fe55_coeff, Co57_coeff, Zn65_coeff, GeGa68_coeff),kTRUE)

# only selection of peaks
#selected_peaks_pdf = RooAddPdf('combined_bckgd_pdf','combined_bckgd_pdf',RooArgList(V49_pdf, Fe55_pdf, Co57_pdf, Zn65_pdf, Ge68_pdf, Ga68_pdf, flat_gamma_bckgd_pdf),RooArgList(V49_coeff, Fe55_coeff, Co57_coeff, Zn65_coeff, Ge68_coeff, Ga68_coeff),kTRUE)

#extended model
#combined_bckgd_pdf = RooAddPdf('combined_bckgd_pdf','combined_bckgd_pdf',RooArgList(V49_pdf, Cr51_pdf, Mn54_pdf, Fe55_pdf, Co57_pdf, Zn65_pdf, Ge68_pdf, Ga68_pdf, flat_gamma_bckgd_pdf),RooArgList(n_V49, n_Cr51, n_Mn54, n_Fe55, n_Co57, n_Zn65, n_Ge68, n_Ge68_small))


# maximum likelihood fit to data
# automatic mode
#FitResults = gamma_bckgd_pdf.fitTo(realdata)


# combine signal and background
signal_ratio = RooRealVar('signal_ratio','signal_ratio',0.0,0.0,1.0)
final_pdf = RooAddPdf('final_pdf','final_pdf',gamma_bckgd_pdf,signal_pdf,signal_ratio)

# manual mode
nll = RooNLLVar('nll_background_only','nll background only',final_pdf,realdata,RooFit.PrintEvalErrors(2))
minuit = RooMinuit(nll)
FitResults = minuit.fit('hvr')


ndf = FitResults.floatParsFinal().getSize()


# Monte Carlo
MC_study = RooMCStudy(gamma_bckgd_pdf,RooArgSet(rec,ion))
MC_study.generateAndFit(10,events,kTRUE)
nllframe = MC_study.plotNLL()


# histogram
gamma_bckgd_hist = gamma_bckgd_pdf.createHistogram('final_gamma_bckgd_hist',rec,RooFit.Binning(int((rec.getMax()-rec.getMin())*10)),RooFit.YVar(ion,RooFit.Binning(int((ion.getMax()-ion.getMin())*10))))


recbins = int((rec.getMax()-rec.getMin())*5)
ionbins = int((ion.getMax()-ion.getMin())*10)

parameters = RooArgSet(energy_correction_rec,energy_correction_ion,Ga68_coeff)

# RooFit frames
recframe = rec.frame()
realdata.plotOn(recframe, RooFit.Name("data"), RooFit.Binning(recbins), RooFit.MarkerColor(kBlack), RooFit.MarkerSize(1.0))
gamma_bckgd_pdf.plotOn(recframe, RooFit.Name('model'), RooFit.LineColor(kBlue), RooFit.LineWidth(2))
gamma_bckgd_pdf.plotOn(recframe, RooFit.Components("V49_rec_pdf"), RooFit.LineColor(kRed), RooFit.LineWidth(2), RooFit.LineStyle(kDashed))
gamma_bckgd_pdf.plotOn(recframe, RooFit.Components("Cr51_rec_pdf"), RooFit.LineColor(kRed), RooFit.LineWidth(2), RooFit.LineStyle(kDashed))
gamma_bckgd_pdf.plotOn(recframe, RooFit.Components("Mn54_rec_pdf"), RooFit.LineColor(kRed), RooFit.LineWidth(2), RooFit.LineStyle(kDashed))
gamma_bckgd_pdf.plotOn(recframe, RooFit.Components("Fe55_rec_pdf"), RooFit.LineColor(kRed), RooFit.LineWidth(2), RooFit.LineStyle(kDashed))
gamma_bckgd_pdf.plotOn(recframe, RooFit.Components("Co57_rec_pdf"), RooFit.LineColor(kRed), RooFit.LineWidth(2), RooFit.LineStyle(kDashed))
gamma_bckgd_pdf.plotOn(recframe, RooFit.Components("Zn65_rec_pdf"), RooFit.LineColor(kRed), RooFit.LineWidth(2), RooFit.LineStyle(kDashed))
gamma_bckgd_pdf.plotOn(recframe, RooFit.Components("Ge68_rec_pdf"), RooFit.LineColor(kRed), RooFit.LineWidth(2), RooFit.LineStyle(kDashed))
gamma_bckgd_pdf.plotOn(recframe, RooFit.Components("Ga68_rec_pdf"), RooFit.LineColor(kRed), RooFit.LineWidth(2), RooFit.LineStyle(kDashed))
#realdata.plotOn(recframe, RooFit.Name("data"), RooFit.Binning(recbins), RooFit.MarkerColor(kBlack), RooFit.MarkerSize(1.0))
gamma_bckgd_pdf.paramOn(recframe,RooFit.Parameters(parameters),RooFit.Format('NEU',RooFit.AutoPrecision(1)),RooFit.Layout(0.12,0.7,0.9),RooFit.ShowConstants(kTRUE))
#realdata.statOn(recframe)
red_chi2_rec = recframe.chiSquare("model", "data", ndf)
print "reduced chi2 for rec:",red_chi2_rec

ionframe = ion.frame()
realdata.plotOn(ionframe, RooFit.Name('data'), RooFit.Binning(ionbins), RooFit.MarkerSize(1.0))
gamma_bckgd_pdf.plotOn(ionframe, RooFit.Name('model'), RooFit.LineColor(kBlue), RooFit.LineWidth(2))
gamma_bckgd_pdf.plotOn(ionframe, RooFit.Components("V49_ion_pdf"), RooFit.LineColor(kRed), RooFit.LineWidth(2), RooFit.LineStyle(kDashed))
gamma_bckgd_pdf.plotOn(ionframe, RooFit.Components("Cr51_ion_pdf"), RooFit.LineColor(kRed), RooFit.LineWidth(2), RooFit.LineStyle(kDashed))
gamma_bckgd_pdf.plotOn(ionframe, RooFit.Components("Mn54_ion_pdf"), RooFit.LineColor(kRed), RooFit.LineWidth(2), RooFit.LineStyle(kDashed))
gamma_bckgd_pdf.plotOn(ionframe, RooFit.Components("Fe55_ion_pdf"), RooFit.LineColor(kRed), RooFit.LineWidth(2), RooFit.LineStyle(kDashed))
gamma_bckgd_pdf.plotOn(ionframe, RooFit.Components("Co57_ion_pdf"), RooFit.LineColor(kRed), RooFit.LineWidth(2), RooFit.LineStyle(kDashed))
gamma_bckgd_pdf.plotOn(ionframe, RooFit.Components("Zn65_ion_pdf"), RooFit.LineColor(kRed), RooFit.LineWidth(2), RooFit.LineStyle(kDashed))
gamma_bckgd_pdf.plotOn(ionframe, RooFit.Components("Ge68_ion_pdf"), RooFit.LineColor(kRed), RooFit.LineWidth(2), RooFit.LineStyle(kDashed))
gamma_bckgd_pdf.plotOn(ionframe, RooFit.Components("Ga68_ion_pdf"), RooFit.LineColor(kRed), RooFit.LineWidth(2), RooFit.LineStyle(kDashed))
#realdata.plotOn(ionframe, RooFit.Name('data'), RooFit.Binning(ionbins), RooFit.MarkerSize(1.0))
gamma_bckgd_pdf.paramOn(ionframe,RooFit.Parameters(parameters),RooFit.Format('NEU',RooFit.AutoPrecision(1)),RooFit.Layout(0.12,0.7,0.9),RooFit.ShowConstants(kTRUE))
red_chi2_ion = ionframe.chiSquare("model", "data", ndf)
red_chi2_ion_label = TPaveLabel()
red_chi2_ion_label.SetLabel('test')
ionframe.addObject(red_chi2_ion_label)
print "reduced chi2 for ion:",red_chi2_ion


# plotting
c0 = TCanvas('c0','gamma band',1000,750)
c0.Divide(2,2)
c0.cd(1)
gamma_bckgd_hist.Draw('COLZ')
gamma_bckgd_hist.SetStats(0)
gamma_bckgd_hist.SetTitle('ID3 WIMP search data + gamma bckgd PDF')
#realdata_scatter.SetMarkerColor(kBlack)
#realdata_scatter.SetMarkerStyle(kPlus)
#realdata_scatter.SetMarkerSize(0.8)
#realdata_scatter.Draw('SAMES')
realdata_graph.SetMarkerColor(kMagenta)
realdata_graph.SetMarkerStyle(kPlus)
realdata_graph.Draw('SAMESP')
c0.cd(1).SetLogz()
ER_centroid.SetLineColor(kBlack)
ER_centroid.SetLineWidth(1)
ER_centroid.DrawCopy('SAME')
NR_centroid.SetLineColor(kBlack)
NR_centroid.SetLineWidth(1)
NR_centroid.Draw('SAME')
c0.cd(2)
ionframe.Draw()
ionframe.SetTitle('Projection in E_{ion}')
ionframe.GetXaxis().SetRangeUser(1,13.5)
c0.cd(3)
recframe.Draw()
recframe.SetTitle('Projection in E_{rec}')
recframe.GetXaxis().SetRangeUser(3,25)
c0.cd(4)
line = TLine()
line.SetLineWidth(2)
line.SetLineColor(kRed)
nllframe.Draw()
nllframe.SetTitle('MC distribution of NLL-values')
gPad.Update()
line.DrawLine(nll.getVal(),gPad.GetUymin(),nll.getVal(),gPad.GetUymax())
