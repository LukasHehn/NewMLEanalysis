#!/usr/bin/env python
from ROOT import *
from Functions import *
from DetectorClass import *





# definition of observables
ion = RooRealVar('ion','E_{ion}',0,15,'keV_{ee}')
rec = RooRealVar('rec','E_{rec}',0,30,'keV_{nr}')
time = RooRealVar('time','time',0.0,1.2,'years')


# detector efficiency
DetectorName = 'ID3'
ID = Detector(DetectorName)
total_efficiency = ID.GetEnergyEfficiency()
efficiency_datahist = RooDataHist('efficiency_datahist','efficiency_datahist',RooArgList(rec,ion),total_efficiency)
efficiency_pdf = RooHistPdf('efficiency_pdf','efficiency_pdf',RooArgSet(rec,ion),efficiency_datahist)


# detector specific parameters
voltage = RooRealVar('voltage','applied voltage',6.4)
voltage.setConstant()

FWHM_heat = 0.71
FWHM_rec = RecoilResolutionFromHeat(FWHM_heat,voltage.getVal(),10)
sigma_rec = RooRealVar('sigma_rec','recoil energy resolution',FWHM_rec/2.35)

FWHM_ion = 0.69
sigma_ion = RooRealVar('sigma_ion','ionization energy resolution',FWHM_ion/2.35)


# dataset
realdata = RooDataSet.read('ID3-eventlist_30keV.txt',RooArgList(time,rec,ion))
realdata_scatter = realdata.createHistogram(rec,ion)
events = int(realdata.numEntries())


# smearing in recoil
mean_rec = RooRealVar('mean_rec','mean_rec mean',0)
recoil_smearing = RooGaussian('recoil_smearing','recoil_smearing',rec,mean_rec,sigma_rec)
rec.setBins(1000,'fft')

# position of bands
centroid = RooFormulaVar('ER_centroid','@0*(1+(0.16*@0^0.18)*(@1/3))/(1+@1/3)',RooArgList(rec,voltage))

# gaussians in ionization energy with shifting mean in recoil energy
ion_gauss_ER = RooGaussian('ion_gauss_ER','gauss with shifted mean',ion,centroid,sigma_ion)

gamma_linear = RooPolynomial('gamma_linear','gamma_linear',ion)
gamma_spectrum = gamma_linear

# -----------------------------------------------------------------------------------------
# here may be dragons
energy_shift = RooRealVar('energy_shift','energy_shift',0.5,1.5)

V49_peak_ion_energy = RooRealVar('V49_peak_ion_energy','V49_peak_ion_energy',4.97)
V49_peak_ion_pos = RooFormulaVar('V49_peak_ion_pos','@0*@1',RooArgList(V49_peak_ion_energy,energy_shift))
V49_peak_ion_sigma = RooRealVar('V49_peak_ion_sigma','V49_peak_ion_sigma',sigma_ion.getVal())
V49_peak_ion = RooGaussian('V49_peak_ion_pdf','V49_peak_ion_pdf with shifted mean',ion,V49_peak_ion_pos,V49_peak_ion_sigma)
V49_peak_rec_pos = RooFormulaVar('V49_peak_rec_pos','@0*(1+@1/3)/(1+0.24*(@1/3))',RooArgList(V49_peak_ion_pos,voltage))
V49_peak_rec_sigma = RooRealVar('V49_peak_rec_sigma','V49_peak_rec_sigma',sigma_rec.getVal())
V49_peak_rec = RooGaussian('V49_peak_rec_pdf','V49_peak_rec_pdf with shifted mean',rec,V49_peak_rec_pos,V49_peak_rec_sigma)
V49_peak_pdf = RooProdPdf('V49_peak_pdf','5keV peak',V49_peak_ion,V49_peak_rec)
V49_peak_coeff = RooRealVar('V49_peak_coeff','scaling factor of 5keV peak',0.5,0.0,1.0)

Cr51_peak_ion_energy = RooRealVar('Cr51_peak_ion_energy','Cr51_peak_ion_energy',5.46)
Cr51_peak_ion_pos = RooFormulaVar('Cr51_peak_ion_pos','@0*@1',RooArgList(Cr51_peak_ion_energy,energy_shift))
Cr51_peak_ion_sigma = RooRealVar('Cr51_peak_ion_sigma','Cr51_peak_ion_sigma',sigma_ion.getVal())
Cr51_peak_ion = RooGaussian('Cr51_peak_ion_pdf','Cr51_peak_ion_pdf with shifted mean',ion,Cr51_peak_ion_pos,Cr51_peak_ion_sigma)
Cr51_peak_rec_pos = RooFormulaVar('Cr51_peak_rec_pos','@0*(1+@1/3)/(1+0.24*(@1/3))',RooArgList(Cr51_peak_ion_pos,voltage))
Cr51_peak_rec_sigma = RooRealVar('Cr51_peak_rec_sigma','Cr51_peak_rec_sigma',sigma_rec.getVal())
Cr51_peak_rec = RooGaussian('Cr51_peak_rec_pdf','Cr51_peak_rec_pdf with shifted mean',rec,Cr51_peak_rec_pos,Cr51_peak_rec_sigma)
Cr51_peak_pdf = RooProdPdf('Cr51_peak_pdf','5keV peak',Cr51_peak_ion,Cr51_peak_rec)
Cr51_peak_coeff = RooRealVar('Cr51_peak_coeff','scaling factor of 5keV peak',0.5,0.0,1.0)

Mn54_peak_ion_energy = RooRealVar('Mn54_peak_ion_energy','Mn54_peak_ion_energy',5.99)
Mn54_peak_ion_pos = RooFormulaVar('Mn54_peak_ion_pos','@0*@1',RooArgList(Mn54_peak_ion_energy,energy_shift))
Mn54_peak_ion_sigma = RooRealVar('Mn54_peak_ion_sigma','Mn54_peak_ion_sigma',sigma_ion.getVal())
Mn54_peak_ion = RooGaussian('Mn54_peak_ion_pdf','Mn54_peak_ion_pdf with shifted mean',ion,Mn54_peak_ion_pos,Mn54_peak_ion_sigma)
Mn54_peak_rec_pos = RooFormulaVar('Mn54_peak_rec_pos','@0*(1+@1/3)/(1+0.24*(@1/3))',RooArgList(Mn54_peak_ion_pos,voltage))
Mn54_peak_rec_sigma = RooRealVar('Mn54_peak_rec_sigma','Mn54_peak_rec_sigma',sigma_rec.getVal())
Mn54_peak_rec = RooGaussian('Mn54_peak_rec_pdf','Mn54_peak_rec_pdf with shifted mean',rec,Mn54_peak_rec_pos,Mn54_peak_rec_sigma)
Mn54_peak_pdf = RooProdPdf('Mn54_peak_pdf','5keV peak',Mn54_peak_ion,Mn54_peak_rec)
Mn54_peak_coeff = RooRealVar('Mn54_peak_coeff','scaling factor of 5keV peak',0.5,0.0,1.0)

Fe55_peak_ion_energy = RooRealVar('Fe55_peak_ion_energy','Fe55_peak_ion_energy',6.54)
Fe55_peak_ion_pos = RooFormulaVar('Fe55_peak_ion_pos','@0*@1',RooArgList(Fe55_peak_ion_energy,energy_shift))
Fe55_peak_ion_sigma = RooRealVar('Fe55_peak_ion_sigma','Fe55_peak_ion_sigma',sigma_ion.getVal())
Fe55_peak_ion = RooGaussian('Fe55_peak_ion_pdf','Fe55_peak_ion_pdf with shifted mean',ion,Fe55_peak_ion_pos,Fe55_peak_ion_sigma)
Fe55_peak_rec_pos = RooFormulaVar('Fe55_peak_rec_pos','@0*(1+@1/3)/(1+0.24*(@1/3))',RooArgList(Fe55_peak_ion_pos,voltage))
Fe55_peak_rec_sigma = RooRealVar('Fe55_peak_rec_sigma','Fe55_peak_rec_sigma',sigma_rec.getVal())
Fe55_peak_rec = RooGaussian('Fe55_peak_rec_pdf','Fe55_peak_rec_pdf with shifted mean',rec,Fe55_peak_rec_pos,Fe55_peak_rec_sigma)
Fe55_peak_pdf = RooProdPdf('Fe55_peak_pdf','6keV peak',Fe55_peak_ion,Fe55_peak_rec)
Fe55_peak_coeff = RooRealVar('Fe55_peak_coeff','scaling factor of 6keV peak',0.5,0.0,1.0)


Co56_peak_ion_energy = RooRealVar('Co56_peak_ion_energy','Co56_peak_ion_energy',7.11)
Co56_peak_ion_pos = RooFormulaVar('Co56_peak_ion_pos','@0*@1',RooArgList(Co56_peak_ion_energy,energy_shift))
Co56_peak_ion_sigma = RooRealVar('Co56_peak_ion_sigma','Co56_peak_ion_sigma',sigma_ion.getVal())
Co56_peak_ion = RooGaussian('Co56_peak_ion_pdf','Co56_peak_ion_pdf with shifted mean',ion,Co56_peak_ion_pos,Co56_peak_ion_sigma)
Co56_peak_rec_pos = RooFormulaVar('Co56_peak_rec_pos','@0*(1+@1/3)/(1+0.24*(@1/3))',RooArgList(Co56_peak_ion_pos,voltage))
Co56_peak_rec_sigma = RooRealVar('Co56_peak_rec_sigma','Co56_peak_rec_sigma',sigma_rec.getVal())
Co56_peak_rec = RooGaussian('Co56_peak_rec_pdf','Co56_peak_rec_pdf with shifted mean',rec,Co56_peak_rec_pos,Co56_peak_rec_sigma)
Co56_peak_pdf = RooProdPdf('Co56_peak_pdf','7keV peak',Co56_peak_ion,Co56_peak_rec)
Co56_peak_coeff = RooRealVar('Co56_peak_coeff','scaling factor of 7keV peak',0.5,0.0,1.0)

Zn65_peak_ion_energy = RooRealVar('Zn65_peak_ion_energy','Zn65_peak_ion_energy',8.98)
Zn65_peak_ion_pos = RooFormulaVar('Zn65_peak_ion_pos','@0*@1',RooArgList(Zn65_peak_ion_energy,energy_shift))
Zn65_peak_ion_sigma = RooRealVar('Zn65_peak_ion_sigma','Zn65_peak_ion_sigma',sigma_ion.getVal())
Zn65_peak_ion = RooGaussian('Zn65_peak_ion_pdf','Zn65_peak_ion_pdf with shifted mean',ion,Zn65_peak_ion_pos,Zn65_peak_ion_sigma)
Zn65_peak_rec_pos = RooFormulaVar('Zn65_peak_rec_pos','@0*(1+@1/3)/(1+0.24*(@1/3))',RooArgList(Zn65_peak_ion_pos,voltage))
Zn65_peak_rec_sigma = RooRealVar('Zn65_peak_rec_sigma','Zn65_peak_rec_sigma',sigma_rec.getVal())
Zn65_peak_rec = RooGaussian('Zn65_peak_rec_pdf','Zn65_peak_rec_pdf with shifted mean',rec,Zn65_peak_rec_pos,Zn65_peak_rec_sigma)
Zn65_peak_pdf = RooProdPdf('Zn65_peak_pdf','9keV peak',Zn65_peak_ion,Zn65_peak_rec)
Zn65_peak_coeff = RooRealVar('Zn65_peak_coeff','scaling factor of 9keV peak',0.5,0.0,1.0)

Ge68_peak_ion_energy = RooRealVar('Ge68_peak_ion_energy','Ge68_peak_ion_energy',10.37)
Ge68_peak_ion_pos = RooFormulaVar('Ge68_peak_ion_pos','@0*@1',RooArgList(Ge68_peak_ion_energy,energy_shift))
Ge68_peak_ion_sigma = RooRealVar('Ge68_peak_ion_sigma','Ge68_peak_ion_sigma',sigma_ion.getVal())
Ge68_peak_ion = RooGaussian('Ge68_peak_ion_pdf','Ge68_peak_ion_pdf with shifted mean',ion,Ge68_peak_ion_pos,Ge68_peak_ion_sigma)
Ge68_peak_rec_pos = RooFormulaVar('Ge68_peak_rec_pos','@0*(1+@1/3)/(1+0.24*(@1/3))',RooArgList(Ge68_peak_ion_pos,voltage))
#Ge68_peak_rec_pos = RooRealVar('Ge68_peak_rec_pos', 'Ge68_peak_rec_pos', GetEnergyRecoilFromEstimator(Ge68_peak_ion_pos.getVal(),6.4))
Ge68_peak_rec_sigma = RooRealVar('Ge68_peak_rec_sigma','Ge68_peak_rec_sigma',sigma_rec.getVal())
Ge68_peak_rec = RooGaussian('Ge68_peak_rec_pdf','Ge68_peak_rec_pdf with shifted mean',rec,Ge68_peak_rec_pos,Ge68_peak_rec_sigma)
Ge68_peak_pdf = RooProdPdf('Ge68_peak_pdf','10keV peak',Ge68_peak_ion,Ge68_peak_rec)
Ge68_peak_coeff = RooRealVar('Ge68_peak_coeff','scaling factor of 10keV peak',0.5,0.0,1.0)

Ge68_2_peak_ion_energy = RooRealVar('Ge68_2_peak_ion_energy','Ge68_2_peak_ion_energy',9.66)
Ge68_2_peak_ion_pos = RooFormulaVar('Ge68_2_peak_ion_pos','@0*@1',RooArgList(Ge68_2_peak_ion_energy,energy_shift))
Ge68_2_peak_ion_sigma = RooRealVar('Ge68_2_peak_ion_sigma','Ge68_2_peak_ion_sigma',sigma_ion.getVal())
Ge68_2_peak_ion = RooGaussian('Ge68_2_peak_ion_pdf','Ge68_2_peak_ion_pdf with shifted mean',ion,Ge68_2_peak_ion_pos,Ge68_2_peak_ion_sigma)
Ge68_2_peak_rec_pos = RooFormulaVar('Ge68_2_peak_rec_pos','@0*(1+@1/3)/(1+0.24*(@1/3))',RooArgList(Ge68_2_peak_ion_pos,voltage))
Ge68_2_peak_rec_sigma = RooRealVar('Ge68_2_peak_rec_sigma','Ge68_2_peak_rec_sigma',sigma_rec.getVal())
Ge68_2_peak_rec = RooGaussian('Ge68_2_peak_rec_pdf','Ge68_2_peak_rec_pdf with shifted mean',rec,Ge68_2_peak_rec_pos,Ge68_2_peak_rec_sigma)
Ge68_2_peak_pdf = RooProdPdf('Ge68_2_peak_pdf','96keV peak',Ge68_2_peak_ion,Ge68_2_peak_rec)
Ge68_2_peak_coeff = RooRealVar('Ge68_2_peak_coeff','scaling factor of 96keV peak',0.5,0.0,1.0)
## -----------------------------------------------------------------------------------------

# 2D model
#gamma_pdf = RooProdPdf('gamma_pdf','gamma_pdf',gamma_spectrum,ion_gauss_ER)
#gamma_pdf_smeared = RooFFTConvPdf('gamma_pdf_smeared','gamma pdf * heat gaussian',rec,gamma_pdf,recoil_smearing)
#gamma_pdf_smeared.setBufferFraction(0.9)
#final_gamma_pdf = RooAddPdf('final_gamma_pdf','gamma band description',RooArgList(peak5_pdf,peak6_pdf,peak7_pdf,gamma_pdf_smeared),RooArgList(peak5_coeff,peak6_coeff,peak7_coeff),kTRUE)
#final_gamma_pdf.fitTo(realdata,RooFit.NumCPU(2))


# 1D model ion
constant_ion = RooPolynomial('constant bckgd ion', 'constant bckgd', ion)
gamma_bckgd_ion = RooAddPdf('gamma_bckgd_ion','gamma_bckgd_ion',RooArgList(Cr51_peak_ion, Fe55_peak_ion, Co56_peak_ion, Zn65_peak_ion, Ge68_peak_ion, Ge68_2_peak_ion,constant_ion),RooArgList(Cr51_peak_coeff, Fe55_peak_coeff, Co56_peak_coeff, Zn65_peak_coeff, Ge68_peak_coeff,Ge68_2_peak_coeff),kTRUE)
#gamma_bckgd_ion = RooAddPdf('gamma_bckgd_ion','gamma_bckgd_ion',RooArgList(V49_peak_ion, Cr51_peak_ion, Mn54_peak_ion, Fe55_peak_ion, Co56_peak_ion, Zn65_peak_ion, Ge68_peak_ion, Ge68_2_peak_ion,constant_ion),RooArgList(V49_peak_coeff, Cr51_peak_coeff, Mn54_peak_coeff, Fe55_peak_coeff, Co56_peak_coeff, Zn65_peak_coeff, Ge68_peak_coeff,Ge68_2_peak_coeff),kTRUE)
gamma_bckgd_ion.fitTo(realdata,RooFit.NumCPU(2))
final_gamma_pdf = RooProdPdf('final_gamma_pdf','final_gamma_pdf',RooArgList(gamma_bckgd_ion, efficiency_pdf))


## 1D model rec
#constant_rec = RooPolynomial('constant bckgd rec', 'constant bckgd', rec)
#gamma_bckgd_rec = RooAddPdf('gamma_bckgd_rec','gamma_bckgd_rec',RooArgList(V49_peak_rec, Cr51_peak_rec, Mn54_peak_rec, Fe55_peak_rec, Co56_peak_rec, Zn65_peak_rec, Ge68_peak_rec, Ge68_2_peak_rec,constant_rec),RooArgList(V49_peak_coeff, Cr51_peak_coeff, Mn54_peak_coeff, Fe55_peak_coeff, Co56_peak_coeff, Zn65_peak_coeff, Ge68_peak_coeff,Ge68_2_peak_coeff),kTRUE)
#gamma_bckgd_rec.fitTo(realdata,RooFit.NumCPU(2))
#final_gamma_pdf = gamma_bckgd_rec


# histogram
gamma_hist = final_gamma_pdf.createHistogram('gamma_hist',rec,RooFit.Binning(int((rec.getMax()-rec.getMin())*20)),RooFit.YVar(ion,RooFit.Binning(int((ion.getMax()-ion.getMin())*20))))


recbins = int((rec.getMax()-rec.getMin())*10)
ionbins = int((ion.getMax()-ion.getMin())*10)



# RooFit frames
#recframe = rec.frame()
#realdata.plotOn(recframe, RooFit.Name("data"), RooFit.Binning(recbins))
#final_gamma_pdf.plotOn(recframe, RooFit.Name('model'), RooFit.LineColor(kBlue))
#final_gamma_pdf.plotOn(recframe, RooFit.Components("V49_peak_rec_pdf"), RooFit.LineColor(kRed), RooFit.LineStyle(kDashed))
#final_gamma_pdf.plotOn(recframe, RooFit.Components("Cr51_peak_rec_pdf"), RooFit.LineColor(kRed), RooFit.LineStyle(kDashed))
#final_gamma_pdf.plotOn(recframe, RooFit.Components("Mn54_peak_rec_pdf"), RooFit.LineColor(kRed), RooFit.LineStyle(kDashed))
#final_gamma_pdf.plotOn(recframe, RooFit.Components("Fe55_peak_rec_pdf"), RooFit.LineColor(kRed), RooFit.LineStyle(kDashed))
#final_gamma_pdf.plotOn(recframe, RooFit.Components("Co56_peak_rec_pdf"), RooFit.LineColor(kRed), RooFit.LineStyle(kDashed))
#final_gamma_pdf.plotOn(recframe, RooFit.Components("Zn65_peak_rec_pdf"), RooFit.LineColor(kRed), RooFit.LineStyle(kDashed))
#final_gamma_pdf.plotOn(recframe, RooFit.Components("Ge68_peak_rec_pdf"), RooFit.LineColor(kRed), RooFit.LineStyle(kDashed))
#final_gamma_pdf.plotOn(recframe, RooFit.Components("Ge68_2_peak_rec_pdf"), RooFit.LineColor(kRed), RooFit.LineStyle(kDashed))
#print recframe.chiSquare("model", "data", 9)

ionframe = ion.frame()
realdata.plotOn(ionframe, RooFit.Name('data'), RooFit.Binning(ionbins))
final_gamma_pdf.plotOn(ionframe, RooFit.Name('model'), RooFit.LineColor(kBlue))
final_gamma_pdf.plotOn(ionframe, RooFit.Components("V49_peak_ion_pdf"), RooFit.LineColor(kRed), RooFit.LineStyle(kDashed))
final_gamma_pdf.plotOn(ionframe, RooFit.Components("Cr51_peak_ion_pdf"), RooFit.LineColor(kRed), RooFit.LineStyle(kDashed))
final_gamma_pdf.plotOn(ionframe, RooFit.Components("Mn54_peak_ion_pdf"), RooFit.LineColor(kRed), RooFit.LineStyle(kDashed))
final_gamma_pdf.plotOn(ionframe, RooFit.Components("Fe55_peak_ion_pdf"), RooFit.LineColor(kRed), RooFit.LineStyle(kDashed))
final_gamma_pdf.plotOn(ionframe, RooFit.Components("Co56_peak_ion_pdf"), RooFit.LineColor(kRed), RooFit.LineStyle(kDashed))
final_gamma_pdf.plotOn(ionframe, RooFit.Components("Zn65_peak_ion_pdf"), RooFit.LineColor(kRed), RooFit.LineStyle(kDashed))
final_gamma_pdf.plotOn(ionframe, RooFit.Components("Ge68_peak_ion_pdf"), RooFit.LineColor(kRed), RooFit.LineStyle(kDashed))
final_gamma_pdf.plotOn(ionframe, RooFit.Components("Ge68_2_peak_ion_pdf"), RooFit.LineColor(kRed), RooFit.LineStyle(kDashed))
print ionframe.chiSquare("model", "data", 9)

c0 = TCanvas('c0','gamma band',800,600)
c0.Divide(1,2)
c0.cd(1)
gamma_hist.Draw('COLZ')
realdata_scatter.Draw('SAMES')
realdata_scatter.SetMarkerColor(kBlack)
ER_centroid.SetParameter(0,voltage.getVal())
ER_centroid.Draw('SAMES')
ER_centroid.SetLineColor(kBlack)
ER_centroid.SetLineWidth(2)
NR_centroid.Draw('SAMES')
NR_centroid.SetLineColor(kBlack)
NR_centroid.SetLineWidth(2)
c0.cd(2)
ionframe.Draw()
#recframe.Draw()
