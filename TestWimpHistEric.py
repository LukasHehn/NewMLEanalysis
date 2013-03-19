#!/usr/bin/env python
from ROOT import *
from Parameters import *
from Functions import *
from DetectorClass import *
import numpy as np
import pyWIMP.DMModels.wimp_model as wimp_model
from pyWIMP.DMModels.flat_model import FlatModel
from pyWIMP.DMModels.base_model import BaseVariables

mass_of_wimp = 10

#ID = Detector('ID3')


def WimpSignal2DEric(sigma_ion,sigma_rec,spectrum):
  denom_i=1./(2*pow(sigma_ion,2))
  denom_r=1./(2*pow(sigma_rec,2))
  #hist = TH2F('densityhist','densityhist',Energy['rec']['bins'],Energy['rec']['min'],Energy['rec']['max'],Energy['ion']['bins'],Energy['ion']['min'],Energy['ion']['max'])
  hist = TH2F('densityhist','densityhist',100,-20,20,100,-10,10)
  for recbin in range(1,hist.GetNbinsX()+1):
    Erec = hist.GetXaxis().GetBinCenter(recbin)
    for ionbin in range(1,hist.GetNbinsY()+1):
      Eion = hist.GetYaxis().GetBinCenter(ionbin)
      summe=0
      for specbin in range(1,spectrum.GetNbinsX()+1):
	ErecSpec = spectrum.GetXaxis().GetBinCenter(specbin)
	Q = LindhardQuenching.Eval(ErecSpec)
	wimprate = spectrum.GetBinContent(specbin)
	kernel=TMath.Exp(-denom_r*pow((Erec-ErecSpec),2)-denom_i*pow((Eion-Q*ErecSpec),2))
	summe+=(kernel*wimprate)
      hist.SetBinContent(recbin,ionbin,0.02002*summe/(2*3.141592*sigma_rec*sigma_ion))
  return hist

FWHM_ion = 0.72
FWHM_heat = 0.82
voltage = 6.4
FWHM_rec = RecoilResolutionFromHeatBaseline(FWHM_heat,voltage,10)
sigma_ion = FWHM_ion/2.35
sigma_rec = FWHM_rec/2.35
threshold_nr = GetEnergyRecoilFromEstimator(1.99,6.4)


def ReadInWimpSpectrum(mass_of_wimp):
  basedir="/kalinka/home/hehn/PhD/LowMassEric/"
  filename=basedir+"wimpsignal_M"+str(mass_of_wimp)+".txt"

  Energies,Rates = [], []

  Integral = 0

  print "reading WIMP spectrum from:",filename
  print "energy  |  rate"

  infile = open(filename,'r')
  for line in infile:
    energy, rate = line.split()
    Integral += float(rate)
    Energies.append(float(energy))
    Rates.append(float(rate))
    print energy, rate
  infile.close()

  DeltaE = Energies[1]-Energies[0]
  Integral *= DeltaE

  Energybins = np.array(Energies, dtype=np.float)

  Hist = TH1F('wimp_spectrum_%sGeV'%mass_of_wimp,'Spectrum for 10^{-6}pb and %sGeV;E_{recoil} [keV];Rate [evts/kg day 0.02keV]'%mass_of_wimp,Energybins.size-1, Energybins.flatten('C'))

  for i in range(len(Energies)):
    Hist.Fill(Energies[i], Rates[i])

  print "integrated rate:",Integral

  return [Energies, Rates], Hist, Integral