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

ID = Detector('ID3')


FWHM_ion = ID.GetAverageValue('fiducialmean') #0.72
FWHM_heat = ID.GetAverageValue('heat') #0.82
voltage = ID.GetAverageValue('voltage') #6.4
threshold_ee = ID.GetAverageValue('threshold')

FWHM_rec = RecoilResolutionFromHeatBaseline(FWHM_heat,voltage,10)
sigma_ion = FWHM_ion/2.35
sigma_rec = FWHM_rec/2.35
threshold_nr = GetEnergyRecoilFromEstimator(threshold_ee,voltage)


List, Spectrum, Integral = ReadInWimpSpectrumEric(mass_of_wimp)
Wimp_hist = WimpSignal2DEric(mass_of_wimp, sigma_ion,sigma_rec,Spectrum)
Wimp_hist.Scale(Integral/Wimp_hist.Integral('WIDTH'))
print "wimp hist integral",Wimp_hist.Integral('WIDTH')
Signal = Wimp_hist.Clone('wimp_signal')
efficiency = ID.GetProjectedEnergyEfficiency()
gammacut_efficiency = ID.GetGammaCutEfficiency()
Signal.Multiply(efficiency)
Signal.Multiply(gammacut_efficiency)
print "signal hist integral",Signal.Integral('WIDTH')
