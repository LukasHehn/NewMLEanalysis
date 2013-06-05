from ROOT import *
from Functions import *

gROOT.LoadMacro("/kalinka/home/hehn/PhD/LowMassEric/WimpDistri.C")

FWHM_heat = 0.82
FWHM_rec = RecoilResolutionFromHeat(FWHM_heat,6.4,10)
FWHM_ion = 0.72

TriggerEfficiency.SetParameter(0, 3.874)
TriggerEfficiency.SetParameter(1, FWHM_rec)

TriggerEfficiency.SetNpx(1000)
efficiency = TriggerEfficiency.GetHistogram()

hist = WimpDistri('10', 'ID3', FWHM_rec, FWHM_ion, efficiency, 0, 0, 0, 6.4, 1)

hist.Draw('CONT0')

ER_centroid.SetLineColor(kMagenta)
ER_centroid.SetLineWidth(2)
ER_centroid.SetParameter(0,6.4)
ER_centroid.DrawCopy('SAME')

NR_centroid.SetLineColor(kMagenta)
NR_centroid.SetLineWidth(2)
NR_centroid.DrawCopy('SAME')
NR_centroid.Draw('SAME')

hist.GetXaxis().SetRangeUser(3.2,12.8)
hist.GetXaxis().SetTitle('E_{rec} (keV_{nr})')
hist.GetYaxis().SetRangeUser(0.8,5.8)
hist.GetYaxis().SetTitle('E_{ion} (keV_{ee})')
