from ROOT import *
from Functions import *

gROOT.LoadMacro("/kalinka/home/hehn/PhD/LowMassEric/WimpDistri.C")

FWHM_heat = 0.82
FWHM_rec = RecoilResolutionFromHeat(FWHM_heat,6.4,10)
sigma_rec = FWHM_rec/2.35
FWHM_ion = 0.72

wimp_mass = '15'
detector = 'ID3'

TriggerEfficiency.SetParameter(0, 3.874)
TriggerEfficiency.SetParameter(1, sigma_rec)

TriggerEfficiency.SetNpx(1000)
efficiency = TriggerEfficiency.GetHistogram()

signal = WimpDistri(wimp_mass, detector, FWHM_rec, FWHM_ion, efficiency, 0, 0, 0, 6.4, 1)

c1 = TCanvas('c1','WIMP signal test',1000,500)
c1.Divide(2)
c1.cd(1)
signal.Draw('CONT0')
signal.SetContour(30)

ER_centroid.SetLineColor(kMagenta)
ER_centroid.SetLineWidth(2)
ER_centroid.SetParameter(0,6.4)
ER_centroid.DrawCopy('SAME')

NR_centroid.SetLineColor(kMagenta)
NR_centroid.SetLineWidth(2)
NR_centroid.DrawCopy('SAME')
NR_centroid.Draw('SAME')

#signal.GetXaxis().SetRangeUser(3.2,12.8)
signal.GetXaxis().SetTitle('E_{rec} (keV_{nr})')
#signal.GetYaxis().SetRangeUser(0.8,5.8)
signal.GetYaxis().SetTitle('E_{ion} (keV_{ee})')

c1.cd(2)
efficiency.Draw()
TriggerEfficiency.Draw('SAME')

print "integrated efficiency corrected signal rate [10^-6pb/kg/day]",signal.Integral('WIDTH')