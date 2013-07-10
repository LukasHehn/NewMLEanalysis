from ROOT import *
from Functions import *

#gROOT.LoadMacro("/kalinka/home/hehn/PhD/LowMassEric/WimpDistriExtended.C")
gROOT.LoadMacro("/kalinka/home/hehn/PhD/LowMassEric/WimpDistriAdapted.C")

FWHM_heat = 0.82
FWHM_rec = RecoilResolutionFromHeat(FWHM_heat,6.4,10)
sigma_rec = FWHM_rec/2.35
FWHM_ion = 0.72

wimp_mass = 10
flag2cut = 1
ioncut = 0
reccut = 0
detector = 'ID3'
voltage = 6.4

TriggerEfficiency.SetParameter(0, 3.874)
TriggerEfficiency.SetParameter(1, sigma_rec)

TriggerEfficiency.SetNpx(200)
TriggerEfficiencyHist = TriggerEfficiency.GetHistogram()

spectrum = ReadInWimpSpectrumEric(str(wimp_mass))

gammacut = 0
signal = WimpDistri(str(wimp_mass), detector, FWHM_rec, FWHM_ion, TriggerEfficiencyHist, ioncut, reccut, gammacut, voltage, flag2cut)

gammacut = 1
signalcut = WimpDistri(str(wimp_mass), detector, FWHM_rec, FWHM_ion, TriggerEfficiencyHist, ioncut, reccut, gammacut, voltage, flag2cut)

c1 = TCanvas('c1','WIMP signal test',1000,500)
c1.Divide(2,1,0.01,0.01)
c1.cd(1)
spectrum.SetTitle('WIMP Spectrum %i GeV'%wimp_mass)
spectrum.SetStats(0)
spectrum.GetXaxis().SetTitleSize(0.05)
spectrum.GetXaxis().SetLabelSize(0.05)
spectrum.GetXaxis().SetTitleOffset(1.2)
spectrum.GetYaxis().SetTitleSize(0.05)
spectrum.GetYaxis().SetLabelSize(0.05)
spectrum.GetYaxis().SetTitleOffset(1.0)
spectrum.Draw()
gPad.Update()

TriggerEfficiencyHist.Scale(gPad.GetUymax())
rightmax = 1.0
TriggerEfficiencyHist.Draw('SAME')

axis1 = TGaxis(gPad.GetUxmax(),gPad.GetUymin(), gPad.GetUxmax(), gPad.GetUymax(),0,1.0,510,"+L")
axis1.SetLineColor(kRed)
axis1.SetLabelColor(kRed)
axis1.SetTitleColor(kRed)
axis1.SetTitle('Trigger Efficiency')
axis1.SetTitleSize(0.05)
axis1.SetLabelSize(0.05)
axis1.SetTitleOffset(0.7)
axis1.Draw()

c1.cd(2)
signal.Draw('CONT0Z')
signal.SetContour(50)

ER_centroid.SetLineColor(kBlack)
ER_centroid.SetLineWidth(2)
ER_centroid.SetParameter(0,6.4)
ER_centroid.Draw('SAME')

NR_centroid.SetLineColor(kBlack)
NR_centroid.SetLineWidth(2)
NR_centroid.Draw('SAME')
#NR_centroid.Draw('SAME')

#signal.GetXaxis().SetRangeUser(3.2,12.8)
#signal.GetYaxis().SetRangeUser(0.8,5.8)

signal.SetTitle('WIMP Signal %i GeV'%wimp_mass)
signal.SetStats(0)
signal.GetXaxis().SetTitle('E_{recoil} [keV_{nr}]')
signal.GetYaxis().SetTitle('E_{ion} [keV_{ee}]')


print 'wimp mass: {0:} GeV'.format(wimp_mass)
print 'integral 1D rate:  {0:.2f} [10^-6pb/kg/day]'.format(spectrum.Integral('WIDTH'))
print 'integral 2D rate:  {0:.2e} [10^-6pb/kg/day]'.format(signal.Integral('WIDTH'))
print 'integral 2D rate (with gamma cut):  {0:.2e} [10^-6pb/kg/day]'.format(signalcut.Integral('WIDTH'))
