#!/usr/bin/env python

####################################################################################################
##
## Calculation of 90% C.L. cross section limits using Maximum Likelihood methods based on RooFit
## Lukas Hehn, 2013
##
####################################################################################################

import ROOT
import functions
import parameters

from ROOT import RooFit as rf
from ROOT import TCanvas


# Global variables used in the script
DETECTOR_NAME = 'ID401'
LIVETIME = parameters.MEASURED_VALUES_ERIC[DETECTOR_NAME]['livetime']  # livetime after cuts in days
DETECTOR_MASS = 0.160  # detector mass in kg

#DATA_FILE = '/kalinka/home/hehn/PhD/LowMassEric/ID3_eventlist.txt'
DATA_FILE = '/kalinka/home/hehn/PhD/NewMLEanalysis/Data/{detector}_eventlist_ion-rec-only.txt'.format(detector=DETECTOR_NAME)


# Energy Binning and binsize
E_ION_MAX = 14.
E_REC_MAX = 25.
BINSIZE = 0.1


WIMP_MASS = 15
NUM_MC_SETS = False  # number of MC toy event sets: 0 means no MC study
(N_SIGNAL_MIN, N_SIGNAL_MAX) = (0.0, 30.0)  # lower border for n_signal parameter


# Flags to control procedures
E_ION_SCALING = True  # switches on fit parameter to scale energy ion
E_REC_SCALING = True  # switches on fit parameter to scale energy rec
BCKGD_MC = False  # switches on MC study where toy events are produces with bckgd model only
N_SIG_TEST = False  # either a False flag or a value for N_signal parameter before starting MC study
SAVE_PLOTS = True  # save all plots as png


# Detector specific parameters like resolutions
VOLTAGE = parameters.MEASURED_VALUES_LUKAS[DETECTOR_NAME]['voltage']
E_THRESH_EE = parameters.MEASURED_VALUES_LUKAS[DETECTOR_NAME]['threshold_ee']
E_THRESH_NR = parameters.MEASURED_VALUES_LUKAS[DETECTOR_NAME]['threshold_nr']
FWHM_HEAT = parameters.ENERGY_RESOLUTIONS_ERIC[DETECTOR_NAME]['Heat']  # 0.82=value@10keV for ID3 adopted from Eric
FWHM_ION = parameters.ENERGY_RESOLUTIONS_ERIC[DETECTOR_NAME]['Fiducial']  # 0.72=value@10keV for ID3 adopted from Eric
FWHM_REC = functions.fwhm_rec_from_heat(FWHM_HEAT, VOLTAGE, 10.)

print '\nVoltage={voltage}V, threshold={threshEE}keVee={threshNR}keVnr, FWHM_ion={ion}keVee, FWHM_heat={heat}keVee, FWHM_rec={rec:.2f}keVnr\n'.format(
    voltage=VOLTAGE, threshEE=E_THRESH_EE, threshNR=E_THRESH_NR, ion=FWHM_ION, heat=FWHM_HEAT, rec=FWHM_REC)


# Definition of maximum likelihood observables in RooFit
ION = ROOT.RooRealVar('ion', 'E_{ion}', 0., E_ION_MAX, 'keVee')
REC = ROOT.RooRealVar('rec', 'E_{rec}', 0., E_REC_MAX, 'keVnr')
SIGMA_REC = ROOT.RooConstVar('sigma_rec', 'recoil energy resolution', FWHM_REC/2.35)
SIGMA_ION = ROOT.RooConstVar('sigma_ion', 'ionization energy resolution', FWHM_ION/2.35)


# Calculation of specific detector efficiency and pdf
total_efficiency = functions.simple_efficiency(DETECTOR_NAME, E_THRESH_NR, FWHM_REC/2.35,
                                               binsize=BINSIZE, rec_max=E_REC_MAX, ion_max=E_ION_MAX
                                               )
total_efficiency_datahist = ROOT.RooDataHist('total_efficiency_datahist', 'total_efficiency_datahist',
                                             ROOT.RooArgList(REC, ION), total_efficiency)
total_efficiency_pdf = ROOT.RooHistPdf('total_efficiency_pdf', 'total_efficiency_pdf',
                                       ROOT.RooArgSet(REC, ION), total_efficiency_datahist)


# Read in of data set and output as scatter/graph
realdata = ROOT.RooDataSet.read(DATA_FILE, ROOT.RooArgList(REC, ION))
realdata_scatter = realdata.createHistogram(REC, ION, int(E_REC_MAX*1./BINSIZE), int(E_ION_MAX*1./BINSIZE))
realdata_graph = functions.tgraph_from_dataset(realdata)  # not significantly more accurate than binned scatter
events = int(realdata.numEntries())


# ER_centroid for calculation of peak position in Erec
functions.ER_CENTROID.FixParameter(0, VOLTAGE)


# Scaling factors allow for correction of gamma peak energy positions.
ion_scaling = ROOT.RooRealVar('ion_scaling', 'scaling factor ionization energy', 1.0, 0.95, 1.05)
rec_scaling = ROOT.RooRealVar('rec_scaling', 'scaling factor recoil energy', 1.0, 0.95, 1.05)
if not E_ION_SCALING:
    ion_scaling.setConstant(ROOT.kTRUE)
if not E_REC_SCALING:
    rec_scaling.setConstant(ROOT.kTRUE)


# Definition of gamma peaks with known Eion and calculated Erec.
V49_ion_energy = ROOT.RooConstVar('V49_ion_energy', 'V49 peak ion energy', 4.97)
V49_ion_pos = ROOT.RooFormulaVar('V49_ion_pos', '@0*@1', ROOT.RooArgList(V49_ion_energy, ion_scaling))
V49_ion_pdf = ROOT.RooGaussian('V49_ion_pdf', 'V49 peak gauss pdf in ion', ION, V49_ion_pos, SIGMA_ION)
V49_rec_energy = ROOT.RooConstVar('V49_rec_energy', 'recoil energy V49 peak', functions.ER_CENTROID.GetX(V49_ion_energy.getVal()))
V49_rec_pos = ROOT.RooFormulaVar('V49_rec_pos', '@0*@1', ROOT.RooArgList(V49_rec_energy, rec_scaling))
V49_rec_pdf = ROOT.RooGaussian('V49_rec_pdf', 'V49 peak pdf in recoil energy', REC, V49_rec_pos, SIGMA_REC)
V49_pdf = ROOT.RooProdPdf('V49_pdf', 'V49 peak pdf', V49_ion_pdf, V49_rec_pdf)
N_V49 = ROOT.RooRealVar('N_V49', 'evts of V49 peak (4.97keV)', 0., 0., events)
V49_ext = ROOT.RooExtendPdf('V49_ext', 'V49_ext', V49_pdf, N_V49)

Cr51_ion_energy = ROOT.RooConstVar('Cr51_ion_energy', 'Cr51 peak ion energy', 5.46)
Cr51_ion_pos = ROOT.RooFormulaVar('Cr51_ion_pos', '@0*@1', ROOT.RooArgList(Cr51_ion_energy, ion_scaling))
Cr51_ion = ROOT.RooGaussian('Cr51_ion_pdf', 'Cr51 peak gauss pdf with shifted mean', ION, Cr51_ion_pos, SIGMA_ION)
Cr51_rec_energy = ROOT.RooConstVar('Cr51_rec_energy', 'v_rec_energy', functions.ER_CENTROID.GetX(Cr51_ion_energy.getVal()))
Cr51_rec_pos = ROOT.RooFormulaVar('Cr51_rec_pos', '@0*@1', ROOT.RooArgList(Cr51_rec_energy, rec_scaling))
Cr51_rec = ROOT.RooGaussian('Cr51_rec_pdf', 'Cr51_rec_pdf with shifted mean', REC, Cr51_rec_pos, SIGMA_REC)
Cr51_pdf = ROOT.RooProdPdf('Cr51_pdf', 'Cr51 peak pdf', Cr51_ion, Cr51_rec)
N_Cr51 = ROOT.RooRealVar('N_Cr51', 'evts of 51Cr peak (5.46keV)', 0., 0., events)
Cr51_ext = ROOT.RooExtendPdf('Cr51_ext', 'Cr51_ext', Cr51_pdf, N_Cr51)

Mn54_ion_energy = ROOT.RooConstVar('Mn54_ion_energy', 'Mn54_ion_energy', 5.99)
Mn54_ion_pos = ROOT.RooFormulaVar('Mn54_ion_pos', '@0*@1', ROOT.RooArgList(Mn54_ion_energy, ion_scaling))
Mn54_ion = ROOT.RooGaussian('Mn54_ion_pdf', 'Mn54_ion_pdf with shifted mean', ION, Mn54_ion_pos, SIGMA_ION)
Mn54_rec_energy = ROOT.RooConstVar('Mn54_rec_energy', 'Mn54_rec_energy', functions.ER_CENTROID.GetX(Mn54_ion_energy.getVal()))
Mn54_rec_pos = ROOT.RooFormulaVar('Mn54_rec_pos', '@0*@1', ROOT.RooArgList(Mn54_rec_energy, rec_scaling))
Mn54_rec = ROOT.RooGaussian('Mn54_rec_pdf', 'Mn54_rec_pdf with shifted mean', REC, Mn54_rec_pos, SIGMA_REC)
Mn54_pdf = ROOT.RooProdPdf('Mn54_pdf', 'Mn54 peak pdf', Mn54_ion, Mn54_rec)
N_Mn54 = ROOT.RooRealVar('N_Mn54', 'evts of 54Mn peak (5.99keV)', 0., 0., events)
Mn54_ext = ROOT.RooExtendPdf('Mn54_ext', 'Mn54_ext', Mn54_pdf, N_Mn54)

Fe55_ion_energy = ROOT.RooConstVar('Fe55_ion_energy', 'Fe55_ion_energy', 6.54)
Fe55_ion_pos = ROOT.RooFormulaVar('Fe55_ion_pos', '@0*@1', ROOT.RooArgList(Fe55_ion_energy, ion_scaling))
Fe55_ion = ROOT.RooGaussian('Fe55_ion_pdf', 'Fe55_ion_pdf with shifted mean', ION, Fe55_ion_pos, SIGMA_ION)
Fe55_rec_energy = ROOT.RooConstVar('Fe55_rec_energy', 'Fe55_rec_energy', functions.ER_CENTROID.GetX(Fe55_ion_energy.getVal()))
Fe55_rec_pos = ROOT.RooFormulaVar('Fe55_rec_pos', '@0*@1', ROOT.RooArgList(Fe55_rec_energy, rec_scaling))
Fe55_rec = ROOT.RooGaussian('Fe55_rec_pdf', 'Fe55_rec_pdf with shifted mean', REC, Fe55_rec_pos, SIGMA_REC)
Fe55_pdf = ROOT.RooProdPdf('Fe55_pdf', 'Fe55 peak pdf', Fe55_ion, Fe55_rec)
N_Fe55 = ROOT.RooRealVar('N_Fe55', 'evts of 55Fe peak (6.54keV)', 0., 0., events)
Fe55_ext = ROOT.RooExtendPdf('Fe55_ext', 'Fe55_ext', Fe55_pdf, N_Fe55)

Co57_ion_energy = ROOT.RooConstVar('Co57_ion_energy', 'Co57_ion_energy', 7.11)
Co57_ion_pos = ROOT.RooFormulaVar('Co57_ion_pos', '@0*@1', ROOT.RooArgList(Co57_ion_energy, ion_scaling))
Co57_ion = ROOT.RooGaussian('Co57_ion_pdf', 'Co57_ion_pdf with shifted mean', ION, Co57_ion_pos, SIGMA_ION)
Co57_rec_energy = ROOT.RooConstVar('Co57_rec_energy', 'Co57_rec_energy', functions.ER_CENTROID.GetX(Co57_ion_energy.getVal()))
Co57_rec_pos = ROOT.RooFormulaVar('Co57_rec_pos', '@0*@1', ROOT.RooArgList(Co57_rec_energy, rec_scaling))
Co57_rec = ROOT.RooGaussian('Co57_rec_pdf', 'Co57_rec_pdf with shifted mean', REC, Co57_rec_pos, SIGMA_REC)
Co57_pdf = ROOT.RooProdPdf('Co57_pdf', 'Co57 peak pdf', Co57_ion, Co57_rec)
N_Co57 = ROOT.RooRealVar('N_Co57', 'evts of 57Co peak (7.11keV)', 0., 0., events)
Co57_ext = ROOT.RooExtendPdf('Co57_ext', 'Co57_ext', Co57_pdf, N_Co57)

Zn65_ion_energy = ROOT.RooConstVar('Zn65_ion_energy', 'Zn65_ion_energy', 8.98)
Zn65_ion_pos = ROOT.RooFormulaVar('Zn65_ion_pos', '@0*@1', ROOT.RooArgList(Zn65_ion_energy, ion_scaling))
Zn65_ion = ROOT.RooGaussian('Zn65_ion_pdf', 'Zn65_ion_pdf with shifted mean', ION, Zn65_ion_pos, SIGMA_ION)
Zn65_rec_energy = ROOT.RooConstVar('Zn65_rec_energy', 'Zn65_rec_energy', functions.ER_CENTROID.GetX(Zn65_ion_energy.getVal()))
Zn65_rec_pos = ROOT.RooFormulaVar('Zn65_rec_pos', '@0*@1', ROOT.RooArgList(Zn65_rec_energy, rec_scaling))
Zn65_rec = ROOT.RooGaussian('Zn65_rec_pdf', 'Zn65_rec_pdf with shifted mean', REC, Zn65_rec_pos, SIGMA_REC)
Zn65_pdf = ROOT.RooProdPdf('Zn65_pdf', 'Zn65 peak pdf', Zn65_ion, Zn65_rec)
N_Zn65 = ROOT.RooRealVar('N_Zn65', 'evts of 65Zn peak (8.98keV)', 0., 0., events)
Zn65_ext = ROOT.RooExtendPdf('Zn65_ext', 'Zn65_ext', Zn65_pdf, N_Zn65)

Ga68_ion_energy = ROOT.RooConstVar('Ga68_ion_energy', 'Ga68_ion_energy', 9.66)
Ga68_ion_pos = ROOT.RooFormulaVar('Ga68_ion_pos', '@0*@1', ROOT.RooArgList(Ga68_ion_energy, ion_scaling))
Ga68_ion = ROOT.RooGaussian('Ga68_ion_pdf', 'Ga68_ion_pdf with shifted mean', ION, Ga68_ion_pos, SIGMA_ION)
Ga68_rec_energy = ROOT.RooConstVar('Ga68_rec_energy', 'Ga68_rec_energy', functions.ER_CENTROID.GetX(Ga68_ion_energy.getVal()))
Ga68_rec_pos = ROOT.RooFormulaVar('Ga68_rec_pos', '@0*@1', ROOT.RooArgList(Ga68_rec_energy, rec_scaling))
Ga68_rec = ROOT.RooGaussian('Ga68_rec_pdf', 'Ga68_rec_pdf with shifted mean', REC, Ga68_rec_pos, SIGMA_REC)
Ga68_pdf = ROOT.RooProdPdf('Ga68_pdf', 'Ga68 peak pdf', Ga68_ion, Ga68_rec)
N_Ga68 = ROOT.RooRealVar('N_Ga68', 'evts of 68Ga peak (9.66keV)', 0., 0., events)
Ga68_ext = ROOT.RooExtendPdf('Ga68_ext', 'Ga68_ext', Ga68_pdf, N_Ga68)

Ge68_ion_energy = ROOT.RooConstVar('Ge68_ion_energy', 'Ge68_ion_energy', 10.37)
Ge68_ion_pos = ROOT.RooFormulaVar('Ge68_ion_pos', '@0*@1', ROOT.RooArgList(Ge68_ion_energy, ion_scaling))
Ge68_ion = ROOT.RooGaussian('Ge68_ion_pdf', 'Ge68_ion_pdf with shifted mean', ION, Ge68_ion_pos, SIGMA_ION)
Ge68_rec_energy = ROOT.RooConstVar('Ge68_rec_energy', 'Ge68_rec_energy', functions.ER_CENTROID.GetX(Ge68_ion_energy.getVal()))
Ge68_rec_pos = ROOT.RooFormulaVar('Ge68_rec_pos', '@0*@1', ROOT.RooArgList(Ge68_rec_energy, rec_scaling))
Ge68_rec = ROOT.RooGaussian('Ge68_rec_pdf', 'Ge68 peak pdf in rec', REC, Ge68_rec_pos, SIGMA_REC)
Ge68_pdf = ROOT.RooProdPdf('Ge68_pdf', 'Ge68 peak pdf', Ge68_ion, Ge68_rec)
N_Ge68 = ROOT.RooRealVar('N_Ge68', 'evts of 68Ge peak (10.37keV)', 0., 0., events)
Ge68_ext = ROOT.RooExtendPdf('Ge68_ext', 'Ge68_ext', Ge68_pdf, N_Ge68)


# Definition of WIMP signal and pdf
signal_hist = functions.wimp_signal(WIMP_MASS, SIGMA_ION.getVal(), SIGMA_REC.getVal(),
                                    binsize=BINSIZE, rec_max=E_REC_MAX, ion_max=E_ION_MAX
                                    )
signal_hist.Multiply(total_efficiency)
signal_datahist = ROOT.RooDataHist('signal_datahist', 'signal_datahist',
                                   ROOT.RooArgList(REC, ION), signal_hist)
signal_pdf = ROOT.RooHistPdf('signal_pdf', 'signal_pdf',
                             ROOT.RooArgSet(REC, ION), signal_datahist)
N_signal = ROOT.RooRealVar('N_signal', 'WIMP signal events', 0., N_SIGNAL_MIN, N_SIGNAL_MAX)
sig_ext = ROOT.RooExtendPdf('sig_ext', 'sig_ext', signal_pdf, N_signal)


# Definition of flat gamma background pdf
flat_gamma_bckgd_hist = functions.flat_gamma_bckgd(SIGMA_ION.getVal(), SIGMA_REC.getVal(), VOLTAGE,
                                                   binsize=BINSIZE, rec_max=E_REC_MAX, ion_max=E_ION_MAX
                                                   )
flat_gamma_bckgd_hist.Multiply(total_efficiency)
flat_gamma_bckgd_datahist = ROOT.RooDataHist('flat_gamma_bckgd_datahist', 'flat_gamma_bckgd_datahist',
                                             ROOT.RooArgList(REC, ION), flat_gamma_bckgd_hist)
flat_gamma_bckgd_pdf = ROOT.RooHistPdf('flat_gamma_bckgd_pdf', 'flat_gamma_bckgd_pdf',
                                       ROOT.RooArgSet(REC, ION), flat_gamma_bckgd_datahist)
N_flat = ROOT.RooRealVar('N_flat', 'bckgd events', 0., 0., events)
flat_ext = ROOT.RooExtendPdf('flat_ext', 'flat_ext', flat_gamma_bckgd_pdf, N_flat)


# Definition of background only as well as background plus signal pdf
bckgd_and_sig_pdf = ROOT.RooAddPdf('bckgd_and_sig_pdf', 'bckgd_and_sig_pdf',
                                   ROOT.RooArgList(flat_ext, V49_ext, Cr51_ext, Mn54_ext, Fe55_ext, Zn65_ext, Ga68_ext, Ge68_ext, sig_ext))  #  , Co57_ext
bckgd_only_pdf = ROOT.RooAddPdf('bckgd_only_pdf', 'bckgd_only_pdf',
                                ROOT.RooArgList(flat_ext, V49_ext, Cr51_ext, Mn54_ext, Fe55_ext, Zn65_ext, Ga68_ext, Ge68_ext))  #  , Co57_ext


# Create negative log likelihood (NLL) object manually and minimize it
nll = ROOT.RooNLLVar('nll', 'nll', bckgd_and_sig_pdf, realdata, rf.Extended(ROOT.kTRUE),
                     rf.PrintEvalErrors(0), rf.Verbose(ROOT.kFALSE), rf.NumCPU(2))
minuit = ROOT.RooMinuit(nll)
minuit.migrad()  # find minimum
minuit.hesse()  # find symmetric errors
minuit.minos()  # find Asymmetric errors
FitResult = minuit.save('realfit', 'fit to real data')
FitParams = FitResult.floatParsFinal()
FitResult.Print('v')


# Create histogram of best fit PDF
bckgd_and_sig_hist = bckgd_and_sig_pdf.createHistogram('bckgd_and_sig_hist', REC, rf.Binning(int(E_REC_MAX*1/BINSIZE)),
                                                       rf.YVar(ION, rf.Binning(int(E_ION_MAX*1/BINSIZE)))
                                                       )
bckgd_and_sig_hist.SetTitle('{detector} best fit PDF Bckgd + {mass}GeV Signal'.format(detector=DETECTOR_NAME, mass=WIMP_MASS))
bckgd_and_sig_hist.SetStats(0)


# Calculate cross section limit from fit results and output of results
signal_parameter = FitResult.floatParsFinal().find('N_signal')
N_sig = signal_parameter.getVal()
sigma_n_sig_low = signal_parameter.getAsymErrorLo()
sigma_n_sig_high = signal_parameter.getAsymErrorHi()
N_max = N_sig + 1.28 * sigma_n_sig_high
rate = signal_hist.Integral('WIDTH')  # option width absolutely necessary
N_wimp = rate * LIVETIME * DETECTOR_MASS * 1e6  # 1e6 factor to scale to pb
xs_max = N_max / N_wimp
result_overview = {'N_sig' : N_sig, 'N_max' : N_max, 'rate' : rate, 'xs_max' : xs_max, 'N_wimp' : N_wimp}

print '{0:9} | {1:10} | {2:8} | {3:8} | {4:10} | {5:10}'.format('WIMP_MASS', 'Rate', 'N_signal',
                                                                'ErrorLow', 'ErrorHigh', 'XS-limit [pb]')
print '-'*80
print '{0:9} | {1:10} | {2:8} | {3:8} | {4:10} | {5:10} ({5:.1e})'.format(WIMP_MASS, rate, N_sig,
                                                                sigma_n_sig_low, sigma_n_sig_high, xs_max)


# Definition of RooFit projections over both energies starting with appropriate binning
recbins = int(E_REC_MAX*2)
ionbins = int(E_ION_MAX*4)

ionframe = ION.frame()
ionframe.SetTitle('PDF component projection in E_{ion}')
realdata.plotOn(ionframe, rf.Name('data'),
                rf.Binning(ionbins), rf.MarkerSize(1.0))
bckgd_and_sig_pdf.plotOn(ionframe, rf.Components("flat_ext"), rf.Normalization(1.0,ROOT.RooAbsReal.RelativeExpected),
                         rf.LineColor(ROOT.kGreen), rf.LineWidth(2), rf.LineStyle(ROOT.kSolid))
bckgd_and_sig_pdf.plotOn(ionframe, rf.Components("V49_ext"), rf.Normalization(1.0,ROOT.RooAbsReal.RelativeExpected),
                         rf.LineColor(ROOT.kRed), rf.LineWidth(2), rf.LineStyle(ROOT.kDashed))
bckgd_and_sig_pdf.plotOn(ionframe, rf.Components("Cr51_ext"), rf.Normalization(1.0,ROOT.RooAbsReal.RelativeExpected),
                         rf.LineColor(ROOT.kRed), rf.LineWidth(2), rf.LineStyle(ROOT.kDashed))
bckgd_and_sig_pdf.plotOn(ionframe, rf.Components("Mn54_ext"), rf.Normalization(1.0,ROOT.RooAbsReal.RelativeExpected),
                         rf.LineColor(ROOT.kRed), rf.LineWidth(2), rf.LineStyle(ROOT.kDashed))
bckgd_and_sig_pdf.plotOn(ionframe, rf.Components("Fe55_ext"), rf.Normalization(1.0,ROOT.RooAbsReal.RelativeExpected),
                         rf.LineColor(ROOT.kRed), rf.LineWidth(2), rf.LineStyle(ROOT.kDashed))
#bckgd_and_sig_pdf.plotOn(ionframe, rf.Components("Co57_ext"), rf.Normalization(1.0,ROOT.RooAbsReal.RelativeExpected),
                 #rf.LineColor(ROOT.kRed), rf.LineWidth(2), rf.LineStyle(ROOT.kDashed))
bckgd_and_sig_pdf.plotOn(ionframe, rf.Components("Zn65_ext"), rf.Normalization(1.0,ROOT.RooAbsReal.RelativeExpected),
                         rf.LineColor(ROOT.kRed), rf.LineWidth(2), rf.LineStyle(ROOT.kDashed))
bckgd_and_sig_pdf.plotOn(ionframe, rf.Components("Ge68_ext"), rf.Normalization(1.0,ROOT.RooAbsReal.RelativeExpected),
                         rf.LineColor(ROOT.kRed), rf.LineWidth(2), rf.LineStyle(ROOT.kDashed))
bckgd_and_sig_pdf.plotOn(ionframe, rf.Components("Ga68_ext"), rf.Normalization(1.0,ROOT.RooAbsReal.RelativeExpected),
                         rf.LineColor(ROOT.kRed), rf.LineWidth(2), rf.LineStyle(ROOT.kDashed))
bckgd_and_sig_pdf.plotOn(ionframe, rf.Components("sig_ext"), rf.Normalization(1.0,ROOT.RooAbsReal.RelativeExpected),
                         rf.LineColor(ROOT.kMagenta), rf.LineWidth(3), rf.LineStyle(ROOT.kSolid))  # normal signal
bckgd_and_sig_pdf.plotOn(ionframe, rf.Components("sig_ext"), rf.Normalization( (N_max/N_sig), ROOT.RooAbsReal.RelativeExpected),
                         rf.LineColor(ROOT.kMagenta), rf.LineWidth(3), rf.LineStyle(ROOT.kSolid))  # N_max
bckgd_and_sig_pdf.plotOn(ionframe, rf.Normalization(1.0,ROOT.RooAbsReal.RelativeExpected),
                         rf.LineColor(ROOT.kBlue), rf.LineWidth(2), rf.LineStyle(ROOT.kSolid))


recframe = REC.frame()
recframe.SetTitle('PDF component projection in E_{rec}')
realdata.plotOn(recframe, rf.Name("data"),
                rf.Binning(recbins), rf.MarkerColor(ROOT.kBlack), rf.MarkerSize(1.0))
bckgd_and_sig_pdf.plotOn(recframe, rf.Components("flat_ext"), rf.Normalization(1.0,ROOT.RooAbsReal.RelativeExpected),
                         rf.LineColor(ROOT.kGreen), rf.LineWidth(2), rf.LineStyle(ROOT.kSolid))
bckgd_and_sig_pdf.plotOn(recframe, rf.Components("V49_ext"), rf.Normalization(1.0,ROOT.RooAbsReal.RelativeExpected),
                         rf.LineColor(ROOT.kRed), rf.LineWidth(2), rf.LineStyle(ROOT.kDashed))
bckgd_and_sig_pdf.plotOn(recframe, rf.Components("Cr51_ext"), rf.Normalization(1.0,ROOT.RooAbsReal.RelativeExpected),
                         rf.LineColor(ROOT.kRed), rf.LineWidth(2), rf.LineStyle(ROOT.kDashed))
bckgd_and_sig_pdf.plotOn(recframe, rf.Components("Mn54_ext"), rf.Normalization(1.0,ROOT.RooAbsReal.RelativeExpected),
                         rf.LineColor(ROOT.kRed), rf.LineWidth(2), rf.LineStyle(ROOT.kDashed))
bckgd_and_sig_pdf.plotOn(recframe, rf.Components("Fe55_ext"), rf.Normalization(1.0,ROOT.RooAbsReal.RelativeExpected),
                         rf.LineColor(ROOT.kRed), rf.LineWidth(2), rf.LineStyle(ROOT.kDashed))
#bckgd_and_sig_pdf.plotOn(recframe, rf.Components("Co57_ext"), rf.Normalization(1.0,ROOT.RooAbsReal.RelativeExpected),
                         #rf.LineColor(ROOT.kRed), rf.LineWidth(2), rf.LineStyle(ROOT.kDashed))
bckgd_and_sig_pdf.plotOn(recframe, rf.Components("Zn65_ext"), rf.Normalization(1.0,ROOT.RooAbsReal.RelativeExpected),
                         rf.LineColor(ROOT.kRed), rf.LineWidth(2), rf.LineStyle(ROOT.kDashed))
bckgd_and_sig_pdf.plotOn(recframe, rf.Components("Ge68_ext"), rf.Normalization(1.0,ROOT.RooAbsReal.RelativeExpected),
                         rf.LineColor(ROOT.kRed), rf.LineWidth(2), rf.LineStyle(ROOT.kDashed))
bckgd_and_sig_pdf.plotOn(recframe, rf.Components("Ga68_ext"), rf.Normalization(1.0,ROOT.RooAbsReal.RelativeExpected),
                         rf.LineColor(ROOT.kRed), rf.LineWidth(2), rf.LineStyle(ROOT.kDashed))
bckgd_and_sig_pdf.plotOn(recframe, rf.Components("sig_ext"), rf.Normalization(1.0,ROOT.RooAbsReal.RelativeExpected),
                         rf.LineColor(ROOT.kMagenta), rf.LineWidth(3), rf.LineStyle(ROOT.kSolid))
bckgd_and_sig_pdf.plotOn(recframe, rf.Components("sig_ext"), rf.Normalization( (N_max/N_sig), ROOT.RooAbsReal.RelativeExpected),
                         rf.LineColor(ROOT.kMagenta), rf.LineWidth(3), rf.LineStyle(ROOT.kSolid))  # N_max
bckgd_and_sig_pdf.plotOn(recframe, rf.Normalization(1.0,ROOT.RooAbsReal.RelativeExpected),
                         rf.LineColor(ROOT.kBlue), rf.LineWidth(2), rf.LineStyle(ROOT.kSolid))
# Additional box with parameter fit values
bckgd_and_sig_pdf.paramOn(recframe, rf.Format('NEU', rf.AutoPrecision(2)),
                          rf.Layout(0.1, 0.5, 0.9), rf.ShowConstants(ROOT.kFALSE))


# NLL frame for N_signal parameter
NsigFit = FitResult.floatParsFinal().find('N_signal')
signalNLLframe = NsigFit.frame()
signalNLLframe.SetTitle('Negative Log-Likelihood for N_signal')
nll.plotOn(signalNLLframe, rf.Precision(1e-5), rf.ShiftToZero())
signalNLLframe.GetXaxis().SetTitle('N_signal')


# Definition of some vertical lines which are used to mark parameter fit-value positions
NllLine = ROOT.TLine()
NllLine.SetLineWidth(2)
NllLine.SetLineStyle(1)
NllLine.SetLineColor(ROOT.kBlue)

ZeroLine = ROOT.TLine()
ZeroLine.SetLineWidth(2)
ZeroLine.SetLineStyle(7)
ZeroLine.SetLineColor(ROOT.kBlack)

ParamLine = ROOT.TLine()
ParamLine.SetLineWidth(2)
ParamLine.SetLineStyle(1)
ParamLine.SetLineColor(ROOT.kMagenta)

ERline = functions.ER_CENTROID
ERline.SetLineColor(ROOT.kBlack)
ERline.SetLineWidth(2)
NRline = functions.NR_CENTROID
NRline.SetLineColor(ROOT.kBlack)
NRline.SetLineWidth(2)


# Plotting of fit results (best fit pdf + data, projections in both energies, NLL function)
c1 = TCanvas('c1', 'Fit result overview for %i GeV'%WIMP_MASS, 1200, 900)
c1.Divide(2, 2)

c1.cd(1)
signal_hist.Draw('CONT LIST')
signal_hist.SetContour(99)
c1.Update()
contours = ROOT.gROOT.GetListOfSpecials().FindObject("contours")
lowcontour = contours.At(10)
lowdummy = lowcontour.First()
lowlevel = lowdummy.Clone()
lowlevel.SetLineColor(ROOT.kMagenta)
lowlevel.SetLineWidth(3)


midcontour = contours.At(50)
middummy = midcontour.First()
midlevel = middummy.Clone()
midlevel.SetLineColor(ROOT.kMagenta)
midlevel.SetLineWidth(3)


highcontour = contours.At(90)
highdummy = highcontour.First()
highlevel = highdummy.Clone()
highlevel.SetLineColor(ROOT.kMagenta)
highlevel.SetLineWidth(3)

c1.cd(1).SetLogz()
bckgd_and_sig_hist.Draw('COLZ')
realdata_graph.SetMarkerColor(ROOT.kBlack)
realdata_graph.SetMarkerStyle(ROOT.kPlus)
realdata_graph.Draw('SAMESP')
ERline.DrawCopy('SAME')
NRline.DrawCopy('SAME')
lowlevel.Draw("SAMEL")
midlevel.Draw("SAMEL")
highlevel.Draw("SAMEL")
c1.cd(3)
signalNLLframe.Draw()
signalNLLframe.SetMinimum(0.)
signalNLLframe.SetMaximum(6.)
ROOT.gPad.Update()
ParamLine.DrawLine(NsigFit.getVal(), ROOT.gPad.GetUymin(), NsigFit.getVal(), ROOT.gPad.GetUymax())
ParamLine.SetLineWidth(3)

c1.cd(2)
ionframe.Draw()

c1.cd(4)
recframe.Draw()


# Set N_signal to different value than best fit to test resulting MC distribution
if not BCKGD_MC and N_SIG_TEST:
    N_signal.setVal(N_SIG_TEST)
    N_signal.setError(0.0)


# Monte Carlo toy event sets and output
if NUM_MC_SETS:
    # Creation of Monte Carlo toy event sets from bckgd only for theoretical limits
    if BCKGD_MC:
        MC_study = ROOT.RooMCStudy(bckgd_only_pdf, ROOT.RooArgSet(REC, ION), rf.Silence(),
                                   rf.Extended(ROOT.kTRUE), rf.FitOptions(rf.Save(ROOT.kTRUE)),
                                   rf.FitModel(bckgd_and_sig_pdf)
                                   )
    else:
        MC_study = ROOT.RooMCStudy(bckgd_and_sig_pdf, ROOT.RooArgSet(REC, ION), rf.Silence(),
                                   rf.Extended(ROOT.kTRUE), rf.FitOptions(rf.Save(ROOT.kTRUE)),
                                   )
    MC_study.generateAndFit(NUM_MC_SETS)


    # plot results of the Monte Carlo toy set fit for all parameters
    c2 = TCanvas('c2', 'MC toy set statistics for {mass} GeV'.format(mass=WIMP_MASS))
    NumFitParams = FitParams.getSize()
    c2.Divide(1, NumFitParams, 0.001, 0.001)

    ParamLine.SetLineColor(ROOT.kRed)

    for i in range(NumFitParams):
        parameter = FitParams[i]
        paramname = parameter.GetName()
        paramvalue = parameter.getVal()

        pad = c2.cd(i+1)
        pad.Divide(3, 1, 0.001, 0.001)

        pad.cd(1)
        ParamNLLFrame = parameter.frame()
        #ParamNLLFrame.SetTitle('NLL fit {name}'.format(name=paramname))
        nll.plotOn(ParamNLLFrame, rf.Precision(1e-5), rf.ShiftToZero())
        ParamNLLFrame.SetMinimum(0.)
        ParamNLLFrame.SetMaximum(100.)
        ParamNLLFrame.Draw()
        ROOT.gPad.Update()
        ParamLine.DrawLine(paramvalue, ROOT.gPad.GetUymin(), paramvalue, ROOT.gPad.GetUymax())

        pad.cd(2)
        ParamDistriFrame = MC_study.plotParam(parameter)
        #ParamDistriFrame.SetTitle('MC distri {name}'.format(name=paramname))
        ParamDistriFrame.Draw()
        #ParamDistriFrame.getHist().Fit('gaus', 'QEM')
        ROOT.gPad.Update()
        ParamLine.DrawLine(paramvalue, ROOT.gPad.GetUymin(), paramvalue, ROOT.gPad.GetUymax())

        #pad.cd(3)
        #ParamPullFrame = MC_study.plotPull(parameter, rf.FitGauss())
        ##ParamPullFrame.SetTitle('MC pull distri {name}'.format(name=paramname))
        #ParamPullFrame.Draw()
        #ROOT.gPad.Update()
        #ZeroLine.DrawLine(0, ROOT.gPad.GetUymin(), 0, ROOT.gPad.GetUymax())

    c2.SetCanvasSize(1200, 2400)


    # Additional canvas with distribution of NLL values for all toy set fits
    RealDataBestFitNLL = nll.getVal()
    c3 = TCanvas('c3', 'NLL distribution for Monte Carlo toy sets {mass}GeV'.format(mass=WIMP_MASS), 800, 600)
    MCnllframe = MC_study.plotNLL()
    MCnllframe.SetTitle('NLL distribution for {mass}GeV'.format(mass=WIMP_MASS))
    MCnllframe.Draw()
    MCnllframe.getHist().Fit('gaus', 'QEM')
    ROOT.gPad.Update()
    NllLine.DrawLine(RealDataBestFitNLL, ROOT.gPad.GetUymin(), RealDataBestFitNLL, ROOT.gPad.GetUymax())


# Canvas for signal parameters only (including MC toy event set fits)
if NUM_MC_SETS:
    signal = FitResult.floatParsFinal().find('N_signal')
    signalvalue = signal.getVal()

    c4 = TCanvas('c4', 'N_signal fit & MC toy set statistics for {mass} GeV'.format(mass=WIMP_MASS), 1000, 500)
    c4.Divide(2)

    c4.cd(1)
    ParamDistriFrame = MC_study.plotParam(signal, rf.Binning(400))
    ParamDistriFrame.SetTitle('MC toy set fits')
    ParamDistriFrame.Draw()
    ParamDistriFrame.getHist().Fit('gaus', 'QEM')
    ROOT.gPad.Update()
    ParamLine.SetLineColor(ROOT.kMagenta)
    ParamLine.DrawLine(signalvalue, ROOT.gPad.GetUymin(), signalvalue, ROOT.gPad.GetUymax())

    c4.cd(2)
    ParamPullFrame = MC_study.plotPull(signal)  #rf.FitGauss(ROOT.kTRUE) not really useful because of non-gaussian features
    ParamPullFrame.SetTitle('MC toy set fits pull')
    ParamPullFrame.Draw()
    ROOT.gPad.Update()
    ZeroLine.DrawLine(0, ROOT.gPad.GetUymin(), 0, ROOT.gPad.GetUymax())


# Alternative calculation of 90% C.L. limit on cross section using Monte Carlo statistics
if NUM_MC_SETS:
    N_sig_list = []
    for i in range(NUM_MC_SETS):
        value = MC_study.fitParams(i).find('N_signal').getVal()
        N_sig_list.append(value)
    N_sig_list.sort()
    result_overview['N_sig_MC'] = N_sig_list[int(0.9*NUM_MC_SETS)]
    result_overview['xs_max_MC'] = result_overview['N_sig_MC']/result_overview['N_wimp']
    print "\n2 different methods for limit estimation:"
    print '{0:15} | {1:15} | {2:10}'.format('method', '90% C.L. N_max','90% C.L. xs')
    print '-'*80
    print '{0:15} | {1:15.3f} | {2:15.2e}'.format('NLL curve',
                                                  result_overview['N_max'],
                                                  result_overview['xs_max'])
    print '{0:15} | {1:15.3f} | {2:15.2e}'.format('MC toy set fits',
                                                  result_overview['N_sig_MC'],
                                                  result_overview['xs_max_MC'])


if SAVE_PLOTS:
    c1.SaveAs('%iGeV_PDFs.png'%WIMP_MASS)
    if NUM_MC_SETS:
        c2.SaveAs('%iGeV_param-stats.png'%WIMP_MASS)
        c3.SaveAs('%iGeV_signal-stats.png'%WIMP_MASS)


# Print overview of fit results for all MC toy sets
if NUM_MC_SETS:
    goodlist, badlist = [], []
    for i in range(NUM_MC_SETS):
        result = MC_study.fitResult(i)
        n_sig = result.floatParsFinal().find('N_signal')
        edm = result.edm()
        covqual = result.covQual()
        status = result.status()
        invalid_nll = result.numInvalidNLL()
        n_sig_error_low = n_sig.getErrorLo()
        n_sig_val = n_sig.getVal()
        print '{i:3<} {n_sig_val:8.2f} {n_sig_error_low:8.2f} {edm:10.2e} {status} {covqual} {invalid_nll}'.format(i=i, edm=edm, covqual=covqual, invalid_nll=invalid_nll, n_sig_error_low=n_sig_error_low, n_sig_val=n_sig_val, status=status)
        if invalid_nll > 0:
            badlist.append(n_sig_val)
        else:
            goodlist.append(n_sig_val)
