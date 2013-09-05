#!/usr/bin/env python

######################################################
#
# Calculation of 90% C.L. cross section limits using Maximum Likelihood methods based on RooFit
# Lukas Hehn, 2013
#
######################################################

import ROOT
import functions
from ROOT import RooFit as rf



# Global variables used in the script
DETECTOR_NAME = 'ID3'
DATA_FILE = '/kalinka/home/hehn/PhD/LowMassEric/ID3_eventlist.txt'
E_ION_MAX = 10.
E_REC_MAX = 20.
WIMP_MASS = 8
NUM_MC_SETS = int(1e3)  # number of MC toy event sets: 0 means no MC study
SAVE_PLOTS = True


# Detector specific parameters like resolutions
VOLTAGE = 6.4
E_THRESH = 3.874
FWHM_HEAT = 0.82  # value for ID3 adopted from Eric for 10keV
FWHM_ION = 0.72  # value for ID3 adopted from Eric for 10keV
FWHM_REC = functions.fwhm_rec_from_heat(FWHM_HEAT, VOLTAGE, 10.)


# Definition of maximum likelihood observables in RooFit
ION = ROOT.RooRealVar('ion', 'E_{ion}', 0., E_ION_MAX, 'keV_{ee}')
REC = ROOT.RooRealVar('rec', 'E_{rec}', 0., E_REC_MAX, 'keV_{nr}')
SIGMA_REC = ROOT.RooRealVar('sigma_rec', 'recoil energy resolution', FWHM_REC/2.35)
SIGMA_ION = ROOT.RooRealVar('sigma_ion', 'ionization energy resolution', FWHM_ION/2.35)


# Calculation of specific detector efficiency and pdf
total_efficiency = functions.simple_efficiency(DETECTOR_NAME, E_THRESH, FWHM_REC/2.35)
total_efficiency_datahist = ROOT.RooDataHist('total_efficiency_datahist', 'total_efficiency_datahist', 
                                             ROOT.RooArgList(REC, ION), total_efficiency)
total_efficiency_pdf = ROOT.RooHistPdf('total_efficiency_pdf', 'total_efficiency_pdf', 
                                       ROOT.RooArgSet(REC, ION), total_efficiency_datahist)


# Read in of data set and output as scatter/graph
realdata = ROOT.RooDataSet.read(DATA_FILE, ROOT.RooArgList(REC, ION))
realdata_scatter = realdata.createHistogram(REC, ION, int(E_REC_MAX*10), int(E_ION_MAX*10))
#realdata_graph = functions.tgraph_from_dataset(realdata)  # tgraph not really more accurate than binned scatter
events = int(realdata.numEntries())


# Definition of gamma peaks
ion_scaling = ROOT.RooRealVar('ion_scaling', 'scaling factor ionization energy', 1.0, 0.95, 1.05)
rec_scaling = ROOT.RooRealVar('rec_scaling', 'scaling factor recoil energy', 1.0, 0.95, 1.05)
ion_scaling.setConstant(ROOT.kTRUE)
rec_scaling.setConstant(ROOT.kTRUE)
functions.ER_CENTROID.FixParameter(0, VOLTAGE)  # ER_centroid for calculation of peak position in Erec

V49_ion_energy = ROOT.RooRealVar('V49_ion_energy', 'V49 peak ion energy', 4.97)
V49_ion_pos = ROOT.RooFormulaVar('V49_ion_pos', '@0*@1', ROOT.RooArgList(V49_ion_energy, ion_scaling))
V49_ion_pdf = ROOT.RooGaussian('V49_ion_pdf', 'V49 peak gauss pdf in ion', ION, V49_ion_pos, SIGMA_ION)
V49_rec_energy = ROOT.RooRealVar('V49_rec_energy', 'recoil energy V49 peak', functions.ER_CENTROID.GetX(V49_ion_energy.getVal()))
V49_rec_pos = ROOT.RooFormulaVar('V49_rec_pos', '@0*@1', ROOT.RooArgList(V49_rec_energy, rec_scaling))
V49_rec_pdf = ROOT.RooGaussian('V49_rec_pdf', 'V49 peak pdf in recoil energy', REC, V49_rec_pos, SIGMA_REC)
V49_pdf = ROOT.RooProdPdf('V49_pdf', 'V49 peak pdf', V49_ion_pdf, V49_rec_pdf)
V49_pdf_eff = ROOT.RooProdPdf('V49_pdf_eff', 'eff corr V49 peak pdf', V49_pdf, total_efficiency_pdf)
N_V49 = ROOT.RooRealVar('N_V49', 'evts of V49 peak (4.97keV)', 16., 0., events)

Cr51_ion_energy = ROOT.RooRealVar('Cr51_ion_energy', 'Cr51 peak ion energy', 5.46)
Cr51_ion_pos = ROOT.RooFormulaVar('Cr51_ion_pos', '@0*@1', ROOT.RooArgList(Cr51_ion_energy, ion_scaling))
Cr51_ion = ROOT.RooGaussian('Cr51_ion_pdf', 'Cr51 peak gauss pdf with shifted mean', ION, Cr51_ion_pos, SIGMA_ION)
Cr51_rec_energy = ROOT.RooRealVar('Cr51_rec_energy', 'v_rec_energy', functions.ER_CENTROID.GetX(Cr51_ion_energy.getVal()))
Cr51_rec_pos = ROOT.RooFormulaVar('Cr51_rec_pos', '@0*@1', ROOT.RooArgList(Cr51_rec_energy, rec_scaling))
Cr51_rec = ROOT.RooGaussian('Cr51_rec_pdf', 'Cr51_rec_pdf with shifted mean', REC, Cr51_rec_pos, SIGMA_REC)
Cr51_pdf = ROOT.RooProdPdf('Cr51_pdf', 'Cr51 peak pdf', Cr51_ion, Cr51_rec)
N_Cr51 = ROOT.RooRealVar('N_Cr51', 'evts of 51Cr peak (5.46keV)', 11., 0., events)

Mn54_ion_energy = ROOT.RooRealVar('Mn54_ion_energy', 'Mn54_ion_energy', 5.99)
Mn54_ion_pos = ROOT.RooFormulaVar('Mn54_ion_pos', '@0*@1', ROOT.RooArgList(Mn54_ion_energy, ion_scaling))
Mn54_ion = ROOT.RooGaussian('Mn54_ion_pdf', 'Mn54_ion_pdf with shifted mean', ION, Mn54_ion_pos, SIGMA_ION)
Mn54_rec_energy = ROOT.RooRealVar('Mn54_rec_energy', 'Mn54_rec_energy', functions.ER_CENTROID.GetX(Mn54_ion_energy.getVal()))
Mn54_rec_pos = ROOT.RooFormulaVar('Mn54_rec_pos', '@0*@1', ROOT.RooArgList(Mn54_rec_energy, rec_scaling))
Mn54_rec = ROOT.RooGaussian('Mn54_rec_pdf', 'Mn54_rec_pdf with shifted mean', REC, Mn54_rec_pos, SIGMA_REC)
Mn54_pdf = ROOT.RooProdPdf('Mn54_pdf', 'Mn54 peak pdf', Mn54_ion, Mn54_rec)
N_Mn54 = ROOT.RooRealVar('N_Mn54', 'evts of 54Mn peak (5.99keV)', 4., 0., events)

Fe55_ion_energy = ROOT.RooRealVar('Fe55_ion_energy', 'Fe55_ion_energy', 6.54)
Fe55_ion_pos = ROOT.RooFormulaVar('Fe55_ion_pos', '@0*@1', ROOT.RooArgList(Fe55_ion_energy, ion_scaling))
Fe55_ion = ROOT.RooGaussian('Fe55_ion_pdf', 'Fe55_ion_pdf with shifted mean', ION, Fe55_ion_pos, SIGMA_ION)
Fe55_rec_energy = ROOT.RooRealVar('Fe55_rec_energy', 'Fe55_rec_energy', functions.ER_CENTROID.GetX(Fe55_ion_energy.getVal()))
Fe55_rec_pos = ROOT.RooFormulaVar('Fe55_rec_pos', '@0*@1', ROOT.RooArgList(Fe55_rec_energy, rec_scaling))
Fe55_rec = ROOT.RooGaussian('Fe55_rec_pdf', 'Fe55_rec_pdf with shifted mean', REC, Fe55_rec_pos, SIGMA_REC)
Fe55_pdf = ROOT.RooProdPdf('Fe55_pdf', 'Fe55 peak pdf', Fe55_ion, Fe55_rec)
N_Fe55 = ROOT.RooRealVar('N_Fe55', 'evts of 55Fe peak (6.54keV)', 31., 0., events)

Co57_ion_energy = ROOT.RooRealVar('Co57_ion_energy', 'Co57_ion_energy', 7.11)
Co57_ion_pos = ROOT.RooFormulaVar('Co57_ion_pos', '@0*@1', ROOT.RooArgList(Co57_ion_energy, ion_scaling))
Co57_ion = ROOT.RooGaussian('Co57_ion_pdf', 'Co57_ion_pdf with shifted mean', ION, Co57_ion_pos, SIGMA_ION)
Co57_rec_energy = ROOT.RooRealVar('Co57_rec_energy', 'Co57_rec_energy', functions.ER_CENTROID.GetX(Co57_ion_energy.getVal()))
Co57_rec_pos = ROOT.RooFormulaVar('Co57_rec_pos', '@0*@1', ROOT.RooArgList(Co57_rec_energy, rec_scaling))
Co57_rec = ROOT.RooGaussian('Co57_rec_pdf', 'Co57_rec_pdf with shifted mean', REC, Co57_rec_pos, SIGMA_REC)
Co57_pdf = ROOT.RooProdPdf('Co57_pdf', 'Co57 peak pdf', Co57_ion, Co57_rec)
N_Co57 = ROOT.RooRealVar('N_Co57', 'evts of 57Co peak (7.11keV)', 2., 0., events)

Zn65_ion_energy = ROOT.RooRealVar('Zn65_ion_energy', 'Zn65_ion_energy', 8.98)
Zn65_ion_pos = ROOT.RooFormulaVar('Zn65_ion_pos', '@0*@1', ROOT.RooArgList(Zn65_ion_energy, ion_scaling))
Zn65_ion = ROOT.RooGaussian('Zn65_ion_pdf', 'Zn65_ion_pdf with shifted mean', ION, Zn65_ion_pos, SIGMA_ION)
Zn65_rec_energy = ROOT.RooRealVar('Zn65_rec_energy', 'Zn65_rec_energy', functions.ER_CENTROID.GetX(Zn65_ion_energy.getVal()))
Zn65_rec_pos = ROOT.RooFormulaVar('Zn65_rec_pos', '@0*@1', ROOT.RooArgList(Zn65_rec_energy, rec_scaling))
Zn65_rec = ROOT.RooGaussian('Zn65_rec_pdf', 'Zn65_rec_pdf with shifted mean', REC, Zn65_rec_pos, SIGMA_REC)
Zn65_pdf = ROOT.RooProdPdf('Zn65_pdf', 'Zn65 peak pdf', Zn65_ion, Zn65_rec)
N_Zn65 = ROOT.RooRealVar('N_Zn65', 'evts of 65Zn peak (8.98keV)', 110., 0., events)

Ga68_ion_energy = ROOT.RooRealVar('Ga68_ion_energy', 'Ga68_ion_energy', 9.66)
Ga68_ion_pos = ROOT.RooFormulaVar('Ga68_ion_pos', '@0*@1', ROOT.RooArgList(Ga68_ion_energy, ion_scaling))
Ga68_ion = ROOT.RooGaussian('Ga68_ion_pdf', 'Ga68_ion_pdf with shifted mean', ION, Ga68_ion_pos, SIGMA_ION)
Ga68_rec_energy = ROOT.RooRealVar('Ga68_rec_energy', 'Ga68_rec_energy', functions.ER_CENTROID.GetX(Ga68_ion_energy.getVal()))
Ga68_rec_pos = ROOT.RooFormulaVar('Ga68_rec_pos', '@0*@1', ROOT.RooArgList(Ga68_rec_energy, rec_scaling))
Ga68_rec = ROOT.RooGaussian('Ga68_rec_pdf', 'Ga68_rec_pdf with shifted mean', REC, Ga68_rec_pos, SIGMA_REC)
Ga68_pdf = ROOT.RooProdPdf('Ga68_pdf', 'Ga68 peak pdf', Ga68_ion, Ga68_rec)
N_Ga68 = ROOT.RooRealVar('N_Ga68', 'evts of 68Ga peak (9.66keV)', 32., 0., events)

Ge68_ion_energy = ROOT.RooRealVar('Ge68_ion_energy', 'Ge68_ion_energy', 10.37)
Ge68_ion_pos = ROOT.RooFormulaVar('Ge68_ion_pos', '@0*@1', ROOT.RooArgList(Ge68_ion_energy, ion_scaling))
Ge68_ion = ROOT.RooGaussian('Ge68_ion_pdf', 'Ge68_ion_pdf with shifted mean', ION, Ge68_ion_pos, SIGMA_ION)
Ge68_rec_energy = ROOT.RooRealVar('Ge68_rec_energy', 'Ge68_rec_energy', functions.ER_CENTROID.GetX(Ge68_ion_energy.getVal()))
Ge68_rec_pos = ROOT.RooFormulaVar('Ge68_rec_pos', '@0*@1', ROOT.RooArgList(Ge68_rec_energy, rec_scaling))
Ge68_rec = ROOT.RooGaussian('Ge68_rec_pdf', 'Ge68 peak pdf in rec', REC, Ge68_rec_pos, SIGMA_REC)
Ge68_pdf = ROOT.RooProdPdf('Ge68_pdf', 'Ge68 peak pdf', Ge68_ion, Ge68_rec)
N_Ge68 = ROOT.RooRealVar('N_Ge68', 'evts of 68Ge peak (10.37keV)', 0., 0., events)


# Definition of WIMP signal and pdf
signal_hist = functions.wimp_signal(WIMP_MASS, SIGMA_ION.getVal(), SIGMA_REC.getVal())
signal_hist.Multiply(total_efficiency)
signal_datahist = ROOT.RooDataHist('signal_datahist', 'signal_datahist', 
                                   ROOT.RooArgList(REC, ION), signal_hist)
signal_pdf = ROOT.RooHistPdf('signal_pdf', 'signal_pdf', 
                             ROOT.RooArgSet(REC, ION), signal_datahist)
N_signal = ROOT.RooRealVar('N_signal', 'WIMP signal events', 0., -10., 10.)


# Definition of flat gamma background pdf
flat_gamma_bckgd_hist = functions.flat_gamma_bckgd(SIGMA_ION.getVal(), SIGMA_REC.getVal())
flat_gamma_bckgd_hist.Multiply(total_efficiency)
flat_gamma_bckgd_datahist = ROOT.RooDataHist('flat_gamma_bckgd_datahist', 'flat_gamma_bckgd_datahist', 
                                             ROOT.RooArgList(REC, ION), flat_gamma_bckgd_hist)
flat_gamma_bckgd_pdf = ROOT.RooHistPdf('flat_gamma_bckgd_pdf', 'flat_gamma_bckgd_pdf', 
                                       ROOT.RooArgSet(REC, ION), flat_gamma_bckgd_datahist)
N_flat = ROOT.RooRealVar('N_flat', 'bckgd events', 70., 0., events)


# Definition of background only as well as background plus signal pdf
bckgd_and_sig_pdfs = ROOT.RooArgList(signal_pdf, flat_gamma_bckgd_pdf, V49_pdf, Cr51_pdf, Co57_pdf, Mn54_pdf, Fe55_pdf, Zn65_pdf, Ga68_pdf)  # Ge68 not in energy range and therefore excluded
bckgd_only_pdfs = ROOT.RooArgList(flat_gamma_bckgd_pdf, V49_pdf, Cr51_pdf, Co57_pdf, Mn54_pdf, Fe55_pdf, Zn65_pdf, Ga68_pdf)

bckgd_and_sig_params = ROOT.RooArgList(N_signal, N_flat, N_V49, N_Cr51, N_Co57, N_Mn54, N_Fe55, N_Zn65, N_Ga68)
bckgd_only_params = ROOT.RooArgList(N_flat, N_V49, N_Cr51, N_Co57, N_Mn54, N_Fe55, N_Zn65, N_Ga68)

bckgd_and_sig_pdf = ROOT.RooAddPdf('bckgd_and_sig_pdf', 'bckgd_and_sig_pdf', bckgd_and_sig_pdfs, bckgd_and_sig_params)
bckgd_only_pdf = ROOT.RooAddPdf('bckgd_only_pdf', 'bckgd_only_pdf', bckgd_only_pdfs, bckgd_only_params)

bckgd_and_sig_hist = bckgd_and_sig_pdf.createHistogram('bckgd_and_sig_hist', REC, rf.Binning(int(E_REC_MAX*10)), 
                                       rf.YVar(ION, rf.Binning(int(E_ION_MAX*10)))
                                       )
bckgd_and_sig_hist.SetTitle('ID3 WIMP search data + best fit PDF')


# Create negative log likelihood (NLL) object manually and minimize it
nll = ROOT.RooNLLVar('nll', 'nll', bckgd_and_sig_pdf, realdata, rf.Extended(ROOT.kTRUE), 
                     rf.PrintEvalErrors(0), rf.Verbose(ROOT.kFALSE))
minuit = ROOT.RooMinuit(nll)
minuit.migrad()  # find minimum
minuit.hesse()  # find symmetric errors
minuit.minos()  # find asymmetric errors
FitResult = minuit.save('realfit', 'fit to real data')
FitParams = FitResult.floatParsFinal()
NumFitParams = FitParams.getSize()
FitResult.Print('v')


# Calculate cross section limit from fit results and output of results
signal_parameter = FitResult.floatParsFinal().find('N_signal')
N_sig = signal_parameter.getVal()
sigma_n_sig_low = signal_parameter.getAsymErrorLo()
sigma_n_sig_high = signal_parameter.getAsymErrorHi()
N_max = N_sig + 1.28 * sigma_n_sig_high
livetime = 197.  # ID3 livetime after cuts in days
mass = 0.160  # ID3 detector mass in kg
rate = signal_hist.Integral('WIDTH')  # option width absolutely necessary
N_wimp = rate * livetime * mass * 1e6
xs_max = N_max / N_wimp
result_overview = {'N_sig' : N_sig, 'N_max' : N_max, 'rate' : rate, 'xs_max' : xs_max, 'N_wimp' : N_wimp}

print '{0:9} | {1:10} | {2:8} | {3:8} | {4:10} | {5:10}'.format('WIMP_MASS', 'Rate', 'N_signal', 
                                                                'ErrorLow', 'ErrorHigh', 'XS-limit [pb]')
print '-'*80
print '{0:9} | {1:10} | {2:8} | {3:8} | {4:10} | {5:10}'.format(WIMP_MASS, rate, N_sig, 
                                                                sigma_n_sig_low, sigma_n_sig_high, xs_max)


# Definition of RooFit projections over both energies starting with appropriate binning
recbins = int(E_REC_MAX*2)
ionbins = int(E_ION_MAX*5)

ionframe = ION.frame()
ionframe.SetTitle('Projection in E_{ion}')
realdata.plotOn(ionframe, rf.Name('data'), 
                rf.Binning(ionbins), rf.MarkerSize(1.0))
bckgd_and_sig_pdf.plotOn(ionframe, rf.Components("flat_gamma_bckgd_pdf"), 
                 rf.LineColor(ROOT.kGreen), rf.LineWidth(2))
bckgd_and_sig_pdf.plotOn(ionframe, rf.Name('model'), 
                 rf.LineColor(ROOT.kBlue), rf.LineWidth(2), rf.LineStyle(ROOT.kSolid))
bckgd_and_sig_pdf.plotOn(ionframe, rf.Components("V49_ion_pdf"), 
                 rf.LineColor(ROOT.kRed), rf.LineWidth(2), rf.LineStyle(ROOT.kDashed))
bckgd_and_sig_pdf.plotOn(ionframe, rf.Components("Cr51_ion_pdf"), 
                 rf.LineColor(ROOT.kRed), rf.LineWidth(2), rf.LineStyle(ROOT.kDashed))
bckgd_and_sig_pdf.plotOn(ionframe, rf.Components("Mn54_ion_pdf"), 
                 rf.LineColor(ROOT.kRed), rf.LineWidth(2), rf.LineStyle(ROOT.kDashed))
bckgd_and_sig_pdf.plotOn(ionframe, rf.Components("Fe55_ion_pdf"), 
                 rf.LineColor(ROOT.kRed), rf.LineWidth(2), rf.LineStyle(ROOT.kDashed))
bckgd_and_sig_pdf.plotOn(ionframe, rf.Components("Co57_ion_pdf"), 
                 rf.LineColor(ROOT.kRed), rf.LineWidth(2), rf.LineStyle(ROOT.kDashed))
bckgd_and_sig_pdf.plotOn(ionframe, rf.Components("Zn65_ion_pdf"), 
                 rf.LineColor(ROOT.kRed), rf.LineWidth(2), rf.LineStyle(ROOT.kDashed))
#bckgd_and_sig_pdf.plotOn(ionframe, rf.Components("Ge68_ion_pdf"), 
                 #rf.LineColor(ROOT.kRed), rf.LineWidth(2), rf.LineStyle(ROOT.kDashed))
bckgd_and_sig_pdf.plotOn(ionframe, rf.Components("Ga68_ion_pdf"), 
                 rf.LineColor(ROOT.kRed), rf.LineWidth(2), rf.LineStyle(ROOT.kDashed))
bckgd_and_sig_pdf.plotOn(ionframe, rf.Components("signal_pdf"), 
                 rf.LineColor(ROOT.kMagenta), rf.LineWidth(3), rf.LineStyle(ROOT.kSolid))
bckgd_and_sig_pdf.paramOn(ionframe, rf.Format('NEU', rf.AutoPrecision(2)), 
                  rf.Layout(0.1, 0.55, 0.9), rf.ShowConstants(ROOT.kFALSE))

recframe = REC.frame()
recframe.SetTitle('Projection in E_{rec}')
realdata.plotOn(recframe, rf.Name("data"), 
                rf.Binning(recbins), rf.MarkerColor(ROOT.kBlack), rf.MarkerSize(1.0))
bckgd_and_sig_pdf.plotOn(recframe, rf.Components("flat_gamma_bckgd_pdf"), 
                 rf.LineColor(ROOT.kGreen), rf.LineWidth(2), rf.LineStyle(ROOT.kSolid))
bckgd_and_sig_pdf.plotOn(recframe, rf.Name('model'), 
                 rf.LineColor(ROOT.kBlue), rf.LineWidth(2), rf.LineStyle(ROOT.kSolid))
bckgd_and_sig_pdf.plotOn(recframe, rf.Components("V49_rec_pdf"), 
                 rf.LineColor(ROOT.kRed), rf.LineWidth(2), rf.LineStyle(ROOT.kDashed))
bckgd_and_sig_pdf.plotOn(recframe, rf.Components("Cr51_rec_pdf"), 
                 rf.LineColor(ROOT.kRed), rf.LineWidth(2), rf.LineStyle(ROOT.kDashed))
bckgd_and_sig_pdf.plotOn(recframe, rf.Components("Mn54_rec_pdf"), 
                 rf.LineColor(ROOT.kRed), rf.LineWidth(2), rf.LineStyle(ROOT.kDashed))
bckgd_and_sig_pdf.plotOn(recframe, rf.Components("Fe55_rec_pdf"), 
                 rf.LineColor(ROOT.kRed), rf.LineWidth(2), rf.LineStyle(ROOT.kDashed))
bckgd_and_sig_pdf.plotOn(recframe, rf.Components("Co57_rec_pdf"), 
                 rf.LineColor(ROOT.kRed), rf.LineWidth(2), rf.LineStyle(ROOT.kDashed))
bckgd_and_sig_pdf.plotOn(recframe, rf.Components("Zn65_rec_pdf"), 
                 rf.LineColor(ROOT.kRed), rf.LineWidth(2), rf.LineStyle(ROOT.kDashed))
#bckgd_and_sig_pdf.plotOn(recframe, rf.Components("Ge68_rec_pdf"), 
                 #rf.LineColor(ROOT.kRed), rf.LineWidth(2), rf.LineStyle(ROOT.kDashed))
bckgd_and_sig_pdf.plotOn(recframe, rf.Components("Ga68_rec_pdf"), 
                 rf.LineColor(ROOT.kRed), rf.LineWidth(2), rf.LineStyle(ROOT.kDashed))
bckgd_and_sig_pdf.plotOn(recframe, rf.Components("signal_pdf"), 
                 rf.LineColor(ROOT.kMagenta), rf.LineWidth(3), rf.LineStyle(ROOT.kSolid))


# Plotting of fit results (PDFs and projected spectra)
c1 = ROOT.TCanvas('c1', 'Fit result overview for %i GeV'%WIMP_MASS, 1000, 750)
c1.Divide(2, 2)
c1.cd(1)
c1.cd(1).SetLogz()
bckgd_and_sig_hist.Draw('COL')
bckgd_and_sig_hist.SetStats(0)
realdata_scatter.SetMarkerColor(ROOT.kBlack)
realdata_scatter.SetMarkerStyle(ROOT.kFullDotLarge)
realdata_scatter.Draw('SAMESP')
functions.ER_CENTROID.SetLineColor(ROOT.kBlack)
functions.ER_CENTROID.SetLineWidth(1)
functions.ER_CENTROID.DrawCopy('SAME')
functions.NR_CENTROID.SetLineColor(ROOT.kBlack)
functions.NR_CENTROID.SetLineWidth(1)
functions.NR_CENTROID.DrawCopy('SAME')
c1.cd(3)
signal_hist.Draw('CONT0')
signal_hist.SetStats(0)
signal_hist.SetContour(30)
realdata_scatter.Draw('SAMESP')
functions.ER_CENTROID.DrawCopy('SAME')
functions.NR_CENTROID.DrawCopy('SAME')
c1.cd(2)
ionframe.Draw()
c1.cd(4)
recframe.Draw()


# Creation of Monte Carlo toy event sets and output
if NUM_MC_SETS:
    MC_study = ROOT.RooMCStudy(bckgd_only_pdf, ROOT.RooArgSet(REC, ION), rf.FitModel(bckgd_and_sig_pdf), rf.Silence(), 
                               rf.Extended(ROOT.kTRUE), rf.FitOptions(rf.Save(ROOT.kTRUE)))  # , rf.FitModel(bckgd_and_sig_pdf)
    MC_study.generateAndFit(NUM_MC_SETS)

    ParamNLLFrameList = []
    ParamDistriFrameList = []
    ParamPullFrameList = []
    ParamLineList = []

    c2 = ROOT.TCanvas('c2', 'Fit + MC toy set statistics for %i GeV'%WIMP_MASS)
    c2.Divide(1, NumFitParams, 0.001, 0.001)

    nll_line = ROOT.TLine()
    nll_line.SetLineWidth(2)
    nll_line.SetLineStyle(7)
    nll_line.SetLineColor(ROOT.kBlue)

    zeroline = ROOT.TLine()
    zeroline.SetLineWidth(2)
    zeroline.SetLineStyle(7)
    zeroline.SetLineColor(ROOT.kBlack)

    for i in range(NumFitParams):
        parameter = FitParams[i]
        paramname = parameter.GetName()
        paramvalue = parameter.getVal()

        paramline = ROOT.TLine()
        paramline.SetLineWidth(2)
        paramline.SetLineStyle(1)
        paramline.SetLineColor(ROOT.kBlue)
        ParamLineList.append(paramline)

        pad = c2.cd(i+1)
        pad.Divide(3, 1, 0.001, 0.001)

        ParamNLLFrame = parameter.frame()
        ParamNLLFrame.SetTitle('NLL fit '+paramname)
        nll.plotOn(ParamNLLFrame, rf.Precision(1e-5), rf.ShiftToZero())
        ParamNLLFrame.SetMaximum(20.)
        ParamNLLFrame.SetMinimum(0.)
        ParamNLLFrameList.append(ParamNLLFrame)
        pad.cd(1)
        ParamNLLFrame.Draw()
        ROOT.gPad.Update()
        paramline.DrawLine(paramvalue, ROOT.gPad.GetUymin(), paramvalue, ROOT.gPad.GetUymax())

        ParamDistriFrame = MC_study.plotParam(parameter)
        ParamDistriFrame.SetTitle('MC distri '+paramname)
        ParamDistriFrameList.append(ParamDistriFrame)
        pad.cd(2)
        ParamDistriFrame.Draw()
        ParamDistriFrame.getHist().Fit('gaus', 'QEM')
        ROOT.gPad.Update()
        paramline.DrawLine(paramvalue, ROOT.gPad.GetUymin(), paramvalue, ROOT.gPad.GetUymax())

        #ParamPullFrame = MC_study.plotPull(parameter, rf.FitGauss())
        #ParamPullFrame.SetTitle('MC pull distri '+paramname)
        #ParamPullFrameList.append(ParamPullFrame)
        #pad.cd(3)
        #ParamPullFrame.Draw()
        #ROOT.gPad.Update()
        #zeroline.DrawLine(0, ROOT.gPad.GetUymin(), 0, ROOT.gPad.GetUymax())

    c2.SetCanvasSize(1200, 2400)


# Extra canvas for signal parameters only (including MC toy event set fits)
if NUM_MC_SETS and WIMP_MASS:
    parameter = FitResult.floatParsFinal().find('N_signal')
    paramname = parameter.GetName()
    paramvalue = parameter.getVal()

    nllvalue = nll.getVal()

    c3 = ROOT.TCanvas('c3', 'N_signal fit & MC toy set statistics for %i GeV'%WIMP_MASS, 1000, 750)
    c3.Divide(2, 2)
    c3.cd(1)
    ParamNLLFrame = parameter.frame()
    ParamNLLFrame.SetTitle('NLL fit '+paramname)
    nll.plotOn(ParamNLLFrame, rf.Precision(1e-5), rf.ShiftToZero())
    #ParamNLLFrame.SetMaximum(20.)
    #ParamNLLFrame.SetMinimum(0.)
    ParamNLLFrame.Draw()
    ROOT.gPad.Update()
    paramline.DrawLine(paramvalue, ROOT.gPad.GetUymin(), paramvalue, ROOT.gPad.GetUymax())
    zeroline.DrawLine(0, ROOT.gPad.GetUymin(), 0, ROOT.gPad.GetUymax())

    c3.cd(2)
    ParamDistriFrame = MC_study.plotParam(parameter)  #, rf.FrameBins(200), rf.FrameRange(0., 10.))
    ParamDistriFrame.SetTitle('MC distri '+paramname)
    ParamDistriFrame.Draw()
    #ParamDistriFrame.getHist().Fit('gaus', 'QEM')
    ROOT.gPad.Update()
    paramline.DrawLine(paramvalue, ROOT.gPad.GetUymin(), paramvalue, ROOT.gPad.GetUymax())
    zeroline.DrawLine(0, ROOT.gPad.GetUymin(), 0, ROOT.gPad.GetUymax())

    #c3.cd(3)
    #ParamPullFrame = MC_study.plotPull(parameter, rf.FrameBins(100), rf.FrameRange(-5, 5), rf.FitGauss(ROOT.kTRUE))
    #ParamPullFrame.SetTitle('MC pull distri '+paramname)
    #ParamPullFrame.Draw()
    #ROOT.gPad.Update()
    #zeroline.DrawLine(0, ROOT.gPad.GetUymin(), 0, ROOT.gPad.GetUymax())

    c3.cd(4)
    MCnllframe = MC_study.plotNLL()
    MCnllframe.SetTitle('NLL distri')
    MCnllframe.Draw()
    MCnllframe.getHist().Fit('gaus', 'QEM')
    ROOT.gPad.Update()
    nll_line.SetLineStyle(1)
    nll_line.DrawLine(nllvalue, ROOT.gPad.GetUymin(), nllvalue, ROOT.gPad.GetUymax())
    paramline.DrawLine(paramvalue, ROOT.gPad.GetUymin(), paramvalue, ROOT.gPad.GetUymax())


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
