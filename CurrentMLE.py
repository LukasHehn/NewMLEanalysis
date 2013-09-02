#!/usr/bin/env python

######################################################
#
# Calculation of 90% C.L. cross section limits using Maximum Likelihood methods based on RooFit
# Lukas Hehn, 2013
#
######################################################

import ROOT
import Functions
from ROOT import RooFit as rf



# Global variables used in the script
DataFile = '/kalinka/home/hehn/PhD/LowMassEric/ID3_eventlist.txt'
EnergyIonMax = 10.
EnergyRecMax = 20.
wimp_mass = 8  # if False, signal is switched off
MC_sets = int(1e4)  # number of MC simulations: 0 means none at all
SavePlots = False


# Definition of maximum likelihood observables
ion = ROOT.RooRealVar('ion', 'E_{ion}', 0., EnergyIonMax, 'keV_{ee}')
rec = ROOT.RooRealVar('rec', 'E_{rec}', 0., EnergyRecMax, 'keV_{nr}')


# Detector specific parameters like resolutions
Voltage = 6.4
E_thresh = 3.874
FWHM_heat = 0.82  # value for ID3 adopted from Eric
FWHM_rec = Functions.RecoilResolutionFromHeat(FWHM_heat, Voltage, 10.)
FWHM_ion = 0.72  # value for ID3 adopted from Eric
sigma_rec = ROOT.RooRealVar('sigma_rec', 'recoil energy resolution', FWHM_rec/2.35)
sigma_ion = ROOT.RooRealVar('sigma_ion', 'ionization energy resolution', FWHM_ion/2.35)


# Calculation of specific detector efficiency and pdf
total_efficiency = Functions.Simple2DEfficiencyID3(E_thresh, FWHM_rec/2.35)
total_efficiency_datahist = ROOT.RooDataHist('total_efficiency_datahist', 'total_efficiency_datahist', 
                                             ROOT.RooArgList(rec, ion), total_efficiency)
total_efficiency_pdf = ROOT.RooHistPdf('total_efficiency_pdf', 'total_efficiency_pdf', 
                                       ROOT.RooArgSet(rec, ion), total_efficiency_datahist)


# Read in of data set and output as scatter/graph
realdata = ROOT.RooDataSet.read(DataFile, ROOT.RooArgList(rec, ion))
realdata_scatter = realdata.createHistogram(rec, ion, int(EnergyRecMax*10), int(EnergyIonMax*10))
realdata_graph = Functions.TGraphFromDataSet(realdata)
events = int(realdata.numEntries())


# Definition of gamma peaks
ion_scaling = ROOT.RooRealVar('ion_scaling', 'scaling factor ionization energy', 1.0, 0.95, 1.05)
rec_scaling = ROOT.RooRealVar('rec_scaling', 'scaling factor recoil energy', 1.0, 0.95, 1.05)
ion_scaling.setConstant(ROOT.kTRUE)
rec_scaling.setConstant(ROOT.kTRUE)
Functions.ER_centroid.SetParameter(0, 6.4)  # ER_centroid for calculation of peak position in Erec

V49_ion_energy = ROOT.RooRealVar('V49_ion_energy', 'V49 peak ion energy', 4.97)
V49_ion_pos = ROOT.RooFormulaVar('V49_ion_pos', '@0*@1', ROOT.RooArgList(V49_ion_energy, ion_scaling))
V49_ion_pdf = ROOT.RooGaussian('V49_ion_pdf', 'V49 peak gauss pdf in ion', ion, V49_ion_pos, sigma_ion)
V49_rec_energy = ROOT.RooRealVar('V49_rec_energy', 'recoil energy V49 peak', Functions.ER_centroid.GetX(V49_ion_energy.getVal()))
V49_rec_pos = ROOT.RooFormulaVar('V49_rec_pos', '@0*@1', ROOT.RooArgList(V49_rec_energy, rec_scaling))
V49_rec_pdf = ROOT.RooGaussian('V49_rec_pdf', 'V49 peak pdf in recoil energy', rec, V49_rec_pos, sigma_rec)
V49_pdf = ROOT.RooProdPdf('V49_pdf', 'V49 peak pdf', V49_ion_pdf, V49_rec_pdf)
V49_pdf_eff = ROOT.RooProdPdf('V49_pdf_eff', 'eff corr V49 peak pdf', V49_pdf, total_efficiency_pdf)
N_V49 = ROOT.RooRealVar('N_V49', 'evts of V49 peak (4.97keV)', 16., 0., events)

Cr51_ion_energy = ROOT.RooRealVar('Cr51_ion_energy', 'Cr51 peak ion energy', 5.46)
Cr51_ion_pos = ROOT.RooFormulaVar('Cr51_ion_pos', '@0*@1', ROOT.RooArgList(Cr51_ion_energy, ion_scaling))
Cr51_ion = ROOT.RooGaussian('Cr51_ion_pdf', 'Cr51 peak gauss pdf with shifted mean', ion, Cr51_ion_pos, sigma_ion)
Cr51_rec_energy = ROOT.RooRealVar('Cr51_rec_energy', 'v_rec_energy', Functions.ER_centroid.GetX(Cr51_ion_energy.getVal()))
Cr51_rec_pos = ROOT.RooFormulaVar('Cr51_rec_pos', '@0*@1', ROOT.RooArgList(Cr51_rec_energy, rec_scaling))
Cr51_rec = ROOT.RooGaussian('Cr51_rec_pdf', 'Cr51_rec_pdf with shifted mean', rec, Cr51_rec_pos, sigma_rec)
Cr51_pdf = ROOT.RooProdPdf('Cr51_pdf', 'Cr51 peak pdf', Cr51_ion, Cr51_rec)
N_Cr51 = ROOT.RooRealVar('N_Cr51', 'evts of 51Cr peak (5.46keV)', 11., 0., events)

Mn54_ion_energy = ROOT.RooRealVar('Mn54_ion_energy', 'Mn54_ion_energy', 5.99)
Mn54_ion_pos = ROOT.RooFormulaVar('Mn54_ion_pos', '@0*@1', ROOT.RooArgList(Mn54_ion_energy, ion_scaling))
Mn54_ion = ROOT.RooGaussian('Mn54_ion_pdf', 'Mn54_ion_pdf with shifted mean', ion, Mn54_ion_pos, sigma_ion)
Mn54_rec_energy = ROOT.RooRealVar('Mn54_rec_energy', 'Mn54_rec_energy', Functions.ER_centroid.GetX(Mn54_ion_energy.getVal()))
Mn54_rec_pos = ROOT.RooFormulaVar('Mn54_rec_pos', '@0*@1', ROOT.RooArgList(Mn54_rec_energy, rec_scaling))
Mn54_rec = ROOT.RooGaussian('Mn54_rec_pdf', 'Mn54_rec_pdf with shifted mean', rec, Mn54_rec_pos, sigma_rec)
Mn54_pdf = ROOT.RooProdPdf('Mn54_pdf', 'Mn54 peak pdf', Mn54_ion, Mn54_rec)
N_Mn54 = ROOT.RooRealVar('N_Mn54', 'evts of 54Mn peak (5.99keV)', 4., 0., events)

Fe55_ion_energy = ROOT.RooRealVar('Fe55_ion_energy', 'Fe55_ion_energy', 6.54)
Fe55_ion_pos = ROOT.RooFormulaVar('Fe55_ion_pos', '@0*@1', ROOT.RooArgList(Fe55_ion_energy, ion_scaling))
Fe55_ion = ROOT.RooGaussian('Fe55_ion_pdf', 'Fe55_ion_pdf with shifted mean', ion, Fe55_ion_pos, sigma_ion)
Fe55_rec_energy = ROOT.RooRealVar('Fe55_rec_energy', 'Fe55_rec_energy', Functions.ER_centroid.GetX(Fe55_ion_energy.getVal()))
Fe55_rec_pos = ROOT.RooFormulaVar('Fe55_rec_pos', '@0*@1', ROOT.RooArgList(Fe55_rec_energy, rec_scaling))
Fe55_rec = ROOT.RooGaussian('Fe55_rec_pdf', 'Fe55_rec_pdf with shifted mean', rec, Fe55_rec_pos, sigma_rec)
Fe55_pdf = ROOT.RooProdPdf('Fe55_pdf', 'Fe55 peak pdf', Fe55_ion, Fe55_rec)
N_Fe55 = ROOT.RooRealVar('N_Fe55', 'evts of 55Fe peak (6.54keV)', 31., 0., events)

Co57_ion_energy = ROOT.RooRealVar('Co57_ion_energy', 'Co57_ion_energy', 7.11)
Co57_ion_pos = ROOT.RooFormulaVar('Co57_ion_pos', '@0*@1', ROOT.RooArgList(Co57_ion_energy, ion_scaling))
Co57_ion = ROOT.RooGaussian('Co57_ion_pdf', 'Co57_ion_pdf with shifted mean', ion, Co57_ion_pos, sigma_ion)
Co57_rec_energy = ROOT.RooRealVar('Co57_rec_energy', 'Co57_rec_energy', Functions.ER_centroid.GetX(Co57_ion_energy.getVal()))
Co57_rec_pos = ROOT.RooFormulaVar('Co57_rec_pos', '@0*@1', ROOT.RooArgList(Co57_rec_energy, rec_scaling))
Co57_rec = ROOT.RooGaussian('Co57_rec_pdf', 'Co57_rec_pdf with shifted mean', rec, Co57_rec_pos, sigma_rec)
Co57_pdf = ROOT.RooProdPdf('Co57_pdf', 'Co57 peak pdf', Co57_ion, Co57_rec)
N_Co57 = ROOT.RooRealVar('N_Co57', 'evts of 57Co peak (7.11keV)', 2., 0., events)

Zn65_ion_energy = ROOT.RooRealVar('Zn65_ion_energy', 'Zn65_ion_energy', 8.98)
Zn65_ion_pos = ROOT.RooFormulaVar('Zn65_ion_pos', '@0*@1', ROOT.RooArgList(Zn65_ion_energy, ion_scaling))
Zn65_ion = ROOT.RooGaussian('Zn65_ion_pdf', 'Zn65_ion_pdf with shifted mean', ion, Zn65_ion_pos, sigma_ion)
Zn65_rec_energy = ROOT.RooRealVar('Zn65_rec_energy', 'Zn65_rec_energy', Functions.ER_centroid.GetX(Zn65_ion_energy.getVal()))
Zn65_rec_pos = ROOT.RooFormulaVar('Zn65_rec_pos', '@0*@1', ROOT.RooArgList(Zn65_rec_energy, rec_scaling))
Zn65_rec = ROOT.RooGaussian('Zn65_rec_pdf', 'Zn65_rec_pdf with shifted mean', rec, Zn65_rec_pos, sigma_rec)
Zn65_pdf = ROOT.RooProdPdf('Zn65_pdf', 'Zn65 peak pdf', Zn65_ion, Zn65_rec)
N_Zn65 = ROOT.RooRealVar('N_Zn65', 'evts of 65Zn peak (8.98keV)', 110., 0., events)

Ga68_ion_energy = ROOT.RooRealVar('Ga68_ion_energy', 'Ga68_ion_energy', 9.66)
Ga68_ion_pos = ROOT.RooFormulaVar('Ga68_ion_pos', '@0*@1', ROOT.RooArgList(Ga68_ion_energy, ion_scaling))
Ga68_ion = ROOT.RooGaussian('Ga68_ion_pdf', 'Ga68_ion_pdf with shifted mean', ion, Ga68_ion_pos, sigma_ion)
Ga68_rec_energy = ROOT.RooRealVar('Ga68_rec_energy', 'Ga68_rec_energy', Functions.ER_centroid.GetX(Ga68_ion_energy.getVal()))
Ga68_rec_pos = ROOT.RooFormulaVar('Ga68_rec_pos', '@0*@1', ROOT.RooArgList(Ga68_rec_energy, rec_scaling))
Ga68_rec = ROOT.RooGaussian('Ga68_rec_pdf', 'Ga68_rec_pdf with shifted mean', rec, Ga68_rec_pos, sigma_rec)
Ga68_pdf = ROOT.RooProdPdf('Ga68_pdf', 'Ga68 peak pdf', Ga68_ion, Ga68_rec)
N_Ga68 = ROOT.RooRealVar('N_Ga68', 'evts of 68Ga peak (9.66keV)', 32., 0., events)

Ge68_ion_energy = ROOT.RooRealVar('Ge68_ion_energy', 'Ge68_ion_energy', 10.37)
Ge68_ion_pos = ROOT.RooFormulaVar('Ge68_ion_pos', '@0*@1', ROOT.RooArgList(Ge68_ion_energy, ion_scaling))
Ge68_ion = ROOT.RooGaussian('Ge68_ion_pdf', 'Ge68_ion_pdf with shifted mean', ion, Ge68_ion_pos, sigma_ion)
Ge68_rec_energy = ROOT.RooRealVar('Ge68_rec_energy', 'Ge68_rec_energy', Functions.ER_centroid.GetX(Ge68_ion_energy.getVal()))
Ge68_rec_pos = ROOT.RooFormulaVar('Ge68_rec_pos', '@0*@1', ROOT.RooArgList(Ge68_rec_energy, rec_scaling))
Ge68_rec = ROOT.RooGaussian('Ge68_rec_pdf', 'Ge68 peak pdf in rec', rec, Ge68_rec_pos, sigma_rec)
Ge68_pdf = ROOT.RooProdPdf('Ge68_pdf', 'Ge68 peak pdf', Ge68_ion, Ge68_rec)
N_Ge68 = ROOT.RooRealVar('N_Ge68', 'evts of 68Ge peak (10.37keV)', 0., 0., events)


# Definition of WIMP signal and pdf
if wimp_mass:
    signal_hist = Functions.WimpSignal2DEric(wimp_mass, sigma_ion.getVal(), sigma_rec.getVal())
    signal_hist.Multiply(total_efficiency)
    signal_datahist = ROOT.RooDataHist('signal_datahist', 'signal_datahist', 
                                       ROOT.RooArgList(rec, ion), signal_hist)
    signal_pdf = ROOT.RooHistPdf('signal_pdf', 'signal_pdf', 
                                 ROOT.RooArgSet(rec, ion), signal_datahist)
    N_signal = ROOT.RooRealVar('N_signal', 'WIMP signal events', 0., 0., 10.)


# Definition of flat gamma background pdf
flat_gamma_bckgd_hist = Functions.FlatGammaBckgd2DEric(sigma_ion.getVal(), sigma_rec.getVal())
flat_gamma_bckgd_hist.Multiply(total_efficiency)
flat_gamma_bckgd_datahist = ROOT.RooDataHist('flat_gamma_bckgd_datahist', 'flat_gamma_bckgd_datahist', 
                                             ROOT.RooArgList(rec, ion), flat_gamma_bckgd_hist)
flat_gamma_bckgd_pdf = ROOT.RooHistPdf('flat_gamma_bckgd_pdf', 'flat_gamma_bckgd_pdf', 
                                       ROOT.RooArgSet(rec, ion), flat_gamma_bckgd_datahist)
N_flat = ROOT.RooRealVar('N_flat', 'bckgd events', 70., 0., events)


# Definition of final pdf depending on wheather WIMP mass is given or not
if wimp_mass:
    pdf_list = ROOT.RooArgList(signal_pdf, flat_gamma_bckgd_pdf, V49_pdf, Cr51_pdf, Co57_pdf, Mn54_pdf, Fe55_pdf, Zn65_pdf, Ga68_pdf)  # Ge68 not in energy range and therefore excluded
    parameter_list = ROOT.RooArgList(N_signal, N_flat, N_V49, N_Cr51, N_Co57, N_Mn54, N_Fe55, N_Zn65, N_Ga68)
    final_pdf = ROOT.RooAddPdf('final_pdf', 'final_pdf', pdf_list, parameter_list)
else:
    pdf_list = ROOT.RooArgList(flat_gamma_bckgd_pdf, V49_pdf, Cr51_pdf, Mn54_pdf, Fe55_pdf, Co57_pdf, Zn65_pdf, Ga68_pdf, Ge68_pdf)
    parameter_list = ROOT.RooArgList(N_flat, N_V49, N_Cr51, N_Mn54, N_Fe55, N_Co57, N_Zn65, N_Ga68, N_Ge68)
    final_pdf = ROOT.RooAddPdf('final_pdf', 'final_pdf', pdf_list, parameter_list)
    
final_hist = final_pdf.createHistogram('final_hist', rec, 
                                       rf.Binning(int(EnergyRecMax*10)), 
                                       rf.YVar(ion, rf.Binning(int(EnergyIonMax*10)))
                                       )


# Create negative log likelihood (NLL) object manually and minimize it
nll = ROOT.RooNLLVar('nll', 'nll', final_pdf, realdata, rf.Extended(ROOT.kTRUE), 
                     rf.PrintEvalErrors(2), rf.Verbose(ROOT.kFALSE))
minuit = ROOT.RooMinuit(nll)
minuit.migrad()  # find minimum
minuit.hesse()  # find symmetric errors
minuit.minos()  # find asymmetric errors
FitResult = minuit.save('realfit', 'fit to real data')
FitParams = FitResult.floatParsFinal()
NumFitParams = FitParams.getSize()
FitResult.Print('v')


# Calculate cross section limit from fit results and output of results
if wimp_mass:
    signal_parameter = FitResult.floatParsFinal().find('N_signal')
    wimp_events = signal_parameter.getVal()
    error_low = signal_parameter.getAsymErrorLo()
    error_high = signal_parameter.getAsymErrorHi()
    wimp_events_limit = wimp_events + 1.28 * error_high
    livetime = 197.  # ID3 livetime after cuts in days
    mass = 0.160  # ID3 detector mass in kg
    rate = signal_hist.Integral('WIDTH')  # option width absolutely necessary
    signal_events = rate * livetime * mass * 1e6
    sigma_limit = wimp_events_limit / signal_events
    result_overview = {'wimp events' : wimp_events,
                       'wimp events limit' : wimp_events_limit,
                       'rate' : rate,
                       'sigma limit' : sigma_limit,
                       }
    
    print '{0:9} | {1:10} | {2:8} | {3:8} | {4:10} | {5:10}'.format('wimp_mass', 'Rate', 'N_signal', 'ErrorLow', 'ErrorHigh', 'XS-limit [pb]')
    print '-'*80
    print '{0:9} | {1:10} | {2:8} | {3:8} | {4:10} | {5:10}'.format(wimp_mass, rate, wimp_events, error_low, error_high, sigma_limit)


# Definition of RooFit projections over both energies starting with appropriate binning
recbins = int(EnergyRecMax*2)
ionbins = int(EnergyIonMax*5)

ionframe = ion.frame()
ionframe.SetTitle('Projection in E_{ion}')
realdata.plotOn(ionframe, rf.Name('data'), 
                rf.Binning(ionbins), rf.MarkerSize(1.0))
final_pdf.plotOn(ionframe, rf.Components("flat_gamma_bckgd_pdf"), 
                 rf.LineColor(ROOT.kGreen), rf.LineWidth(2))
final_pdf.plotOn(ionframe, rf.Name('model'), 
                 rf.LineColor(ROOT.kBlue), rf.LineWidth(2), rf.LineStyle(ROOT.kSolid))
final_pdf.plotOn(ionframe, rf.Components("V49_ion_pdf"), 
                 rf.LineColor(ROOT.kRed), rf.LineWidth(2), rf.LineStyle(ROOT.kDashed))
final_pdf.plotOn(ionframe, rf.Components("Cr51_ion_pdf"), 
                 rf.LineColor(ROOT.kRed), rf.LineWidth(2), rf.LineStyle(ROOT.kDashed))
final_pdf.plotOn(ionframe, rf.Components("Mn54_ion_pdf"), 
                 rf.LineColor(ROOT.kRed), rf.LineWidth(2), rf.LineStyle(ROOT.kDashed))
final_pdf.plotOn(ionframe, rf.Components("Fe55_ion_pdf"), 
                 rf.LineColor(ROOT.kRed), rf.LineWidth(2), rf.LineStyle(ROOT.kDashed))
#final_pdf.plotOn(ionframe, rf.Components("Co57_ion_pdf"), 
                 #rf.LineColor(ROOT.kRed), rf.LineWidth(2), rf.LineStyle(ROOT.kDashed))
final_pdf.plotOn(ionframe, rf.Components("Zn65_ion_pdf"), 
                 rf.LineColor(ROOT.kRed), rf.LineWidth(2), rf.LineStyle(ROOT.kDashed))
final_pdf.plotOn(ionframe, rf.Components("Ge68_ion_pdf"), 
                 rf.LineColor(ROOT.kRed), rf.LineWidth(2), rf.LineStyle(ROOT.kDashed))
final_pdf.plotOn(ionframe, rf.Components("Ga68_ion_pdf"), 
                 rf.LineColor(ROOT.kRed), rf.LineWidth(2), rf.LineStyle(ROOT.kDashed))
final_pdf.plotOn(ionframe, rf.Components("signal_pdf"), 
                 rf.LineColor(ROOT.kMagenta), rf.LineWidth(3), rf.LineStyle(ROOT.kSolid))
final_pdf.paramOn(ionframe, rf.Format('NEU', rf.AutoPrecision(2)), 
                  rf.Layout(0.1, 0.55, 0.9), rf.ShowConstants(ROOT.kFALSE))

recframe = rec.frame()
recframe.SetTitle('Projection in E_{rec}')
realdata.plotOn(recframe, rf.Name("data"), 
                rf.Binning(recbins), rf.MarkerColor(ROOT.kBlack), rf.MarkerSize(1.0))
final_pdf.plotOn(recframe, rf.Components("flat_gamma_bckgd_pdf"), 
                 rf.LineColor(ROOT.kGreen), rf.LineWidth(2), rf.LineStyle(ROOT.kSolid))
final_pdf.plotOn(recframe, rf.Name('model'), 
                 rf.LineColor(ROOT.kBlue), rf.LineWidth(2), rf.LineStyle(ROOT.kSolid))
final_pdf.plotOn(recframe, rf.Components("V49_rec_pdf"), 
                 rf.LineColor(ROOT.kRed), rf.LineWidth(2), rf.LineStyle(ROOT.kDashed))
final_pdf.plotOn(recframe, rf.Components("Cr51_rec_pdf"), 
                 rf.LineColor(ROOT.kRed), rf.LineWidth(2), rf.LineStyle(ROOT.kDashed))
final_pdf.plotOn(recframe, rf.Components("Mn54_rec_pdf"), 
                 rf.LineColor(ROOT.kRed), rf.LineWidth(2), rf.LineStyle(ROOT.kDashed))
final_pdf.plotOn(recframe, rf.Components("Fe55_rec_pdf"), 
                 rf.LineColor(ROOT.kRed), rf.LineWidth(2), rf.LineStyle(ROOT.kDashed))
#final_pdf.plotOn(recframe, rf.Components("Co57_rec_pdf"), 
                 #rf.LineColor(ROOT.kRed), rf.LineWidth(2), rf.LineStyle(ROOT.kDashed))
final_pdf.plotOn(recframe, rf.Components("Zn65_rec_pdf"), 
                 rf.LineColor(ROOT.kRed), rf.LineWidth(2), rf.LineStyle(ROOT.kDashed))
final_pdf.plotOn(recframe, rf.Components("Ge68_rec_pdf"), 
                 rf.LineColor(ROOT.kRed), rf.LineWidth(2), rf.LineStyle(ROOT.kDashed))
final_pdf.plotOn(recframe, rf.Components("Ga68_rec_pdf"), 
                 rf.LineColor(ROOT.kRed), rf.LineWidth(2), rf.LineStyle(ROOT.kDashed))
final_pdf.plotOn(recframe, rf.Components("signal_pdf"), 
                 rf.LineColor(ROOT.kMagenta), rf.LineWidth(3), rf.LineStyle(ROOT.kSolid))


# Plotting of fit results (PDFs and projected spectra)
c1 = ROOT.TCanvas('c1', 'Fit Results', 1000, 750)
c1.Divide(2, 2)
# final pdf and data
c1.cd(1)
c1.cd(1).SetLogz()
final_hist.Draw('COL')
final_hist.SetStats(0)
final_hist.SetTitle('ID3 WIMP search data + best fit PDF')
realdata_graph.SetMarkerColor(ROOT.kBlack)
realdata_graph.SetMarkerStyle(ROOT.kFullDotLarge)
realdata_graph.Draw('SAMESP')
Functions.ER_centroid.SetLineColor(ROOT.kBlack)
Functions.ER_centroid.SetLineWidth(1)
Functions.ER_centroid.DrawCopy('SAME')
Functions.NR_centroid.SetLineColor(ROOT.kBlack)
Functions.NR_centroid.SetLineWidth(1)
Functions.NR_centroid.DrawCopy('SAME')
if wimp_mass:
    c1.cd(3)
    signal_hist.Draw('CONT0')
    signal_hist.SetStats(0)
    signal_hist.SetContour(30)
    realdata_graph.Draw('SAMESP')
    Functions.ER_centroid.DrawCopy('SAME')
    Functions.NR_centroid.DrawCopy('SAME')
c1.cd(2)
ionframe.Draw()
c1.cd(4)
recframe.Draw()


# Creation of Monte Carlo toy event sets and output
if MC_sets:
    MC_study = ROOT.RooMCStudy(final_pdf, ROOT.RooArgSet(rec, ion), rf.Silence(), 
                               rf.Extended(ROOT.kTRUE), rf.FitOptions(rf.Save(ROOT.kTRUE)))
    MC_study.generateAndFit(MC_sets)

    ParamNLLFrameList = []
    ParamDistriFrameList = []
    ParamPullFrameList = []
    ParamLineList = []

    c2 = ROOT.TCanvas('c2', 'Statistical interpretation (MC study)')
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

        paramnllframe = parameter.frame()
        paramnllframe.SetTitle('NLL fit '+paramname)
        nll.plotOn(paramnllframe, rf.Precision(1e-5), rf.ShiftToZero())
        paramnllframe.SetMaximum(20.)
        paramnllframe.SetMinimum(0.)
        ParamNLLFrameList.append(paramnllframe)
        pad.cd(1)
        paramnllframe.Draw()
        ROOT.gPad.Update()
        paramline.DrawLine(paramvalue, ROOT.gPad.GetUymin(), paramvalue, ROOT.gPad.GetUymax())

        paramdistriframe = MC_study.plotParam(parameter)
        paramdistriframe.SetTitle('MC distri '+paramname)
        ParamDistriFrameList.append(paramdistriframe)
        pad.cd(2)
        paramdistriframe.Draw()
        paramdistriframe.getHist().Fit('gaus', 'QEM')
        ROOT.gPad.Update()
        paramline.DrawLine(paramvalue, ROOT.gPad.GetUymin(), paramvalue, ROOT.gPad.GetUymax())

        parampullframe = MC_study.plotPull(parameter, rf.FitGauss())
        parampullframe.SetTitle('MC pull distri '+paramname)
        ParamPullFrameList.append(parampullframe)
        pad.cd(3)
        parampullframe.Draw()
        ROOT.gPad.Update()
        zeroline.DrawLine(0, ROOT.gPad.GetUymin(), 0, ROOT.gPad.GetUymax())

    c2.SetCanvasSize(1200, 2400)

# Extra canvas for signal parameters only (including MC toy event set fits)
if MC_sets and wimp_mass:
    parameter = FitResult.floatParsFinal().find('N_signal')
    paramname = parameter.GetName()
    paramvalue = parameter.getVal()

    nllvalue = nll.getVal()

    c3 = ROOT.TCanvas('c3', 'Fit Statistics WIMP Signal', 1000, 750)
    c3.Divide(2, 2)
    c3.cd(1)
    paramnllframe = parameter.frame()
    paramnllframe.SetTitle('NLL fit '+paramname)
    nll.plotOn(paramnllframe, rf.Precision(1e-5), rf.ShiftToZero())
    paramnllframe.SetMaximum(20.)
    paramnllframe.SetMinimum(0.)
    paramnllframe.Draw()
    ROOT.gPad.Update()
    paramline.DrawLine(paramvalue, ROOT.gPad.GetUymin(), paramvalue, ROOT.gPad.GetUymax())

    c3.cd(2)
    paramdistriframe = MC_study.plotParam(parameter, rf.FrameBins(200), rf.FrameRange(0., 10.))
    paramdistriframe.SetTitle('MC distri '+paramname)
    paramdistriframe.Draw()
    #paramdistriframe.getHist().Fit('gaus', 'QEM')
    ROOT.gPad.Update()
    paramline.DrawLine(paramvalue, ROOT.gPad.GetUymin(), paramvalue, ROOT.gPad.GetUymax())

    c3.cd(3)
    parampullframe = MC_study.plotPull(parameter, rf.FrameBins(100), rf.FrameRange(-5, 5), rf.FitGauss(ROOT.kTRUE))
    parampullframe.SetTitle('MC pull distri '+paramname)
    parampullframe.Draw()
    ROOT.gPad.Update()
    zeroline.DrawLine(0, ROOT.gPad.GetUymin(), 0, ROOT.gPad.GetUymax())

    c3.cd(4)
    MCnllframe = MC_study.plotNLL()
    MCnllframe.Draw()
    MCnllframe.getHist().Fit('gaus', 'QEM')
    ROOT.gPad.Update()
    nll_line.SetLineStyle(1)
    nll_line.DrawLine(nllvalue, ROOT.gPad.GetUymin(), nllvalue, ROOT.gPad.GetUymax())
    paramline.DrawLine(paramvalue, ROOT.gPad.GetUymin(), paramvalue, ROOT.gPad.GetUymax())


# Alternative calculation of 90% C.L. limit on cross section using Monte Carlo statistics
if MC_sets:
    N_sig_list = []
    for i in range(MC_sets):
        value = MC_study.fitParams(i).find('N_signal').getVal()
        N_sig_list.append(value)
    N_sig_list.sort()
    limit = N_sig_list[int(0.9*MC_sets)]
    result_overview['wimp events limit MC'] = limit
    print '-'*80
    print "N_signal limits with 2 different methods:"
    print '{0:30} {1:10}'.format('from NLL curve:', result_overview['wimp events limit'])
    print '{0:30} {1:10}'.format('from 90% from MC toy set fits:', result_overview['wimp events limit MC'])


if SavePlots:
    c1.SaveAs('%iGeV_PDFs.png'%wimp_mass)
    if MC_sets:
        c2.SaveAs('%iGeV_param-stats.png'%wimp_mass)
        c3.SaveAs('%iGeV_signal-stats.png'%wimp_mass)
