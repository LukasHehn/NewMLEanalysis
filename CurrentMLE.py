#!/usr/bin/env python

######################################################
#
# calculation of 90% C.L. cross section limits using Maximum Likelihood methods based on ROOT.RooFit
# Lukas Hehn, 2013
#
######################################################

import ROOT
import Functions
#from DetectorClass import *


# define maximal energies
EnergyIonMax = 10.
EnergyRecMax = 20.


#switches and input parameters to control script
wimp_mass = 8 #set wimp mass or switch signal of entirely (with False)
MC_sets = 10 #int(1e4) #set number of MC simulations: 0 means none at all
SavePlots = False #flag decides whether plots are saved

DataFile = '/kalinka/home/hehn/PhD/LowMassEric/ID3_eventlist.txt'


# definition of observables
ion = ROOT.RooRealVar('ion', 'E_{ion}', 0., EnergyIonMax, 'keV_{ee}')
rec = ROOT.RooRealVar('rec', 'E_{rec}', 0., EnergyRecMax, 'keV_{nr}')


# detector specific parameters
voltage = ROOT.RooRealVar('voltage', 'applied voltage', 6.4)

E_thresh = 3.874
FWHM_heat = 0.82 #Eric@10keV
FWHM_rec = Functions.RecoilResolutionFromHeat(FWHM_heat, voltage.getVal(), 10)
sigma_rec = ROOT.RooRealVar('sigma_rec', 'recoil energy resolution', FWHM_rec/2.35)

FWHM_ion = 0.72 #Eric
sigma_ion = ROOT.RooRealVar('sigma_ion', 'ionization energy resolution', FWHM_ion/2.35)


# detector efficiency
total_efficiency = Functions.Simple2DEfficiencyID3(E_thresh, FWHM_rec/2.35)
total_efficiency_datahist = ROOT.RooDataHist('total_efficiency_datahist', 'total_efficiency_datahist', ROOT.RooArgList(rec, ion), total_efficiency)
total_efficiency_pdf = ROOT.RooHistPdf('total_efficiency_pdf', 'total_efficiency_pdf', ROOT.RooArgSet(rec, ion), total_efficiency_datahist)


# dataset
realdata = ROOT.RooDataSet.read(DataFile, ROOT.RooArgList(rec, ion))
realdata_scatter = realdata.createHistogram(rec, ion, int(EnergyRecMax*10), int(EnergyIonMax*10))
realdata_graph = Functions.TGraphFromDataSet(realdata)
events = int(realdata.numEntries())


# -----------------------------------------------------------------------------------------
# definition of gamma peaks
ion_scaling = ROOT.RooRealVar('ion_scaling', 'scaling factor ionization energy', 1.0, 0.95, 1.05)
rec_scaling = ROOT.RooRealVar('rec_scaling', 'scaling factor recoil energy', 1.0, 0.95, 1.05)
ion_scaling.setConstant(ROOT.kTRUE)
rec_scaling.setConstant(ROOT.kTRUE)
ER_centroid.SetParameter(0, 6.4) #ER_centroid used for calculation of peak position in Erec

V49_ion_energy = ROOT.RooRealVar('V49_ion_energy', 'V49 peak ion energy', 4.97)
V49_ion_pos = ROOT.RooFormulaVar('V49_ion_pos', '@0*@1', ROOT.RooArgList(V49_ion_energy, ion_scaling))
V49_ion_pdf = ROOT.RooGaussian('V49_ion_pdf', 'V49 peak gauss pdf in ion', ion, V49_ion_pos, sigma_ion)
V49_rec_energy = ROOT.RooRealVar('V49_rec_energy', 'recoil energy V49 peak', ER_centroid.GetX(V49_ion_energy.getVal()))
V49_rec_pos = ROOT.RooFormulaVar('V49_rec_pos', '@0*@1', ROOT.RooArgList(V49_rec_energy, rec_scaling))
V49_rec_pdf = ROOT.RooGaussian('V49_rec_pdf', 'V49 peak pdf in recoil energy', rec, V49_rec_pos, sigma_rec)
V49_pdf = ROOT.RooProdPdf('V49_pdf', 'V49 peak pdf', V49_ion_pdf, V49_rec_pdf)
V49_pdf_eff = ROOT.RooProdPdf('V49_pdf_eff', 'eff corr V49 peak pdf', V49_pdf, total_efficiency_pdf)
N_V49 = ROOT.RooRealVar('N_V49', 'evts of V49 peak (4.97keV)', 16., 0., events)


Cr51_ion_energy = ROOT.RooRealVar('Cr51_ion_energy', 'Cr51 peak ion energy', 5.46)
Cr51_ion_pos = ROOT.RooFormulaVar('Cr51_ion_pos', '@0*@1', ROOT.RooArgList(Cr51_ion_energy, ion_scaling))
Cr51_ion = ROOT.RooGaussian('Cr51_ion_pdf', 'Cr51 peak gauss pdf with shifted mean', ion, Cr51_ion_pos, sigma_ion)
Cr51_rec_energy = ROOT.RooRealVar('Cr51_rec_energy', 'v_rec_energy', ER_centroid.GetX(Cr51_ion_energy.getVal()))
Cr51_rec_pos = ROOT.RooFormulaVar('Cr51_rec_pos', '@0*@1', ROOT.RooArgList(Cr51_rec_energy, rec_scaling))
Cr51_rec = ROOT.RooGaussian('Cr51_rec_pdf', 'Cr51_rec_pdf with shifted mean', rec, Cr51_rec_pos, sigma_rec)
Cr51_pdf = ROOT.RooProdPdf('Cr51_pdf', 'Cr51 peak pdf', Cr51_ion, Cr51_rec)
N_Cr51 = ROOT.RooRealVar('N_Cr51', 'evts of 51Cr peak (5.46keV)', 11., 0., events)


Mn54_ion_energy = ROOT.RooRealVar('Mn54_ion_energy', 'Mn54_ion_energy', 5.99)
Mn54_ion_pos = ROOT.RooFormulaVar('Mn54_ion_pos', '@0*@1', ROOT.RooArgList(Mn54_ion_energy, ion_scaling))
Mn54_ion = ROOT.RooGaussian('Mn54_ion_pdf', 'Mn54_ion_pdf with shifted mean', ion, Mn54_ion_pos, sigma_ion)
Mn54_rec_energy = ROOT.RooRealVar('Mn54_rec_energy', 'Mn54_rec_energy', ER_centroid.GetX(Mn54_ion_energy.getVal()))
Mn54_rec_pos = ROOT.RooFormulaVar('Mn54_rec_pos', '@0*@1', ROOT.RooArgList(Mn54_rec_energy, rec_scaling))
Mn54_rec = ROOT.RooGaussian('Mn54_rec_pdf', 'Mn54_rec_pdf with shifted mean', rec, Mn54_rec_pos, sigma_rec)
Mn54_pdf = ROOT.RooProdPdf('Mn54_pdf', 'Mn54 peak pdf', Mn54_ion, Mn54_rec)
N_Mn54 = ROOT.RooRealVar('N_Mn54', 'evts of 54Mn peak (5.99keV)', 4., 0., events)


Fe55_ion_energy = ROOT.RooRealVar('Fe55_ion_energy', 'Fe55_ion_energy', 6.54)
Fe55_ion_pos = ROOT.RooFormulaVar('Fe55_ion_pos', '@0*@1', ROOT.RooArgList(Fe55_ion_energy, ion_scaling))
Fe55_ion = ROOT.RooGaussian('Fe55_ion_pdf', 'Fe55_ion_pdf with shifted mean', ion, Fe55_ion_pos, sigma_ion)
Fe55_rec_energy = ROOT.RooRealVar('Fe55_rec_energy', 'Fe55_rec_energy', ER_centroid.GetX(Fe55_ion_energy.getVal()))
Fe55_rec_pos = ROOT.RooFormulaVar('Fe55_rec_pos', '@0*@1', ROOT.RooArgList(Fe55_rec_energy, rec_scaling))
Fe55_rec = ROOT.RooGaussian('Fe55_rec_pdf', 'Fe55_rec_pdf with shifted mean', rec, Fe55_rec_pos, sigma_rec)
Fe55_pdf = ROOT.RooProdPdf('Fe55_pdf', 'Fe55 peak pdf', Fe55_ion, Fe55_rec)
N_Fe55 = ROOT.RooRealVar('N_Fe55', 'evts of 55Fe peak (6.54keV)', 31., 0., events)


Co57_ion_energy = ROOT.RooRealVar('Co57_ion_energy', 'Co57_ion_energy', 7.11)
Co57_ion_pos = ROOT.RooFormulaVar('Co57_ion_pos', '@0*@1', ROOT.RooArgList(Co57_ion_energy, ion_scaling))
Co57_ion = ROOT.RooGaussian('Co57_ion_pdf', 'Co57_ion_pdf with shifted mean', ion, Co57_ion_pos, sigma_ion)
Co57_rec_energy = ROOT.RooRealVar('Co57_rec_energy', 'Co57_rec_energy', ER_centroid.GetX(Co57_ion_energy.getVal()))
Co57_rec_pos = ROOT.RooFormulaVar('Co57_rec_pos', '@0*@1', ROOT.RooArgList(Co57_rec_energy, rec_scaling))
Co57_rec = ROOT.RooGaussian('Co57_rec_pdf', 'Co57_rec_pdf with shifted mean', rec, Co57_rec_pos, sigma_rec)
Co57_pdf = ROOT.RooProdPdf('Co57_pdf', 'Co57 peak pdf', Co57_ion, Co57_rec)
N_Co57 = ROOT.RooRealVar('N_Co57', 'evts of 57Co peak (7.11keV)', 2., 0., events)


Zn65_ion_energy = ROOT.RooRealVar('Zn65_ion_energy', 'Zn65_ion_energy', 8.98)
Zn65_ion_pos = ROOT.RooFormulaVar('Zn65_ion_pos', '@0*@1', ROOT.RooArgList(Zn65_ion_energy, ion_scaling))
Zn65_ion = ROOT.RooGaussian('Zn65_ion_pdf', 'Zn65_ion_pdf with shifted mean', ion, Zn65_ion_pos, sigma_ion)
Zn65_rec_energy = ROOT.RooRealVar('Zn65_rec_energy', 'Zn65_rec_energy', ER_centroid.GetX(Zn65_ion_energy.getVal()))
Zn65_rec_pos = ROOT.RooFormulaVar('Zn65_rec_pos', '@0*@1', ROOT.RooArgList(Zn65_rec_energy, rec_scaling))
Zn65_rec = ROOT.RooGaussian('Zn65_rec_pdf', 'Zn65_rec_pdf with shifted mean', rec, Zn65_rec_pos, sigma_rec)
Zn65_pdf = ROOT.RooProdPdf('Zn65_pdf', 'Zn65 peak pdf', Zn65_ion, Zn65_rec)
N_Zn65 = ROOT.RooRealVar('N_Zn65', 'evts of 65Zn peak (8.98keV)', 110., 0., events)


Ga68_ion_energy = ROOT.RooRealVar('Ga68_ion_energy', 'Ga68_ion_energy', 9.66)
Ga68_ion_pos = ROOT.RooFormulaVar('Ga68_ion_pos', '@0*@1', ROOT.RooArgList(Ga68_ion_energy, ion_scaling))
Ga68_ion = ROOT.RooGaussian('Ga68_ion_pdf', 'Ga68_ion_pdf with shifted mean', ion, Ga68_ion_pos, sigma_ion)
Ga68_rec_energy = ROOT.RooRealVar('Ga68_rec_energy', 'Ga68_rec_energy', ER_centroid.GetX(Ga68_ion_energy.getVal()))
Ga68_rec_pos = ROOT.RooFormulaVar('Ga68_rec_pos', '@0*@1', ROOT.RooArgList(Ga68_rec_energy, rec_scaling))
Ga68_rec = ROOT.RooGaussian('Ga68_rec_pdf', 'Ga68_rec_pdf with shifted mean', rec, Ga68_rec_pos, sigma_rec)
Ga68_pdf = ROOT.RooProdPdf('Ga68_pdf', 'Ga68 peak pdf', Ga68_ion, Ga68_rec)
N_Ga68 = ROOT.RooRealVar('N_Ga68', 'evts of 68Ga peak (9.66keV)', 32., 0., events)


Ge68_ion_energy = ROOT.RooRealVar('Ge68_ion_energy', 'Ge68_ion_energy', 10.37)
Ge68_ion_pos = ROOT.RooFormulaVar('Ge68_ion_pos', '@0*@1', ROOT.RooArgList(Ge68_ion_energy, ion_scaling))
Ge68_ion = ROOT.RooGaussian('Ge68_ion_pdf', 'Ge68_ion_pdf with shifted mean', ion, Ge68_ion_pos, sigma_ion)
Ge68_rec_energy = ROOT.RooRealVar('Ge68_rec_energy', 'Ge68_rec_energy', ER_centroid.GetX(Ge68_ion_energy.getVal()))
Ge68_rec_pos = ROOT.RooFormulaVar('Ge68_rec_pos', '@0*@1', ROOT.RooArgList(Ge68_rec_energy, rec_scaling))
Ge68_rec = ROOT.RooGaussian('Ge68_rec_pdf', 'Ge68 peak pdf in rec', rec, Ge68_rec_pos, sigma_rec)
Ge68_pdf = ROOT.RooProdPdf('Ge68_pdf', 'Ge68 peak pdf', Ge68_ion, Ge68_rec)
N_Ge68 = ROOT.RooRealVar('N_Ge68', 'evts of 68Ge peak (10.37keV)', 0., 0., events)
# -----------------------------------------------------------------------------------------

# wimp signal
if wimp_mass:
    # read in WIMP spectrum
    signal_hist = Functions.WimpSignal2DEric(wimp_mass, sigma_ion.getVal(), sigma_rec.getVal())
    signal_hist.Multiply(total_efficiency)
    signal_datahist = ROOT.RooDataHist('signal_datahist', 'signal_datahist', 
                                       ROOT.RooArgList(rec, ion), signal_hist)
    signal_pdf = ROOT.RooHistPdf('signal_pdf', 'signal_pdf', 
                                 ROOT.RooArgSet(rec, ion), signal_datahist)
    N_signal = ROOT.RooRealVar('N_signal', 'WIMP signal events', 0., 0., 10.)


# flat gamma background component
flat_gamma_bckgd_hist = Functions.FlatGammaBckgd2DEric(sigma_ion.getVal(), sigma_rec.getVal())
flat_gamma_bckgd_hist.Multiply(total_efficiency)
flat_gamma_bckgd_datahist = ROOT.RooDataHist('flat_gamma_bckgd_datahist', 'flat_gamma_bckgd_datahist', 
                                             ROOT.RooArgList(rec, ion), flat_gamma_bckgd_hist)
flat_gamma_bckgd_pdf = ROOT.RooHistPdf('flat_gamma_bckgd_pdf', 'flat_gamma_bckgd_pdf', 
                                       ROOT.RooArgSet(rec, ion), flat_gamma_bckgd_datahist)
N_flat = ROOT.RooRealVar('N_flat', 'bckgd events', 70., 0., events)


if wimp_mass:
    # definition of pdf to be fitted
    final_pdf = ROOT.RooAddPdf('final_pdf', 'final_pdf', 
                             ROOT.RooArgList(signal_pdf, flat_gamma_bckgd_pdf, V49_pdf, Cr51_pdf, Mn54_pdf, Fe55_pdf, Co57_pdf, Zn65_pdf, Ga68_pdf), 
                             ROOT.RooArgList(N_signal, N_flat, N_V49, N_Cr51, N_Mn54, N_Fe55, N_Co57, N_Zn65, N_Ga68))
else:
    final_pdf = ROOT.RooAddPdf('final_pdf', 'final_pdf', 
                               ROOT.RooArgList(flat_gamma_bckgd_pdf, V49_pdf, Cr51_pdf, Mn54_pdf, Fe55_pdf, Co57_pdf, Zn65_pdf, Ga68_pdf, Ge68_pdf), 
                               ROOT.RooArgList(N_flat, N_V49, N_Cr51, N_Mn54, N_Fe55, N_Co57, N_Zn65, N_Ga68, N_Ge68))


# manual mode
nll = ROOT.RooNLLVar('nll', 'nll', final_pdf, realdata, 
                     ROOT.RooFit.Extended(ROOT.kTRUE), 
                     ROOT.RooFit.PrintEvalErrors(2), 
                     ROOT.RooFit.Verbose(ROOT.kFALSE))
minuit = ROOT.RooMinuit(nll)
minuit.migrad() #find minimum
minuit.hesse() #symmetric errors
minuit.minos() #asymmetric errors
FitResult = minuit.save('realfit', 'fit to real data')
ndf = FitResult.floatParsFinal().getSize()

FitResult.Print('v')


# calculate and print cross section limit for this wimp mass:
if wimp_mass:
    signal_parameter = FitResult.floatParsFinal().find('N_signal')
    wimp_events = signal_parameter.getVal()
    error_low = signal_parameter.getAsymErrorLo()
    error_high = signal_parameter.getAsymErrorHi()
    wimp_events_limit = wimp_events + 1.28 * error_high
    livetime = 197. #wimp search livetime in days
    mass = 0.160 #detector mass in kg
    rate = signal_hist.Integral('WIDTH')
    signal_events = rate * livetime * mass * 1e6
    sigma_limit = wimp_events_limit / signal_events
    
    print '{0:9} | {1:10} | {2:8} | {3:8} | {4:10} | {5:10}'.format('wimp_mass', 'Rate', 'N_signal', 'ErrorLow', 'ErrorHigh', 'XS-limit [pb]')
    print "-----------------------------------------------------------------------------------"
    print '{0:9} | {1:10} | {2:8} | {3:8} | {4:10} | {5:10}'.format(wimp_mass, rate, wimp_events, error_low, error_high, sigma_limit)



# histogram
final_hist = final_pdf.createHistogram('final_hist', rec, ROOT.RooFit.Binning(int((rec.getMax()-rec.getMin())*10)), ROOT.RooFit.YVar(ion, ROOT.RooFit.Binning(int((ion.getMax()-ion.getMin())*10))))


# RooFit frames
recbins = int((rec.getMax()-rec.getMin())*2)
ionbins = int((ion.getMax()-ion.getMin())*5)


ionframe = ion.frame()
ionframe.SetTitle('Projection in E_{ion}')
realdata.plotOn(ionframe, ROOT.RooFit.Name('data'), ROOT.RooFit.Binning(ionbins), ROOT.RooFit.MarkerSize(1.0))
final_pdf.plotOn(ionframe, ROOT.RooFit.Components("flat_gamma_bckgd_pdf"), ROOT.RooFit.LineColor(kGreen), ROOT.RooFit.LineWidth(2))
final_pdf.plotOn(ionframe, ROOT.RooFit.Name('model'), ROOT.RooFit.LineColor(kBlue), ROOT.RooFit.LineWidth(2))
final_pdf.plotOn(ionframe, ROOT.RooFit.Components("V49_ion_pdf"), ROOT.RooFit.LineColor(kRed), ROOT.RooFit.LineWidth(2), ROOT.RooFit.LineStyle(kDashed))
final_pdf.plotOn(ionframe, ROOT.RooFit.Components("Cr51_ion_pdf"), ROOT.RooFit.LineColor(kRed), ROOT.RooFit.LineWidth(2), ROOT.RooFit.LineStyle(kDashed))
final_pdf.plotOn(ionframe, ROOT.RooFit.Components("Mn54_ion_pdf"), ROOT.RooFit.LineColor(kRed), ROOT.RooFit.LineWidth(2), ROOT.RooFit.LineStyle(kDashed))
final_pdf.plotOn(ionframe, ROOT.RooFit.Components("Fe55_ion_pdf"), ROOT.RooFit.LineColor(kRed), ROOT.RooFit.LineWidth(2), ROOT.RooFit.LineStyle(kDashed))
final_pdf.plotOn(ionframe, ROOT.RooFit.Components("Co57_ion_pdf"), ROOT.RooFit.LineColor(kRed), ROOT.RooFit.LineWidth(2), ROOT.RooFit.LineStyle(kDashed))
final_pdf.plotOn(ionframe, ROOT.RooFit.Components("Zn65_ion_pdf"), ROOT.RooFit.LineColor(kRed), ROOT.RooFit.LineWidth(2), ROOT.RooFit.LineStyle(kDashed))
final_pdf.plotOn(ionframe, ROOT.RooFit.Components("Ge68_ion_pdf"), ROOT.RooFit.LineColor(kRed), ROOT.RooFit.LineWidth(2), ROOT.RooFit.LineStyle(kDashed))
final_pdf.plotOn(ionframe, ROOT.RooFit.Components("Ga68_ion_pdf"), ROOT.RooFit.LineColor(kRed), ROOT.RooFit.LineWidth(2), ROOT.RooFit.LineStyle(kDashed))
final_pdf.plotOn(ionframe, ROOT.RooFit.Components("signal_pdf"), ROOT.RooFit.LineColor(kMagenta), ROOT.RooFit.LineWidth(3))
final_pdf.paramOn(ionframe, ROOT.RooFit.Format('NEU', ROOT.RooFit.AutoPrecision(2)), ROOT.RooFit.Layout(0.1, 0.55, 0.9), ROOT.RooFit.ShowConstants(ROOT.kFALSE))


recframe = rec.frame()
recframe.SetTitle('Projection in E_{rec}')
realdata.plotOn(recframe, ROOT.RooFit.Name("data"), ROOT.RooFit.Binning(recbins), ROOT.RooFit.MarkerColor(kBlack), ROOT.RooFit.MarkerSize(1.0))
final_pdf.plotOn(recframe, ROOT.RooFit.Components("flat_gamma_bckgd_pdf"), ROOT.RooFit.LineColor(kGreen), ROOT.RooFit.LineWidth(2))
final_pdf.plotOn(recframe, ROOT.RooFit.Name('model'), ROOT.RooFit.LineColor(kBlue), ROOT.RooFit.LineWidth(2))
final_pdf.plotOn(recframe, ROOT.RooFit.Components("V49_rec_pdf"), ROOT.RooFit.LineColor(kRed), ROOT.RooFit.LineWidth(2), ROOT.RooFit.LineStyle(kDashed))
final_pdf.plotOn(recframe, ROOT.RooFit.Components("Cr51_rec_pdf"), ROOT.RooFit.LineColor(kRed), ROOT.RooFit.LineWidth(2), ROOT.RooFit.LineStyle(kDashed))
final_pdf.plotOn(recframe, ROOT.RooFit.Components("Mn54_rec_pdf"), ROOT.RooFit.LineColor(kRed), ROOT.RooFit.LineWidth(2), ROOT.RooFit.LineStyle(kDashed))
final_pdf.plotOn(recframe, ROOT.RooFit.Components("Fe55_rec_pdf"), ROOT.RooFit.LineColor(kRed), ROOT.RooFit.LineWidth(2), ROOT.RooFit.LineStyle(kDashed))
final_pdf.plotOn(recframe, ROOT.RooFit.Components("Co57_rec_pdf"), ROOT.RooFit.LineColor(kRed), ROOT.RooFit.LineWidth(2), ROOT.RooFit.LineStyle(kDashed))
final_pdf.plotOn(recframe, ROOT.RooFit.Components("Zn65_rec_pdf"), ROOT.RooFit.LineColor(kRed), ROOT.RooFit.LineWidth(2), ROOT.RooFit.LineStyle(kDashed))
final_pdf.plotOn(recframe, ROOT.RooFit.Components("Ge68_rec_pdf"), ROOT.RooFit.LineColor(kRed), ROOT.RooFit.LineWidth(2), ROOT.RooFit.LineStyle(kDashed))
final_pdf.plotOn(recframe, ROOT.RooFit.Components("Ga68_rec_pdf"), ROOT.RooFit.LineColor(kRed), ROOT.RooFit.LineWidth(2), ROOT.RooFit.LineStyle(kDashed))
final_pdf.plotOn(recframe, ROOT.RooFit.Components("signal_pdf"), ROOT.RooFit.LineColor(kMagenta), ROOT.RooFit.LineWidth(3))


# Plotting of fit results (PDFs and projected spectra)
c1 = ROOT.TCanvas('c1', 'Fit Results', 1000, 750)
c1.Divide(2, 2)
# final pdf and data
c1.cd(1)
c1.cd(1).SetLogz()
final_hist.Draw('COL')
final_hist.SetStats(0)
final_hist.SetTitle('ID3 WIMP search data + best fit PDF')
realdata_graph.SetMarkerColor(kMagenta)
realdata_graph.SetMarkerStyle(kPlus)
realdata_graph.Draw('SAMESP')
ER_centroid.SetLineColor(kBlack)
ER_centroid.SetLineWidth(1)
ER_centroid.DrawCopy('SAME')
NR_centroid.SetLineColor(kBlack)
NR_centroid.SetLineWidth(1)
NR_centroid.DrawCopy('SAME')
# wimp signal and data only
if wimp_mass:
    c1.cd(3)
    signal_hist.Draw('CONT0')
    signal_hist.SetStats(0)
    signal_hist.SetContour(30)
    realdata_graph.Draw('SAMESP')
    ER_centroid.DrawCopy('SAME')
    NR_centroid.DrawCopy('SAME')
# data and fit in ionization energy
c1.cd(2)
ionframe.Draw()
# data and fit in recoil energy
c1.cd(4)
recframe.Draw()


# Monte Carlo statistics output
if MC_sets:
    MC_study = ROOT.RooMCStudy(final_pdf, ROOT.RooArgSet(rec, ion), ROOT.RooFit.Silence(), ROOT.RooFit.Extended(ROOT.kTRUE), ROOT.RooFit.FitOptions(ROOT.RooFit.Save(ROOT.kTRUE)))
    MC_study.generateAndFit(MC_sets)
  
    N_sig_list = []
    for i in range(MC_sets):
        value = MC_study.fitParams(i).find('N_signal').getVal()
        print i, value 
        N_sig_list.append(value)
    N_sig_array = np.array(N_sig_list, float)
    N_sig_array.sort()
    c = N_sig_array[0.9*MC_sets]
    print N_sig_array[0.9*MC_sets]
  
    FitParams = FitResult.floatParsFinal()
    NumFitParams = FitParams.getSize()
  
    ParamNLLFrameList = []
    ParamDistriFrameList = []
    ParamPullFrameList = []
    ParamLineList = []
  
    c2 = ROOT.TCanvas('c2', 'Statistical interpretation (MC study)')
    c2.Divide(1, NumFitParams, 0.001, 0.001)
  
    nll_line = TLine()
    nll_line.SetLineWidth(2)
    nll_line.SetLineStyle(7)
    nll_line.SetLineColor(kBlue)
  
    zeroline = TLine()
    zeroline.SetLineWidth(2)
    zeroline.SetLineStyle(7)
    zeroline.SetLineColor(kBlack)
  
    for i in range(NumFitParams):
        parameter = FitParams[i]
        paramname = parameter.GetName()
        paramvalue = parameter.getVal()
    
        paramline = TLine()
        paramline.SetLineWidth(2)
        paramline.SetLineStyle(1)
        paramline.SetLineColor(kBlue)
        ParamLineList.append(paramline)
    
        pad = c2.cd(i+1)
        pad.Divide(3, 1, 0.001, 0.001)
    
        if paramname == 'N_signal':
            pad.SetFillColor(kYellow-9)
            pad.SetFillStyle(4100)
  
        paramnllframe = parameter.frame(ROOT.RooFit.Title('NLL fit '+paramname))
        #paramnllframe.SetTitle('NLL fit '+paramname)
        nll.plotOn(paramnllframe, ROOT.RooFit.Precision(1e-5), ROOT.RooFit.ShiftToZero())
        paramnllframe.SetMaximum(20.)
        paramnllframe.SetMinimum(0.)
        ParamNLLFrameList.append(paramnllframe)
        pad.cd(1)
        paramnllframe.Draw()
        ROOT.gPad.Update()
        paramline.DrawLine(paramvalue, ROOT.gPad.GetUymin(), paramvalue, ROOT.gPad.GetUymax())
    
        paramdistriframe = MC_study.plotParam(parameter, ROOT.RooFit.Title('MC distri '+paramname))
        #paramdistriframe.SetTitle('MC distri '+paramname)
        ParamDistriFrameList.append(paramdistriframe)
        pad.cd(2)
        paramdistriframe.Draw()
        paramdistriframe.getHist().Fit('gaus', 'QEM')
        ROOT.gPad.Update()
        paramline.DrawLine(paramvalue, ROOT.gPad.GetUymin(), paramvalue, ROOT.gPad.GetUymax())
    
        parampullframe = MC_study.plotPull(parameter, ROOT.RooFit.FitGauss())
        parampullframe.SetTitle('MC pull distri '+paramname)
        ParamPullFrameList.append(parampullframe)
        pad.cd(3)
        parampullframe.Draw()
        ROOT.gPad.Update()
        zeroline.DrawLine(0, ROOT.gPad.GetUymin(), 0, ROOT.gPad.GetUymax())
  
    c2.SetCanvasSize(1200, 2400)
  
    # extra canvas for signal parameter only
    if wimp_mass:
        parameter = FitResult.floatParsFinal().find('N_signal')
        paramname = parameter.GetName()
        paramvalue = parameter.getVal()
    
        nllvalue = nll.getVal()
    
        c3 = ROOT.TCanvas('c3', 'Fit Statistics WIMP Signal', 1000, 750)
        c3.Divide(2, 2)
        c3.cd(1)
        paramnllframe = parameter.frame()
        paramnllframe.SetTitle('NLL fit '+paramname)
        nll.plotOn(paramnllframe, ROOT.RooFit.Precision(1e-5), ROOT.RooFit.ShiftToZero())
        paramnllframe.SetMaximum(20.)
        paramnllframe.SetMinimum(0.)
        paramnllframe.Draw()
        ROOT.gPad.Update()
        paramline.DrawLine(paramvalue, ROOT.gPad.GetUymin(), paramvalue, ROOT.gPad.GetUymax())
    
        c3.cd(2)
        paramdistriframe = MC_study.plotParam(parameter, ROOT.RooFit.FrameBins(200), ROOT.RooFit.FrameRange(0., 10.))
        paramdistriframe.SetTitle('MC distri '+paramname)
        paramdistriframe.Draw()
        #paramdistriframe.getHist().Fit('gaus', 'QEM')
        ROOT.gPad.Update()
        paramline.DrawLine(paramvalue, ROOT.gPad.GetUymin(), paramvalue, ROOT.gPad.GetUymax())
    
        c3.cd(3)
        parampullframe = MC_study.plotPull(parameter, ROOT.RooFit.FrameBins(100), ROOT.RooFit.FrameRange(-5, 5), ROOT.RooFit.FitGauss(ROOT.kTRUE))
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


if SavePlots:
    c1.SaveAs('%iGeV_PDFs.png'%wimp_mass)
    if MC_sets:
        c2.SaveAs('%iGeV_param-stats.png'%wimp_mass)
        c3.SaveAs('%iGeV_signal-stats.png'%wimp_mass)


if MC_sets:
    N_sig_list = []
    for i in range(MC_sets):
        value = cc.fitParams(i).find('N_signal').getVal()
        #print i, value 
        N_sig_list.append(value)
    N_sig_list.sort()
    limit = N_sig_list[int(0.9*MC_sets)]
    print "N_signal limit from 90% of MC toy set fits", limit
  