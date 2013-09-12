#!/usr/bin/env python

####################################################################################################
##
## Collection of python functions and ROOT functions
## Lukas Hehn, 2013
##
####################################################################################################

import datetime as dt
import calendar
import numpy as np
import couchdbkit
import ROOT

from couchdbkit import Server, Database
from ROOT import TH1F, TH2F, TF1, TMath, TGraph
from parameters import *



# global
global inpath; inpath = 'Run12PeriodInformation/'


# date and time binning functions
March1st = dt.datetime(2009, 3, 1, 0, 0, 0)
OffsetJules = dt.timedelta(hours=-1)
global dayzero; DayZero = March1st + OffsetJules


def unixtime_to_year(unixtime):
    date = dt.datetime.utcfromtimestamp(unixtime)
    delta = date - DayZero
    return (delta.days/365.24)+(delta.seconds/(365.24*86400))


def unixtime_to_day(unixtime):
    date = dt.datetime.utcfromtimestamp(unixtime)
    delta = date - DayZero
    return delta.days


def unixtime_to_dayhour(unixtime):
    date = dt.datetime.utcfromtimestamp(unixtime)
    delta = date - DayZero
    return delta.days + delta.seconds/86400.


def unixtime_to_hour(unixtime):
    date = dt.datetime.utcfromtimestamp(unixtime)
    delta = date - DayZero
    return float(delta.seconds)/3600


def unixtime_to_date(unixtime):
    date = dt.datetime.utcfromtimestamp(unixtime)
    return date


def day_to_date(days):
    delta = dt.timedelta(days)
    date = DayZero + delta
    return date


def date_to_day(date):
    delta = date - DayZero
    return delta.days


def day_to_unixtime(day):
    date = day_to_date(day)
    unixtime = calendar.timegm(date.timetuple())
    return unixtime

def date_to_unixtime(date):
    unixtime = calendar.timegm(date.timetuple())
    return unixtime


# root functions
# trigger efficiency from DAQ trigger threshold and resolution on heat channel
TRIGGER_EFFICIENCY = TF1('trigger_efficiency', '0.5*(1+ROOT::Math::erf(((x-[0])/([1]*sqrt(2)))))', 0, 25)
TRIGGER_EFFICIENCY.SetNpx(2500)
TRIGGER_EFFICIENCY.SetParName(0, 'Threshold')
TRIGGER_EFFICIENCY.SetParName(1, 'Resolution')
TRIGGER_EFFICIENCY.SetTitle('Trigger Efficiency;E_{Heat} (keV);Efficiency')


# measured fiducial efficiency (division of neutron calibration data histogram with/without fiducial cut)
def fiducial_efficiency_func(x, par):
    xx =x[0]
    if xx <= par[2]:
        f = 0.
    elif xx > par[2]:
        f = par[0]*(1-TMath.exp(par[1]*(xx-par[2])))
    return f
FIDUCIAL_EFFICIENCY = TF1('fiducial_efficiency', fiducial_efficiency_func, 0, 25, 3)
FIDUCIAL_EFFICIENCY.SetNpx(2500)
FIDUCIAL_EFFICIENCY.SetParName(0, 'Maximum')
FIDUCIAL_EFFICIENCY.SetParName(1, 'Slope')
FIDUCIAL_EFFICIENCY.SetParName(2, 'Cutoff')
FIDUCIAL_EFFICIENCY.SetTitle('Fiducial Efficiency;E_{ion} (keVee);Efficiency')


# Lindhard quenching relation (nuclear recoil)
LINDHARD_QUENCHING = TF1('lindhard_quenching', '[0]*(x^[1])', 0, 25)
LINDHARD_QUENCHING.SetNpx(2500)
LINDHARD_QUENCHING.SetParName(0, 'a')
LINDHARD_QUENCHING.SetParName(1, 'b')
LINDHARD_QUENCHING.FixParameter(0, 0.16)
LINDHARD_QUENCHING.FixParameter(1, 0.18)
LINDHARD_QUENCHING.SetTitle('Lindhard Quenching for Nuclear Recoils;E_{Recoil} [keV];Q(E_{Rec})')


# Recoil energy estimator for nuclear recoils
RECOIL_ESTIMATOR = TF1('recoil_energy_estimator', '(x/(1+[0]/[1]))*(1+[0]/[1]*0.16*x^0.18)', 0, 25)
RECOIL_ESTIMATOR.SetNpx(2500)
RECOIL_ESTIMATOR.SetParName(0, 'Voltage')
RECOIL_ESTIMATOR.SetParName(1, 'Creation Potential')
RECOIL_ESTIMATOR.FixParameter(1, 3.0)
RECOIL_ESTIMATOR.SetTitle('E_{Rec} estimator from E_{Heat};E_{Rec} [keVnr];E_{Heat} [keVee]')


# centroids of ER and NR band
def ER_CENTROID_REAL_FUNC(x, par):
    xx =x[0]
    f = xx*(1+(0.16*xx**(0.18+0j))*(par[0]/3))/(1+par[0]/3)
    return f.real
ER_CENTROID_REAL = TF1("ER_CENTROID_REAL", ER_CENTROID_REAL_FUNC, -5, 25, 1)
ER_CENTROID_REAL.SetNpx(3000)
ER_CENTROID_REAL.SetParameter(0, 6.4)
ER_CENTROID_REAL.SetTitle('Electron Recoil Centroid (Real Part Only);E_{Rec} [keVnr];E_{Heat} [keVee]')


NR_CENTROID = TF1('NR_CENTROID', '0.16*x^1.18', 0, 25)
NR_CENTROID.SetNpx(2500)
NR_CENTROID.SetTitle('Nuclear Recoil Centroid;E_{Rec} [keVnr];E_{Heat} [keVee]')


ER_CENTROID = RECOIL_ESTIMATOR.Clone('ER_CENTROID')
ER_CENTROID.SetTitle('Electron Recoil Centroid;E_{Rec} [keVnr];E_{Heat} [keVee]')


GAMMA_CUT = TF1('gamma_cut', 'x*(1+(0.16*x^0.18)*([0]/[1]))/(1+[0]/[1])-[4]*sqrt([2]**2+([3]*(1+(0.16*x^0.18)*([0]/3))/(1+([0]/3)))**2)', 0, 25)
GAMMA_CUT.SetNpx(2500)
GAMMA_CUT.SetParName(0, 'Voltage')
GAMMA_CUT.SetParName(1, 'Creation Potential')
GAMMA_CUT.FixParameter(1, 3.0)
GAMMA_CUT.SetParName(2, 'Sigma_Ion')
GAMMA_CUT.SetParName(3, 'Sigma_Rec')
GAMMA_CUT.SetParName(4, 'Gamma Cut')
GAMMA_CUT.SetTitle('Cut on ER Centroid;E_{Rec} [keVnr];E_{Heat} [keVee]')


def fwhm_rec_from_heat(FWHM_heat, voltage, Erec):
    Q = LINDHARD_QUENCHING.Eval(Erec)
    FWHM_rec = FWHM_heat * ((1+voltage/3.0)/(1+1.18*Q*voltage/3.0))
    return FWHM_rec


def recoil_energy_estimator(energy_ee, voltage):
    function = RECOIL_ESTIMATOR
    function.SetParameter(0, voltage)
    function.FixParameter(1, 3.0) #electron-hole creation potential
    energy_recoil = function.GetX(energy_ee)
    return energy_recoil


def voltage_flag_list(Detector, Runtype):
    RunName, VoltageFlag = [], []
    infile = open(inpath+Detector+'_voltage_'+Runtype+'.txt', 'r')
    for line in infile:
        runname, voltage = line.split()
        RunName.append(str(runname))
        VoltageFlag.append(int(voltage))
    infile.close()
    return [RunName, VoltageFlag]


def voltage_flag(InList, RunName):
    try:
        flag = InList[1][InList[0].index(RunName)]
        return flag
    except ValueError:
        print "error finding run", RunName
        return 0


def samba_number(Detector):
    if Detector in ['ID5', 'ID403', 'FID401']:
        return 'S1'
    elif Detector in ['ID3', 'ID401', 'FID402']:
        return 'S2'
    elif Detector in ['ID4', 'ID6', 'ID404']:
        return 'S3'
    elif Detector in ['ID2', 'ID402', 'ID405']:
        return 'S4'


def era_constants(Detector):
    s = couchdbkit.Server('https://edelweissuser:edwdbw1mp@edelweiss.cloudant.com')
    db = s['analysis']
    vr = db.view('constants/run12', include_docs=True)
    constantDoc = {}
    for row in vr:
        constantDoc[row['key']] = row['doc']
    return constantDoc[Detector]


def wimp_spectrum_eric(WIMP_MASS):
    basedir="/kalinka/home/hehn/PhD/LowMassEric/"
    filename=basedir+"wimpsignal_M"+str(WIMP_MASS)+".txt"

    Energies, Rates = [], []

    print "reading WIMP spectrum from:", filename

    infile = open(filename, 'r')
    for line in infile:
        energy, rate = line.split()
        Energies.append(float(energy))
        Rates.append(float(rate))
    infile.close()

    Energybins = np.array(Energies, dtype=np.float)

    Hist = TH1F('wimp_spectrum_%sGeV'%WIMP_MASS, 'Spectrum Eric %sGeV;E_{recoil} [keV];Rate [evts/kg/day/0.02keV]'%WIMP_MASS, Energybins.size-1, Energybins.flatten('C'))

    for i in range(len(Energies)):
        Hist.Fill(Energies[i], Rates[i])

    return Hist


def wimp_signal(WIMP_MASS, SIGMA_ION, SIGMA_REC,
                rec_bins=200, rec_min=0., rec_max=20., ion_bins=100, ion_min=0., ion_max=10.):
    spectrum = wimp_spectrum_eric(WIMP_MASS)

    denom_i=1./(2*SIGMA_ION**2)
    denom_r=1./(2*SIGMA_REC**2)
    hist = TH2F('wimp_signal_%sGeV'%WIMP_MASS,
              'WIMP signal %2iGeV;E_{rec} (keVnr);E_{ion} (keVee);Rate (cts/kg*day)'%WIMP_MASS,
              rec_bins, rec_min, rec_max, ion_bins, ion_min, ion_max)
    for recbin in range(1, hist.GetNbinsX()+1):
        Erec = hist.GetXaxis().GetBinCenter(recbin)
        for ionbin in range(1, hist.GetNbinsY()+1):
            Eion = hist.GetYaxis().GetBinCenter(ionbin)
            summe=0
            for specbin in range(1, spectrum.GetNbinsX()+1):
                ErecSpec = spectrum.GetXaxis().GetBinCenter(specbin)
                Q = LINDHARD_QUENCHING.Eval(ErecSpec)
                wimprate = spectrum.GetBinContent(specbin)
                tutu = denom_r * pow((Erec - ErecSpec), 2)
                kernel = TMath.exp(-tutu-denom_i * (Eion-Q*ErecSpec)**2)
                summe += (kernel * wimprate)
            hist.SetBinContent(recbin, ionbin, 0.02002*summe/(2*3.141592*SIGMA_REC*SIGMA_ION))
    return hist


def flat_gamma_bckgd(SIGMA_ION, SIGMA_REC, rec_bins=200, rec_min=0., rec_max=20., ion_bins=100, ion_min=0., ion_max=10.):
    denom_i=1./(2*SIGMA_ION**2)
    denom_r=1./(2*SIGMA_REC**2)
    spectrum = TH1F('flat_gamma_spectrum', 'flat gamma spectrum',
                    rec_bins, rec_min, rec_max)  # this spectrum could actually contain a varying bckgd rate!
    hist = TH2F('flat_gamma_bckgd', 'Flat Gamma Bckgd;E_{rec} (keVnr);E_{ion} (keVee);Rate (cts/kg*day)',
                rec_bins, rec_min, rec_max, ion_bins, ion_min, ion_max)
    for recbin in range(1, hist.GetNbinsX()+1):
        Erec = hist.GetXaxis().GetBinCenter(recbin)
        for ionbin in range(1, hist.GetNbinsY()+1):
            Eion = hist.GetYaxis().GetBinCenter(ionbin)
            summe=0
            for specbin in range(1, spectrum.GetNbinsX()+1):
                ErecSpec = spectrum.GetXaxis().GetBinCenter(specbin)
                EionSpec = ER_CENTROID_REAL.Eval(ErecSpec)
                tutu = denom_r * (Erec - ErecSpec)**2
                kernel = TMath.exp(-tutu-denom_i * (Eion - EionSpec)**2)
                summe += kernel
            hist.SetBinContent(recbin, ionbin, summe)
    return hist


def simple_efficiency(detector_name, e_thresh, sigma_rec,
                      rec_bins=200, rec_min=0., rec_max=20., ion_bins=100, ion_min=0., ion_max=10.):
    hist = TH2F('total_efficiency', 'Total Efficiency;E_{rec} (keVnr);E_{ion} (keVee);Efficiency',
                rec_bins, rec_min, rec_max, ion_bins, ion_min, ion_max)

    trigger_eff = TF1('trigger_efficiency', '0.5*(1+ROOT::Math::erf(((x-[0])/([1]*sqrt(2)))))', rec_min, rec_max)
    trigger_eff.FixParameter(0, e_thresh)
    trigger_eff.FixParameter(1, sigma_rec)
    trigger_eff_points = int((rec_max - rec_min) * 100)  # 100 points per keV
    trigger_eff.SetNpx(trigger_eff_points)

    fiducial_eff = TF1('fiducial_efficiency', fiducial_efficiency_func, ion_min, ion_max, 3)
    fiducial_eff_points = int((ion_max - ion_min) * 100)  # 100 points per keV
    fiducial_eff.SetNpx(fiducial_eff_points)

    if detector_name == 'ID2':
        fiducial_eff.FixParameter(0, 0.99)
        fiducial_eff.FixParameter(1, -1.73)
        fiducial_eff.FixParameter(2, 1.74)
    elif detector_name == 'ID3':
        fiducial_eff.FixParameter(0, 0.947)
        fiducial_eff.FixParameter(1, -1.876)
        fiducial_eff.FixParameter(2, 1.247)
    elif detector_name == 'ID6':
        fiducial_eff.FixParameter(0, 0.96)
        fiducial_eff.FixParameter(1, -1.23)
        fiducial_eff.FixParameter(2, 1.43)
    elif detector_name == 'ID401':
        fiducial_eff.FixParameter(0, 0.97)
        fiducial_eff.FixParameter(1, -2.386)
        fiducial_eff.FixParameter(2, 1.394)
    elif detector_name == 'ID404':
        fiducial_eff.FixParameter(0, 0.97)
        fiducial_eff.FixParameter(1, -6.34)
        fiducial_eff.FixParameter(2, 1.90)
    else:
        print "Detector not implemented!"
        return None

    for recbin in range(1, hist.GetNbinsX()+1):
        Erec = hist.GetXaxis().GetBinCenter(recbin)
        eff_rec = trigger_eff.Eval(Erec)
        for ionbin in range(1, hist.GetNbinsY()+1):
            Eion = hist.GetYaxis().GetBinCenter(ionbin)
            eff_ion = fiducial_eff.Eval(Eion)

            if eff_rec >= 0. and eff_ion >= 0.:
                eff_total = eff_rec * eff_ion
            else:
                eff_total = 0.
            hist.SetBinContent(recbin, ionbin, eff_total)
            hist.SetBinError(recbin, ionbin, 0.)
    return hist


def tgraph_from_dataset(dataset):
    entries = dataset.numEntries()
    tgraph = TGraph()
    for event in range(entries):
        time = dataset.get(event).getRealValue('time')
        rec = dataset.get(event).getRealValue('rec')
        ion = dataset.get(event).getRealValue('ion')
        tgraph.SetPoint(event, rec, ion)
    return tgraph


def cut_wimp_signal(hist):
    for recbin in range(1, hist.GetNbinsX()+1):
        for ionbin in range(1, hist.GetNbinsY()+1):
            rate = hist.GetBinContent(recbin, ionbin)
            if rate < 0.1e-3:
                hist.SetBinContent(recbin, ionbin, 0.0)
    return "Histogram cut at 0.1"
