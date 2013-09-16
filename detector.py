#!/usr/bin/env python

####################################################################################################
##
## The Detector class stores all the relevant baseline information for a detector and Run12 and
## automatically calculated efficiency histograms
## Lukas Hehn, 2013
##
####################################################################################################


import ROOT
import parameters
import functions
import numpy as np

from ROOT import TH1F, TH2F, TH3F, TGraph
from parameters import ENERGY_BINNING as binning


class Detector:
    def __init__(self, name):
        self.fName = name

        self.fMass = self._AssignMass()

        self.fBadDays = []

        self.fEventList = [[], [], []]  # time & energy list of added events

        self.fRunType = 'Bckgd'

        self.fVoltageFlagList = functions.voltage_flag_list(self.fName, self.fRunType)

        self.fERAConstants = functions.era_constants(self.fName)

        self.fGoodUnixPeriods = self._ReadInGoodPeriods()

        self.fGoodUnixPeriodsGap = self._InsertGaps(self.fGoodUnixPeriods)  # fit gaps between periods


        # Variable binning lists for all variables (time, recoil energy, ionization energy)
        self.fTimeBinning = self._CalcYearBins(self.fGoodUnixPeriodsGap)

        self.fEnergyRecBinning = np.array([binning['rec']['min']+i*binning['rec']['binsize'] for i in range(0, binning['rec']['bins']+1)], dtype=np.float)

        self.fEnergyIonBinning = np.array([binning['ion']['min']+i*binning['ion']['binsize'] for i in range(0, binning['ion']['bins']+1)], dtype=np.float)


        # ROOT histograms with baseline values over variable time bins for all channels
        self.fHeatBaseline = TH1F(self.fName+'_heat_baseline',
                                  self.fName+' heat baseline;Time (years);Energy (keV)',
                                  self.fTimeBinning.size-1, self.fTimeBinning.flatten('C')
                                  )

        self.fHeatThreshold = TH1F(self.fName+'_heat_threshold',
                                   self.fName+' heat threshold;Time (years);Energy (keV)',
                                   self.fTimeBinning.size-1, self.fTimeBinning.flatten('C')
                                   )

        self.fFiducialTopBaseline = TH1F(self.fName+'_fiducial_top_baseline',
                                         self.fName+' fiducial top baseline;Time (years);Energy (keV)',
                                         self.fTimeBinning.size-1, self.fTimeBinning.flatten('C')
                                         )

        self.fFiducialBottomBaseline = TH1F(self.fName+'_fiducial_bottom_baseline',
                                            self.fName+' fiducial bottom baseline;Time (years);Energy (keV)',
                                            self.fTimeBinning.size-1, self.fTimeBinning.flatten('C')
                                            )

        self.fVetoTopBaseline = TH1F(self.fName+'_veto_top_baseline',
                                     self.fName+' veto top baseline;Time (years);Energy (keV)',
                                     self.fTimeBinning.size-1, self.fTimeBinning.flatten('C')
                                     )

        self.fVetoBottomBaseline = TH1F(self.fName+'_veto_bottom_baseline',
                                        self.fName+' veto bottom baseline;Time (years);Energy (keV)',
                                        self.fTimeBinning.size-1, self.fTimeBinning.flatten('C')
                                        )

        self.fGuardTopBaseline = TH1F(self.fName+'_guard_top_baseline',
                                      self.fName+' guard top baseline;Time (years);Energy (keV)',
                                      self.fTimeBinning.size-1, self.fTimeBinning.flatten('C')
                                      )

        self.fGuardBottomBaseline = TH1F(self.fName+'_guard_bottom_baseline',
                                         self.fName+' guard bottom baseline;Time (years);Energy (keV)',
                                         self.fTimeBinning.size-1, self.fTimeBinning.flatten('C')
                                         )

        self.fFiducialMeanBaseline = TH1F(self.fName+'_fiducial_mean',
                                          self.fName+' fiducial mean baseline;Time (years);Energy (keV)',
                                         self.fTimeBinning.size-1, self.fTimeBinning.flatten('C')
                                         )

        self.fVoltage = TH1F(self.fName+'_voltage',
                             self.fName+' voltage;Time (years);Voltage (V)',
                             self.fTimeBinning.size-1, self.fTimeBinning.flatten('C')
                             )

        # Histograms with efficiencies
        self.fLivetimeEfficiency = TH1F(self.fName+'_livetime_efficiency',
                                        self.fName+' Livetime Efficiency;Time (years);Efficiency',
                                        self.fTimeBinning.size-1, self.fTimeBinning.flatten('C')
                                        )

        self.fTriggerEfficiency = TH2F(self.fName+'_trigger_efficiency',
                                       self.fName+' Trigger Efficiency;Time (years);E_{Rec} (keVnr);Efficiency',
                                       self.fTimeBinning.size-1, self.fTimeBinning.flatten('C'),
                                       self.fEnergyRecBinning.size-1, self.fEnergyRecBinning.flatten('C')
                                       )

        self.fFiducialEfficiency = TH2F(self.fName+'_fiducial_efficiency',
                                        self.fName+' Fiducial Efficiency;Time (years);E_{Ion} (keVee);Efficiency',
                                        self.fTimeBinning.size-1, self.fTimeBinning.flatten('C'),
                                        self.fEnergyIonBinning.size-1, self.fEnergyIonBinning.flatten('C')
                                        )

        self.fTotalEfficiency = TH3F(self.fName+'_total_efficiency',
                                     self.fName+' Total Efficiency; Time (years);E_{Rec} (keVnr);E_{Ion} (keVee);Efficiency',
                                     self.fTimeBinning.size-1, self.fTimeBinning.flatten('C'),
                                     self.fEnergyRecBinning.size-1, self.fEnergyRecBinning.flatten('C'),
                                     self.fEnergyIonBinning.size-1, self.fEnergyIonBinning.flatten('C')
                                     )

        self.fProjectedTriggerEfficiency = TH1F(self.fName+'_projected_trigger_efficiency',
                                                self.fName+' Projected Trigger Efficiency;E_{Rec} (keVnr);Efficiency',
                                                self.fEnergyRecBinning.size-1, self.fEnergyRecBinning.flatten('C')
                                                )

        self.fProjectedFiducialEfficiency = TH1F(self.fName+'_projected_fiducial_efficiency',
                                                 self.fName+' Projected Fiducial Efficiency;E_{Ion} (keVee);Efficiency',
                                                 self.fEnergyIonBinning.size-1, self.fEnergyIonBinning.flatten('C')
                                                 )

        self.fEnergyEfficiency = TH2F(self.fName+'_energy_efficiency',
                                      self.fName+' Energy Efficiency;E_{Rec} (keVnr);E_{Ion} (keVee);Efficiency',
                                      self.fEnergyRecBinning.size-1, self.fEnergyRecBinning.flatten('C'),
                                      self.fEnergyIonBinning.size-1, self.fEnergyIonBinning.flatten('C')
                                      )


        # Calculations to fill all efficiency histograms
        self._FillHistograms(self.fGoodUnixPeriods)
        self._CalcTriggerEfficiency()
        self._CalcFiducialEfficiency()
        self._CalcTotalEfficiency()
        self._CalcEnergyEfficiency()


    def _AssignMass(self):
        mass = parameters.DETECTOR_MASSES[self.fName]['Mass']
        return mass['Fiducial']


    def _ReadInGoodPeriods(self):
        inpath = 'Run12PeriodInformation/'

        UnixStart, UnixEnd, Run, Heat_Baseline, Heat_Threshold, Fiducial_Top, Fiducial_Bottom, Veto_Top, Veto_Bottom, Guard_Top, Guard_Bottom = [], [], [], [], [], [], [], [], [], [], []

        infile = open(inpath+self.fName+'_Bckgd.txt', 'r')

        ELEC = parameters.BASELINE_CUTS_ERIC[self.fName]
        for line in infile:
            unixstart, unixend, run, periodflag, heat_baseline, heat_threshold, fiducial_top, fiducial_bottom, veto_top, veto_bottom, guard_top, guard_bottom = line.split()

            startday = functions.unixtime_to_day(int(unixstart))
            endday = functions.unixtime_to_day(int(unixend))
            if functions.voltage_flag(self.fVoltageFlagList, str(run)) >= 0 \
               and startday not in self.fBadDays \
               and endday not in self.fBadDays \
               and ELEC['Coll1']['Min'] < float(fiducial_top) < ELEC['Coll1']['Max'] \
               and ELEC['Coll2']['Min'] < float(fiducial_bottom) < ELEC['Coll2']['Max'] \
               and ELEC['Veto1']['Min'] < float(veto_top) < ELEC['Veto1']['Max'] \
               and ELEC['Veto2']['Min'] < float(veto_bottom) < ELEC['Veto2']['Max'] \
               and ELEC['Guard1']['Min'] < float(guard_top) < ELEC['Guard1']['Max'] \
               and ELEC['Guard2']['Min'] < float(guard_bottom) < ELEC['Guard2']['Max'] \
               and ELEC['Heat']['Min'] < float(heat_baseline) < ELEC['Heat']['Max']:
                UnixStart.append(int(unixstart))
                UnixEnd.append(int(unixend))
                Run.append(str(run))
                Heat_Baseline.append(float(heat_baseline))
                Heat_Threshold.append(float(heat_threshold))
                Fiducial_Top.append(float(fiducial_top))
                Fiducial_Bottom.append(float(fiducial_bottom))
                Veto_Top.append(float(veto_top))
                Veto_Bottom.append(float(veto_bottom))
                Guard_Top.append(float(guard_top))
                Guard_Bottom.append(float(guard_bottom))
        infile.close()

        return [UnixStart, UnixEnd, Run, Heat_Baseline, Heat_Threshold, Fiducial_Top,
                Fiducial_Bottom, Veto_Top, Veto_Bottom, Guard_Top, Guard_Bottom]


    # Put gaps in between periods if not adjacent
    def _InsertGaps(self, InList):
        OutList = [ [ InList[i][0] ] for i in range(len(InList))] #copy first row of values to new list
        for i in range(1, len(InList[0])):
            unixstart = InList[0][i]
            unixstop = InList[1][i]
            lastunixstop = InList[1][i-1]
            if unixstart != lastunixstop:  # if gap between two periods
                OutList[0].append(lastunixstop)  # insert last unixtime end as start
                OutList[1].append(unixstart)  # insert next unixtime start as end
                for j in range(2, len(InList)):
                    OutList[j].append(0.)  # insert '0' for other values
            OutList[0].append(unixstart)
            OutList[1].append(unixstop)
            for j in range(2, len(InList)):
                OutList[j].append( InList[j][i] )
        return OutList


    # Calculate time binning in units of years
    def _CalcYearBins(self, InList):
        UnixTimeStart = InList[0]
        OutList = []
        for i in range (len(UnixTimeStart)):
            unixtime = UnixTimeStart[i]
            year = functions.unixtime_to_year(unixtime)
            OutList.append(year)
        return np.array(OutList, dtype=np.float)


    def _FillHistograms(self, InList):
        for i in range(len(InList[0])):
            unixtimestart = InList[0][i]
            unixtimeend = InList[1][i]
            run = InList[2][i]
            heat_baseline = InList[3][i]
            heat_threshold = InList[4][i]
            fiducial_top = InList[5][i]
            fiducial_bottom = InList[6][i]
            veto_top = InList[7][i]
            veto_bottom = InList[8][i]
            guard_top = InList[9][i]
            guard_bottom = InList[10][i]

            voltageflag = functions.voltage_flag(self.fVoltageFlagList, run)
            voltage = self.fERAConstants['gVolts'][voltageflag]

            time = functions.unixtime_to_year((unixtimeend + unixtimestart) / 2.)

            self.fHeatBaseline.Fill(time, heat_baseline)
            self.fHeatThreshold.Fill(time, heat_threshold)
            self.fFiducialTopBaseline.Fill(time, fiducial_top)
            self.fFiducialBottomBaseline.Fill(time, fiducial_bottom)
            self.fVetoTopBaseline.Fill(time, veto_top)
            self.fVetoBottomBaseline.Fill(time, veto_bottom)
            self.fGuardTopBaseline.Fill(time, guard_top)
            self.fGuardBottomBaseline.Fill(time, guard_bottom)

            fiducial_mean = (fiducial_top * fiducial_bottom) / ROOT.sqrt(fiducial_top**2 + fiducial_bottom**2)
            self.fFiducialMeanBaseline.Fill(time, fiducial_mean)

            self.fVoltage.Fill(time, voltage)
            self.fLivetimeEfficiency.Fill(time, 1.0)


    # Calculate TH2 histogram with trigger efficiency in (x,y)=(time,rec)
    def _CalcTriggerEfficiency(self):
        trigger_eff_hist = self.fTriggerEfficiency

        # Loop over time bins
        for xbin in range(1, trigger_eff_hist.GetNbinsX()+1):
            FWHM_heat = self.fHeatBaseline.GetBinContent(xbin)
            threshhold_heat = self.fHeatThreshold.GetBinContent(xbin)
            voltage = self.fVoltage.GetBinContent(xbin)

            sigma_heat = FWHM_heat/2.35

            threshhold_rec = functions.recoil_energy_estimator(threshhold_heat, voltage)
            sigma_rec = functions.fwhm_rec_from_heat(sigma_heat, voltage, threshhold_rec) #resolution at threshold energy

            EfficiencyCurve = functions.TRIGGER_EFFICIENCY
            EfficiencyCurve.FixParameter(0, threshhold_rec)
            EfficiencyCurve.FixParameter(1, sigma_rec)

            # Loop over recoil energy bins
            for ybin in range(1, trigger_eff_hist.GetNbinsY()+1):
                mean_energy = trigger_eff_hist.GetYaxis().GetBinCenter(ybin)
                value = EfficiencyCurve.Eval(mean_energy)

                if value >= 0.0 and FWHM_heat != 0.0 and threshhold_heat != 0.0:
                    efficiency = value
                else:
                    efficiency = 0.0

                trigger_eff_hist.SetBinContent(xbin, ybin, efficiency)
                trigger_eff_hist.SetBinError(xbin, ybin, 0.0)
        return None


    # Calculate TH2 histogram with fiducial efficiency in (x,y)=(time,ion)
    def _CalcFiducialEfficiency(self):
        fiducial_eff_hist = self.fFiducialEfficiency

        for xbin in range(1, fiducial_eff_hist.GetNbinsX()+1):  # time bin loop
            fiducial_baseline = self.fFiducialMeanBaseline.GetBinContent(xbin)

            EfficiencyCurve = functions.FIDUCIAL_EFFICIENCY

            max_eff = parameters.FIDUCIAL_EFFICIENCY_PARAMETERS[self.fName]['max']
            slope = parameters.FIDUCIAL_EFFICIENCY_PARAMETERS[self.fName]['slope']
            cutoff = parameters.FIDUCIAL_EFFICIENCY_PARAMETERS[self.fName]['cutoff']

            EfficiencyCurve.FixParameter(0, max_eff)
            EfficiencyCurve.FixParameter(1, slope)
            EfficiencyCurve.FixParameter(2, cutoff)

            for ybin in range(1, fiducial_eff_hist.GetNbinsY()+1):  # energy bin loop
                mean_energy = fiducial_eff_hist.GetYaxis().GetBinCenter(ybin)

                value = EfficiencyCurve.Eval(mean_energy)
                if value >= 0.0 and fiducial_baseline != 0.0:
                    efficiency = value
                else:
                    efficiency = 0.0

                fiducial_eff_hist.SetBinContent(xbin, ybin, efficiency)
                fiducial_eff_hist.SetBinError(xbin, ybin, 0.0)
        return None


    # Calculate the total energy efficiency TH3 histgram in (x,y,z)=(time,rec,ion)
    def _CalcTotalEfficiency(self):
        total_eff_hist = self.fTotalEfficiency
        livetime_eff_hist = self.fLivetimeEfficiency
        trigger_eff_hist = self.fTriggerEfficiency
        fiducial_eff_hist = self.fFiducialEfficiency

        # Loop over time bins
        for xbin in range(1, total_eff_hist.GetNbinsX()+1):
            livetime_eff = livetime_eff_hist.GetBinContent(xbin)

            # Loop over recoil energy bins
            for ybin in range(1, total_eff_hist.GetNbinsY()+1):
                trigger_eff = trigger_eff_hist.GetBinContent(xbin, ybin)

                # Loop over ionization energy bins
                for zbin in range(1, total_eff_hist.GetNbinsZ()+1):
                    fiducial_eff = fiducial_eff_hist.GetBinContent(xbin, zbin)

                    if livetime_eff != 0. and trigger_eff != 0. and fiducial_eff != 0.:
                        total_eff = livetime_eff * trigger_eff * fiducial_eff
                    else:
                        total_eff = 0.0

                    total_eff_hist.SetBinContent(xbin, ybin, zbin, total_eff)
                    total_eff_hist.SetBinError(xbin, ybin, zbin, 0.0)
        return None


    # Calculate the livetime weighted energy efficiency TH2 histogram in (x,y)=(rec,ion)
    def _CalcProjectedEfficiency(self, efficiency):
        if efficiency == 'trigger':
            inhist = self.fTriggerEfficiency
            outhist = self.fProjectedTriggerEfficiency
        elif efficiency == 'fiducial':
            inhist = self.fFiducialEfficiency
            outhist = self.fProjectedFiducialEfficiency

        # Loop over energy bins
        for ybin in range(1, inhist.GetNbinsY()+1):
            temp = 0.0

            # Loop over time bins
            for xbin in range(1, inhist.GetNbinsX()+1):
                timewidth = inhist.GetXaxis().GetBinWidth(xbin)
                efficiency = inhist.GetBinContent(xbin, ybin)
                temp += timewidth * efficiency

            temp /= self.GetLivetime()
            outhist.SetBinContent(ybin, temp)
            outhist.SetBinError(ybin, 0.0)
        return None


    def WriteEventsToFile(self, OutfileName):
        outfile = open(OutfileName, 'w')
        for i in range(len(self.fEventList[0])):
            UnixTime = self.fEventList[0][i]
            YearTime = unixtime_to_year(UnixTime)
            EnergyRec = self.fEventList[1][i]
            EnergyIon = self.fEventList[2][i]
            outfile.write("%s %s %s \n" % (YearTime, EnergyRec, EnergyIon))
        outfile.close()
        print 'Outfile %s written to disc' % OutfileName
        return True


    # Check if event fullfills energy binning and is in valid livetime period
    def IsGoodEvent(self, UnixTime, EnergyRec, EnergyIon):
        if binning['rec']['min'] <= EnergyRec <= binning['rec']['max'] and binning['ion']['min'] <= EnergyIon <= binning['ion']['max']:
            for entry in range(len(self.fGoodUnixPeriods[0])):  # loop over the list with good periods
                if UnixTime >= self.fGoodUnixPeriods[0][entry] and UnixTime <= self.fGoodUnixPeriods[1][entry]:
                    return True
                return false
        else:
            return False


    # Add time and energy information of event to list
    def AddEvent(self, Time, EnergyRec, EnergyIon):
        self.fEventList[0].append(Time)
        self.fEventList[1].append(EnergyRec)
        self.fEventList[2].append(EnergyIon)


    # Returns the name of the detector
    def GetName(self):
        return self.fName


    # Returns the mass of the detector in kg
    def GetMass(self):
        return self.fMass


    # Returns the live time of the detector in years
    def GetLivetime(self):
        return self.fLivetimeEfficiency.Integral('WIDTH')


    # Returns the time span of the detector data in years
    def GetRuntime(self):
        StartTime = self.fLivetimeEfficiency.GetXaxis().GetBinLowEdge(1)
        EndTime = self.fLivetimeEfficiency.GetXaxis().GetBinLowEdge(self.fLivetimeEfficiency.GetNbinsX()+1)
        return EndTime-StartTime


    # Returns detector exposure after cuts in kg.days
    def GetExposure(self):
        return self.GetLivetime()*365.*self.fMass


    # Returns TH1F histogram of livetime efficiency 0 = off, 1 = on with variable binning
    def GetLivetimeEfficiency(self):
        return self.fLivetimeEfficiency


    # Returns TH2F histogram of fiducial cut efficiency in (x,y)=(time,ion)
    def GetFiducialEfficiency(self):
        return self.fFiducialEfficiency


    # Returns TH2F histogram of trigger efficiency in (x,y)=(time,rec)
    def GetTriggerEfficiency(self):
        return self.fTriggerEfficiency


    # Returns TH3F histogram of total efficiency in (x,y,z)=(time,rec,ion)
    def GetTotalEfficiency(self):
        return self.fTotalEfficiency


    def GetNuclearRecoilEnergy(self, energy_ee, VoltageFlag):
        voltage = self.fERAConstants['gVolts'][VoltageFlag]
        energy = functions.recoil_energy_estimator(energy_ee, voltage)
        return energy


    # Return list of list with all baseline values
    def GetAllBaselines(self):
        outlist = [
            self.fHeatBaseline,
            self.fHeatThreshold,
            self.fFiducialTopBaseline,
            self.fFiducialBottomBaseline,
            self.fVetoTopBaseline,
            self.fVetoBottomBaseline,
            self.fGuardTopBaseline,
            self.fGuardBottomBaseline
            ]
        return outlist


    # Return variable time binning in years
    def GetTimeBinning(self):
        return self.fTimeBinning


    # Return efficiency in recoil or ionization energy projected over livetime
    def GetProjectedEfficiency(self, name):
        if name == 'trigger':
            if self.fProjectedTriggerEfficiency.GetEntries() == 0.0:
                self._CalcProjectedEfficiency('trigger')
            return self.fProjectedTriggerEfficiency
        elif name == 'fiducial':
            if self.fProjectedFiducialEfficiency.GetEntries() == 0.0:
                self._CalcProjectedEfficiency('fiducial')
            return self.fProjectedFiducialEfficiency


    #def GetWeightedAverage(self, name):
        #if name == 'voltage': hist = self.fVoltage
        #elif name == 'heat': hist = self.fHeatBaseline
        #elif name == 'threshold': hist = self.fHeatThreshold
        #elif name == 'fiducialtop': hist = self.fFiducialTopBaseline
        #elif name == 'fiducialbottom': hist = self.fFiducialBottomBaseline
        #elif name == 'fiducialmean': hist = self.fFiducialMeanBaseline
        #elif name == 'vetotop': hist = self.fVetoTopBaseline
        #elif name == 'vetobottom': hist = self.fVetoBottomBaseline
        #elif name == 'guardtop': hist = self.fGuardTopBaseline
        #elif name == 'guardbottom': hist = self.fGuardBottomBaseline

        #Temp = 0.0
        #for xbin in range(1, hist.GetNbinsX()+1):
            #timewidth = hist.GetXaxis().GetBinWidth(xbin)
            #value = hist.GetBinContent(xbin)
            #Temp += timewidth * value
        #Temp /= self.GetLivetime()
        #return Temp


    # Return livetime averaged value for baselines or voltage
    def GetWeightedAverage(self, name):
        if name == 'voltage': hist = self.fVoltage
        elif name == 'heat': hist = self.fHeatBaseline
        elif name == 'threshold': hist = self.fHeatThreshold
        elif name == 'fiducialtop': hist = self.fFiducialTopBaseline
        elif name == 'fiducialbottom': hist = self.fFiducialBottomBaseline
        elif name == 'fiducialmean': hist = self.fFiducialMeanBaseline
        elif name == 'vetotop': hist = self.fVetoTopBaseline
        elif name == 'vetobottom': hist = self.fVetoBottomBaseline
        elif name == 'guardtop': hist = self.fGuardTopBaseline
        elif name == 'guardbottom': hist = self.fGuardBottomBaseline

        avg_value = hist.Integral('WIDTH') / self.GetLivetime()

        return avg_value


    def GetEventGraphEnergy(self):
        events = len(self.fEventList[0])
        time = np.array(self.fEventList[0], dtype=np.float)
        E_rec = np.array(self.fEventList[1], dtype=np.float)
        E_ion = np.array(self.fEventList[2], dtype=np.float)
        graph = TGraph(events, E_rec.flatten('C'), E_ion.flatten('C'))
        graph.SetTitle('Events for '+self.fName)
        return graph


    # Calculates TH2F histogram of time weighted energy efficiency in (x,y)=(rec,ion)
    def _CalcEnergyEfficiency(self):
        inhist = self.fTotalEfficiency
        outhist = self.fEnergyEfficiency

        # Loop over ion bins
        for zbin in range(1, inhist.GetNbinsZ()+1):
            # Loop over rec bins
            for ybin in range(1, inhist.GetNbinsY()+1):
                Temp = 0.
                # Loop over time bins
                for xbin in range(1, inhist.GetNbinsX()+1):
                    timewidth = inhist.GetXaxis().GetBinWidth(xbin)
                    total_efficiency = inhist.GetBinContent(xbin, ybin, zbin)
                    Temp += timewidth * total_efficiency
                Temp /= self.GetLivetime()
                outhist.SetBinContent(ybin, zbin, Temp)
                outhist.SetBinError(ybin, zbin, 0.0)
        return None


    # Returns TH2F histogram of time weighted energy efficiency in (x,y)=(rec,ion)
    def GetEnergyEfficiency(self):
        return self.fEnergyEfficiency
