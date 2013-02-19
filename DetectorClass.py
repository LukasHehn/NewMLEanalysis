#!/usr/bin/env python
from ROOT import *
from Parameters import *
from Functions import *
import numpy as np


class Detector:
  def __init__(self, name):
    self.fName = name

    self.fMass = self._AssignMass()

    self.fBadDays = []

    self.fEventList = [[], [], []] #time & energy list of added events

    self.fRunType = 'Bckgd'

    self.fVoltageFlagList = GetVoltageFlagList(self.fName, self.fRunType)

    self.fERAConstants = GetERAConstants(self.fName)

    self.fGoodUnixPeriods = self._ReadInGoodPeriods()

    self.fGoodUnixPeriodsGap = self._InsertGaps(self.fGoodUnixPeriods) #fit gaps between periods

    self.fTimeBinning = self._CalcYearBins(self.fGoodUnixPeriodsGap) #get numpy array with time bins in years

    self.fEnergyRecBinning = np.array([Energy['rec']['min']+i*Energy['rec']['binsize'] for i in range(0,Energy['rec']['bins']+1)], dtype=np.float)

    self.fEnergyIonBinning = np.array([Energy['ion']['min']+i*Energy['ion']['binsize'] for i in range(0,Energy['ion']['bins']+1)], dtype=np.float)

    # histograms
    self.fHeatBaseline = TH1F(self.fName+'_heat_baseline', self.fName+' heat baseline;Time (years);Energy (keV)', \
                                  self.fTimeBinning.size-1, self.fTimeBinning.flatten('C'))
    self.fHeatThreshold = TH1F(self.fName+'_heat_threshold', self.fName+' heat threshold;Time (years);Energy (keV)', \
                                   self.fTimeBinning.size-1, self.fTimeBinning.flatten('C'))
    self.fFiducialTopBaseline = TH1F(self.fName+'_fiducial_top_baseline', self.fName+' fiducial top baseline;Time (years);Energy (keV)', \
                                         self.fTimeBinning.size-1, self.fTimeBinning.flatten('C'))
    self.fFiducialBottomBaseline = TH1F(self.fName+'_fiducial_bottom_baseline', self.fName+' fiducial bottom baseline;Time (years);Energy (keV)', \
                                            self.fTimeBinning.size-1, self.fTimeBinning.flatten('C'))
    self.fVetoTopBaseline = TH1F(self.fName+'_veto_top_baseline', self.fName+' veto top baseline;Time (years);Energy (keV)', \
                                     self.fTimeBinning.size-1, self.fTimeBinning.flatten('C'))
    self.fVetoBottomBaseline = TH1F(self.fName+'_veto_bottom_baseline', self.fName+' veto bottom baseline;Time (years);Energy (keV)', \
                                        self.fTimeBinning.size-1, self.fTimeBinning.flatten('C'))
    self.fGuardTopBaseline = TH1F(self.fName+'_guard_top_baseline', self.fName+' guard top baseline;Time (years);Energy (keV)', \
                                      self.fTimeBinning.size-1, self.fTimeBinning.flatten('C'))
    self.fGuardBottomBaseline = TH1F(self.fName+'_guard_bottom_baseline', self.fName+' guard bottom baseline;Time (years);Energy (keV)', \
                                         self.fTimeBinning.size-1, self.fTimeBinning.flatten('C'))
    self.fFiducialMeanBaseline = TH1F(self.fName+'_fiducial_mean', self.fName+' fiducial mean baseline;Time (years);Energy (keV)', \
                                         self.fTimeBinning.size-1, self.fTimeBinning.flatten('C'))

    self.fVoltage = TH1F(self.fName+'_voltage', self.fName+' voltage;Time (years);Voltage (V)', \
                             self.fTimeBinning.size-1, self.fTimeBinning.flatten('C'))
    self.fLivetimeEfficiency = TH1F(self.fName+'_livetime_efficiency', self.fName+' Livetime Efficiency;Time (years) ;Efficiency', \
                                        self.fTimeBinning.size-1, self.fTimeBinning.flatten('C'))

    self.fTriggerEfficiency = TH2F(self.fName+'_trigger_efficiency', self.fName+' Trigger Efficiency;Time (years); E_{Rec} (keV_{nr}); Efficiency', \
                                   self.fTimeBinning.size-1, self.fTimeBinning.flatten('C'), self.fEnergyRecBinning.size-1, self.fEnergyRecBinning.flatten('C'))

    self.fFiducialEfficiency = TH2F(self.fName+'_fiducial_efficiency', self.fName+' Fiducial Efficiency;Time (years); E_{Ion} (keV_{ee}); Efficiency', \
                                    self.fTimeBinning.size-1, self.fTimeBinning.flatten('C'), self.fEnergyIonBinning.size-1, self.fEnergyIonBinning.flatten('C'))

    self.fTotalEfficiency = TH3F(self.fName+'_total_efficiency', self.fName+' Total Efficiency;Time (years); E_{Rec} (keV_{nr}); E_{Ion} (keV_{ee}); Efficiency', \
                                     self.fTimeBinning.size-1, self.fTimeBinning.flatten('C'), self.fEnergyRecBinning.size-1, self.fEnergyRecBinning.flatten('C'), self.fEnergyIonBinning.size-1, self.fEnergyIonBinning.flatten('C'))


    self.fProjectedTriggerEfficiency = TH1F(self.fName+'_projected_trigger_efficiency', self.fName+' Projected Trigger Efficiency;E_{Rec} (keV_{nr});Efficiency', \
                                      self.fEnergyRecBinning.size-1, self.fEnergyRecBinning.flatten('C'))
    self.fProjectedFiducialEfficiency = TH1F(self.fName+'_projected_fiducial_efficiency', self.fName+' Projected Fiducial Efficiency;E_{Ion} (keV_{ee});Efficiency', \
                                      self.fEnergyIonBinning.size-1, self.fEnergyIonBinning.flatten('C'))
    self.fProjectedTotalEfficiency = TH2F(self.fName+'_projected_total_efficiency', self.fName+' Projected Total Efficiency;E_{Rec} (keV_{nr});E_{Ion} (keV_{ee});Efficiency', \
                                      self.fEnergyRecBinning.size-1, self.fEnergyRecBinning.flatten('C'), self.fEnergyIonBinning.size-1, self.fEnergyIonBinning.flatten('C'))


    self._FillHistograms(self.fGoodUnixPeriods)
    self._CalcTriggerEfficiency()
    self._CalcFiducialEfficiency()
    self._CalcTotalEfficiency()
    self._CalcProjectedEfficiency('trigger')
    self._CalcProjectedEfficiency('fiducial')
    self._CalcProjectedEfficiency('total')



  def _AssignMass(self):
    mass = Masses[self.fName]['Mass']
    return mass['Fiducial']


  def _ReadInGoodPeriods(self):
    inpath = 'Run12PeriodInformation/'

    UnixStart, UnixEnd, Run, Heat_Baseline, Heat_Threshold, Fiducial_Top, Fiducial_Bottom, Veto_Top, Veto_Bottom, Guard_Top, Guard_Bottom = [], [], [], [], [], [], [], [], [], [], []

    infile = open(inpath+self.fName+'_Bckgd.txt','r')

    ELEC = EricsLowEnergyCuts[self.fName]
    for line in infile:
      unixstart, unixend, run, periodflag, heat_baseline, heat_threshold, fiducial_top, fiducial_bottom, veto_top, veto_bottom, guard_top, guard_bottom = line.split()

      startday = unixtime_to_day(int(unixstart))
      endday = unixtime_to_day(int(unixend))
      if GetVoltageFlag(self.fVoltageFlagList, str(run)) >= 0 \
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

    return [UnixStart, UnixEnd, Run, Heat_Baseline, Heat_Threshold, Fiducial_Top, Fiducial_Bottom, Veto_Top, Veto_Bottom, Guard_Top, Guard_Bottom]


  def _InsertGaps(self, InList): #put gaps in between periods if not adjacent
    OutList = [ [ InList[i][0] ]  for i in range(len(InList))] #copy first row of values to new list
    for i in range(1,len(InList[0])):
      unixstart = InList[0][i]
      unixstop = InList[1][i]
      lastunixstop = InList[1][i-1]
      if unixstart != lastunixstop: #if gap between two periods
        OutList[0].append(lastunixstop) #insert last unixtime end as start
        OutList[1].append(unixstart) #insert next unixtime start as end
        for j in range(2,len(InList)): OutList[j].append(0.) #insert '0' for other values
      OutList[0].append(unixstart)
      OutList[1].append(unixstop)
      for j in range(2,len(InList)): OutList[j].append( InList[j][i] )
    return OutList


  def _CalcYearBins(self, InList):
    UnixTimeStart = InList[0]
    OutList = []
    for i in range (len(UnixTimeStart)):
      unixtime = UnixTimeStart[i]
      year = unixtime_to_year(unixtime)
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

      voltageflag = GetVoltageFlag(self.fVoltageFlagList, run)
      voltage = self.fERAConstants['gVolts'][voltageflag]

      time = unixtime_to_year((unixtimeend+unixtimestart)/2)

      self.fHeatBaseline.Fill(time, heat_baseline)
      self.fHeatThreshold.Fill(time, heat_threshold)
      self.fFiducialTopBaseline.Fill(time, fiducial_top)
      self.fFiducialBottomBaseline.Fill(time, fiducial_bottom)
      self.fVetoTopBaseline.Fill(time, veto_top)
      self.fVetoBottomBaseline.Fill(time, veto_bottom)
      self.fGuardTopBaseline.Fill(time, guard_top)
      self.fGuardBottomBaseline.Fill(time, guard_bottom)

      fiducial_mean = (fiducial_top * fiducial_bottom)/sqrt(pow(fiducial_top, 2) + pow(fiducial_bottom, 2))
      self.fFiducialMeanBaseline.Fill(time, fiducial_mean)

      self.fVoltage.Fill(time, voltage)
      self.fLivetimeEfficiency.Fill(time, 1)


  def _CalcTriggerEfficiency(self):
    trigger_eff_hist = self.fTriggerEfficiency

    for xbin in range(1,trigger_eff_hist.GetNbinsX()+1):
      heat_baseline_ee = self.fHeatBaseline.GetBinContent(xbin)
      heat_threshold_ee = self.fHeatThreshold.GetBinContent(xbin)
      voltage = self.fVoltage.GetBinContent(xbin)

      heat_baseline_nr = GetEnergyRecoilFromEstimator(heat_baseline_ee, voltage)
      heat_threshold_nr = GetEnergyRecoilFromEstimator(heat_threshold_ee, voltage)

      heat_sigma_nr = heat_baseline_nr / (2*sqrt(2*log(2)))

      EfficiencyCurve = TriggerEfficiency
      EfficiencyCurve.SetParameter(0, heat_threshold_nr)
      EfficiencyCurve.SetParameter(1, heat_sigma_nr)

      for ybin in range(1, trigger_eff_hist.GetNbinsY()+1):
	mean_energy = trigger_eff_hist.GetYaxis().GetBinCenter(ybin)
	value = EfficiencyCurve.Eval(mean_energy)

	if value >= 0 and heat_baseline_ee != 0 and heat_threshold_ee != 0:
	  efficiency = value
	else:
	  efficiency = 0

	trigger_eff_hist.SetBinContent(xbin, ybin, efficiency)
	trigger_eff_hist.SetBinError(xbin, ybin, 0)
    return None


  def _CalcFiducialEfficiency(self):
    fiducial_eff_hist = self.fFiducialEfficiency

    for xbin in range(1,fiducial_eff_hist.GetNbinsX()+1): #time bin loop
      fiducial_baseline = self.fFiducialMeanBaseline.GetBinContent(xbin)

      EfficiencyCurve = FiducialEfficiency
      EfficiencyCurve.SetParameter(0, -1.87)
      #EfficiencyCurve.SetParameter(1, 1.25)
      EfficiencyCurve.SetParameter(1, 2*fiducial_baseline)

      for ybin in range(1, fiducial_eff_hist.GetNbinsY()+1): #energy bin loop
	mean_energy = fiducial_eff_hist.GetYaxis().GetBinCenter(ybin)

	value = EfficiencyCurve.Eval(mean_energy)
	if value >= 0 and fiducial_baseline != 0:
	  efficiency = value
	else:
	  efficiency = 0.0

	fiducial_eff_hist.SetBinContent(xbin, ybin, efficiency)
	fiducial_eff_hist.SetBinError(xbin, ybin, 0.0)
    return None


  def _CalcTotalEfficiency(self):
    total_eff_hist = self.fTotalEfficiency
    livetime_eff_hist = self.fLivetimeEfficiency
    trigger_eff_hist = self.fTriggerEfficiency
    fiducial_eff_hist = self.fFiducialEfficiency

    for xbin in range(1, total_eff_hist.GetNbinsX()+1):
      livetime_eff = livetime_eff_hist.GetBinContent(xbin)

      for ybin in range(1, total_eff_hist.GetNbinsY()+1):
	trigger_eff = trigger_eff_hist.GetBinContent(xbin, ybin)

	for zbin in range(1, total_eff_hist.GetNbinsZ()+1):
	  fiducial_eff = fiducial_eff_hist.GetBinContent(xbin, zbin)

	  if livetime_eff != 0 and trigger_eff != 0 and fiducial_eff != 0:
	    total_eff = livetime_eff * trigger_eff * fiducial_eff
	  else:
	    total_eff = 0

	  total_eff_hist.SetBinContent(xbin, ybin, zbin, total_eff)
	  total_eff_hist.SetBinError(xbin, ybin, zbin, 0)
    return None


  def _CalcProjectedEfficiency(self, efficiency): #calculate livetime-weighted average energy efficiency hist
    if efficiency in ['trigger', 'fiducial']:
      if efficiency == 'trigger':
	inhist = self.fTriggerEfficiency
	outhist = self.fProjectedTriggerEfficiency
      elif efficiency == 'fiducial':
	inhist = self.fFiducialEfficiency
	outhist = self.fProjectedFiducialEfficiency

      for ybin in range(1, inhist.GetNbinsY()+1):
	Temp = 0
	for xbin in range(1,inhist.GetNbinsX()+1):
	  timewidth = inhist.GetXaxis().GetBinWidth(xbin)
	  efficiency = inhist.GetBinContent(xbin, ybin)
	  Temp += timewidth * efficiency
	outhist.SetBinContent(ybin, Temp)
	outhist.SetBinError(ybin, 0)
    elif efficiency == 'total':
      inhist = self.fTotalEfficiency
      outhist = self.fProjectedTotalEfficiency

      for zbin in range(1, inhist.GetNbinsZ()+1):
        for ybin in range(1, inhist.GetNbinsY()+1):
          Temp = 0
          for xbin in range(1,inhist.GetNbinsX()+1):
            timewidth = inhist.GetXaxis().GetBinWidth(xbin)
            efficiency = inhist.GetBinContent(xbin, ybin, zbin)
            Temp += timewidth * efficiency
          outhist.SetBinContent(ybin, zbin, Temp)
          outhist.SetBinError(ybin, zbin, 0)
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


  def IsGoodEvent(self, UnixTime, EnergyRec, EnergyIon):
    if Energy['rec']['min'] <= EnergyRec <= Energy['rec']['max'] and Energy['ion']['min'] <= EnergyIon <= Energy['ion']['max']:
      for entry in range(len(self.fGoodUnixPeriods[0])): #loop over the list with good periods
	if UnixTime >= self.fGoodUnixPeriods[0][entry] and UnixTime <= self.fGoodUnixPeriods[1][entry]:
	  return True
      return false
    else:
      return False


  def AddEvent(self, Time, EnergyRec, EnergyIon):
    self.fEventList[0].append(Time)
    self.fEventList[1].append(EnergyRec)
    self.fEventList[2].append(EnergyIon)


  def GetName(self):
    return self.fName


  def GetMass(self):
    return self.fMass


  def GetLivetime(self):
    return self.fLivetimeEfficiency.Integral('WIDTH')


  def GetRuntime(self):
    StartTime = self.fLivetimeEfficiency.GetXaxis().GetBinLowEdge(1)
    EndTime = self.fLivetimeEfficiency.GetXaxis().GetBinLowEdge(self.fLivetimeEfficiency.GetNbinsX()+1)
    return EndTime-StartTime


  def GetExposure(self): # in kg*days
    return self.GetLivetime()*365*self.fMass


  def GetLivetimeEfficiency(self):
    return self.fLivetimeEfficiency


  def GetFiducialEfficiency(self):
    return self.fFiducialEfficiency


  def GetTriggerEfficiency(self):
    return self.fTriggerEfficiency


  def GetTotalEfficiency(self):
    return self.fTotalEfficiency


  def GetNuclearRecoilEnergy(self, energy_ee, VoltageFlag):
    voltage = self.fERAConstants['gVolts'][VoltageFlag]
    energy = GetEnergyRecoilFromEstimator(energy_ee, voltage)
    return energy


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


  def GetTimeBinning(self):
    return self.fTimeBinning


  def GetProjectedEfficiency(self, name):
    if name == 'trigger':
      return self.fProjectedTriggerEfficiency
    elif name == 'fiducial':
      return self.fProjectedFiducialEfficiency
    elif name == 'total':
      return self.fProjectedTotalEfficiency


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

    TotalTime = self.fLivetimeEfficiency.Integral('WIDTH')
    Weighted_Value_Sum = 0
    for xbin in range(1,hist.GetNbinsX()+1):
      time = hist.GetXaxis().GetBinWidth(xbin)
      value = hist.GetBinContent(xbin)
      weighted_value = time * value
      Weighted_Value_Sum += weighted_value
    average_value = Weighted_Value_Sum / TotalTime
    return average_value

  def GetEventGraphEnergy(self):
    events = len(self.fEventList[0])
    time = np.array(self.fEventList[0], dtype=np.float)
    E_rec = np.array(self.fEventList[1], dtype=np.float)
    E_ion = np.array(self.fEventList[2], dtype=np.float)
    graph = TGraph(events,E_rec.flatten('C'),E_ion.flatten('C'))
    graph.SetTitle('Events for '+self.fName)
    return graph
