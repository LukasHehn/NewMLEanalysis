#!/usr/bin/env python
ENERGY_BINSIZE = 0.1
ENERGY_BINNING = {
  'ion' : {'Min' : 0., 'Max' : 10., 'Binsize' : ENERGY_BINSIZE},\
  'rec' : {'Min' : 0., 'Max' : 20., 'Binsize' : ENERGY_BINSIZE}
}
ENERGY_BINNING['ion']['Bins'] = int((ENERGY_BINNING['ion']['Max']-ENERGY_BINNING['ion']['Min'])/ENERGY_BINNING['ion']['Binsize'])
ENERGY_BINNING['rec']['Bins'] = int((ENERGY_BINNING['rec']['Max']-ENERGY_BINNING['rec']['Min'])/ENERGY_BINNING['rec']['Binsize'])


# measured fiducial and physical masses for all ID detectors from Run 12
DETECTOR_MASSES = {
   'ID2' : {'Mass' : {'Fiducial' : 0.158, 'Total': 0.3723}},\
   'ID3' : {'Mass' : {'Fiducial' : 0.158, 'Total': 0.3658}},\
   'ID4' : {'Mass' : {'Fiducial' : 0.158, 'Total': 0.3590}},\
   'ID5' : {'Mass' : {'Fiducial' : 0.158, 'Total': 0.3640}},\
   'ID6' : {'Mass' : {'Fiducial' : 0.148, 'Total': 0.3718}},\
   'ID401' : {'Mass' : {'Fiducial' : 0.168, 'Total': 0.410}},\
   'ID402' : {'Mass' : {'Fiducial' : 0.168, 'Total': 0.410}},\
   'ID403' : {'Mass' : {'Fiducial' : 0.168, 'Total': 0.410}},\
   'ID404' : {'Mass' : {'Fiducial' : 0.168, 'Total': 0.410}},\
   'ID405' : {'Mass' : {'Fiducial' : 0.168, 'Total': 0.410}}
   }


BASELINE_CUTS_ERIC = {
   'ID3' : {
      'Coll1' : {'Min' : 0.1, 'Max' : 1.0},\
      'Coll2' : {'Min' : 0.1, 'Max' : 1.8},\
      'Veto1' : {'Min' : 0.1, 'Max' : 1.4},\
      'Veto2' : {'Min' : 0.1, 'Max' : 1.2},\
      'Guard1' : {'Min' : 0.1, 'Max' : 2.5},\
      'Guard2' : {'Min' : 0.1, 'Max' : 2.0},\
      'Heat' : {'Min' : 0.1, 'Max' : 1.0},\
     },\
   'ID5' : {
      'Coll1' : {'Min' : 0.1, 'Max' : 2.5},\
      'Coll2' : {'Min' : 0.1, 'Max' : 1.1},\
      'Veto1' : {'Min' : 0.1, 'Max' : 2.9},\
      'Veto2' : {'Min' : 0.1, 'Max' : 1.7},\
      'Guard1' : {'Min' : 0.1, 'Max' : 1.3},\
      'Guard2' : {'Min' : -999, 'Max' : 999},  # no data from Eric available\
      'Heat' : {'Min' : 0.1, 'Max' : 1.4},\
     },\
   'ID6' : {
      'Coll1' : {'Min' : 0.1, 'Max' : 1.3},\
      'Coll2' : {'Min' : 0.1, 'Max' : 2.0},\
      'Veto1' : {'Min' : 0.1, 'Max' : 1.3},\
      'Veto2' : {'Min' : 0.1, 'Max' : 2.5},\
      'Guard1' : {'Min' : -999, 'Max' : 999},  # no data from Eric available\
      'Guard2' : {'Min' : 0.1, 'Max' : 1.2},\
      'Heat' : {'Min' : 0.1, 'Max' : 0.65},\
     },\
   'ID401' : {
      'Coll1' : {'Min' : 0.1, 'Max' : 1.1},\
      'Coll2' : {'Min' : 0.1, 'Max' : 1.7},\
      'Veto1' : {'Min' : 0.1, 'Max' : 1.6},\
      'Veto2' : {'Min' : 0.1, 'Max' : 1.6},\
      'Guard1' : {'Min' : 0.1, 'Max' : 7.0},\
      'Guard2' : {'Min' : 0.1, 'Max' : 1.6},\
      'Heat' : {'Min' : 0.1, 'Max' : 1.4},\
      },\
   'ID404' : {
      'Coll1' : {'Min' : 0.1, 'Max' : 1.7},\
      'Coll2' : {'Min' : 0.1, 'Max' : 2.3},\
      'Veto1' : {'Min' : 0.1, 'Max' : 2.0},\
      'Veto2' : {'Min' : 0.1, 'Max' : 2.0},\
      'Guard1' : {'Min' : 0.1, 'Max' : 2.2},\
      'Guard2' : {'Min' : 0.1, 'Max' : 1.6},\
      'Heat' : {'Min' : 0.1, 'Max' : 1.3},\
      }\
   }

LOW_MASS_RESULTS_ERIC = {
  'ID3' : {
    7 : 7.90E-4,\
    8 : 1.17e-4,\
    10 : 1.35e-5,\
    12 : 4.17e-6,\
    15 : 1.55e-6,\
    20 : 7.34e-7,\
    25 : 5.47e-7,\
    30 : 4.85e-7
    }
}