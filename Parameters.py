#!/usr/bin/env python
binsize = 0.1
Energy = {
  'ion' : {'min' : 0, 'max' : 10, 'binsize' : binsize},\
  'rec' : {'min' : 0, 'max' : 20, 'binsize' : binsize}
}
Energy['ion']['bins'] = int((Energy['ion']['max']-Energy['ion']['min'])/Energy['ion']['binsize'])
Energy['rec']['bins'] = int((Energy['rec']['max']-Energy['rec']['min'])/Energy['rec']['binsize'])


# measured fiducial and physical masses for all ID detectors from Run 12
Masses = {
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


EricsLowEnergyCuts = {
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
      'Guard2' : {'Min' : -999, 'Max' : 999},#no data from Eric\
      'Heat' : {'Min' : 0.1, 'Max' : 1.4},\
     },\
   'ID6' : {
      'Coll1' : {'Min' : 0.1, 'Max' : 1.3},\
      'Coll2' : {'Min' : 0.1, 'Max' : 2.0},\
      'Veto1' : {'Min' : 0.1, 'Max' : 1.3},\
      'Veto2' : {'Min' : 0.1, 'Max' : 2.5},\
      'Guard1' : {'Min' : -999, 'Max' : 999},#no data from Eric\
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

EricsResults = {
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

Eric = [[7,8,10,12,15,20,25,30],[7.90E-4,1.17E-4,1.35E-5,4.17E-6,1.55E-6,7.34E-7,5.47E-7,4.85E-7]]