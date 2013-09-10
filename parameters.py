#!/usr/bin/env python

####################################################################################################
#
# Frequently used parameters like detector masses, low mass paper cross section limit results and
# energy binning choice.
# Lukas Hehn, 2013
#
####################################################################################################


# Choice of range for ionization and recoil energy as well as binsize and number of bins
ENERGY_BINSIZE = 0.1
ENERGY_BINNING = {
    'ion' : {'Min' : 0., 'Max' : 10., 'Binsize' : ENERGY_BINSIZE},\
    'rec' : {'Min' : 0., 'Max' : 20., 'Binsize' : ENERGY_BINSIZE}
}
ENERGY_BINNING['ion']['Bins'] = int((ENERGY_BINNING['ion']['Max']-ENERGY_BINNING['ion']['Min'])/ENERGY_BINNING['ion']['Binsize'])
ENERGY_BINNING['rec']['Bins'] = int((ENERGY_BINNING['rec']['Max']-ENERGY_BINNING['rec']['Min'])/ENERGY_BINNING['rec']['Binsize'])


# Measured fiducial and physical masses for all ID detectors used in Run 12
DETECTOR_MASSES = {
    'ID2' : {'Mass' : {'Fiducial' : 0.158, 'Total': 0.3723}},
    'ID3' : {'Mass' : {'Fiducial' : 0.158, 'Total': 0.3658}},
    'ID4' : {'Mass' : {'Fiducial' : 0.158, 'Total': 0.3590}},
    'ID5' : {'Mass' : {'Fiducial' : 0.158, 'Total': 0.3640}},
    'ID6' : {'Mass' : {'Fiducial' : 0.148, 'Total': 0.3718}},
    'ID401' : {'Mass' : {'Fiducial' : 0.168, 'Total': 0.410}},
    'ID402' : {'Mass' : {'Fiducial' : 0.168, 'Total': 0.410}},
    'ID403' : {'Mass' : {'Fiducial' : 0.168, 'Total': 0.410}},
    'ID404' : {'Mass' : {'Fiducial' : 0.168, 'Total': 0.410}},
    'ID405' : {'Mass' : {'Fiducial' : 0.168, 'Total': 0.410}}
}


# Selection of baseline cuts used in the published EDW low mass analysis
BASELINE_CUTS_ERIC = {
    'ID3' : {
        'Coll1' : {'Min' : 0.1, 'Max' : 1.0},
        'Coll2' : {'Min' : 0.1, 'Max' : 1.8},
        'Veto1' : {'Min' : 0.1, 'Max' : 1.4},
        'Veto2' : {'Min' : 0.1, 'Max' : 1.2},
        'Guard1' : {'Min' : 0.1, 'Max' : 2.5},
        'Guard2' : {'Min' : 0.1, 'Max' : 2.0},
        'Heat' : {'Min' : 0.1, 'Max' : 1.0},
        },
    'ID5' : {
        'Coll1' : {'Min' : 0.1, 'Max' : 2.5},
        'Coll2' : {'Min' : 0.1, 'Max' : 1.1},
        'Veto1' : {'Min' : 0.1, 'Max' : 2.9},
        'Veto2' : {'Min' : 0.1, 'Max' : 1.7},
        'Guard1' : {'Min' : 0.1, 'Max' : 1.3},
        'Guard2' : {'Min' : -999, 'Max' : 999},  # no data from Eric available\
        'Heat' : {'Min' : 0.1, 'Max' : 1.4}
        },
    'ID6' : {
        'Coll1' : {'Min' : 0.1, 'Max' : 1.3},
        'Coll2' : {'Min' : 0.1, 'Max' : 2.0},
        'Veto1' : {'Min' : 0.1, 'Max' : 1.3},
        'Veto2' : {'Min' : 0.1, 'Max' : 2.5},
        'Guard1' : {'Min' : -999, 'Max' : 999},  # no data from Eric available\
        'Guard2' : {'Min' : 0.1, 'Max' : 1.2},
        'Heat' : {'Min' : 0.1, 'Max' : 0.65}
        },
    'ID401' : {
        'Coll1' : {'Min' : 0.1, 'Max' : 1.1},
        'Coll2' : {'Min' : 0.1, 'Max' : 1.7},
        'Veto1' : {'Min' : 0.1, 'Max' : 1.6},
        'Veto2' : {'Min' : 0.1, 'Max' : 1.6},
        'Guard1' : {'Min' : 0.1, 'Max' : 7.0},
        'Guard2' : {'Min' : 0.1, 'Max' : 1.6},
        'Heat' : {'Min' : 0.1, 'Max' : 1.4}
        },
    'ID404' : {
        'Coll1' : {'Min' : 0.1, 'Max' : 1.7},
        'Coll2' : {'Min' : 0.1, 'Max' : 2.3},
        'Veto1' : {'Min' : 0.1, 'Max' : 2.0},
        'Veto2' : {'Min' : 0.1, 'Max' : 2.0},
        'Guard1' : {'Min' : 0.1, 'Max' : 2.2},
        'Guard2' : {'Min' : 0.1, 'Max' : 1.6},
        'Heat' : {'Min' : 0.1, 'Max' : 1.3}
        }
}


# Upper limits for SI WIMP-nucleon scattering cross sections from published low mass analysis
LOW_MASS_RESULTS_ERIC = {
    'ID3' : {
        7 : 7.90E-4,
        8 : 1.17e-4,
        10 : 1.35e-5,
        12 : 4.17e-6,
        15 : 1.55e-6,
        20 : 7.34e-7,
        25 : 5.47e-7,
        30 : 4.85e-7
    }
}


# measured heat and combined fiducial ionization resolutions (@10keV?) taken from technical/internal 
# complements to the low mass EDW2 paper: "Search for low-energy nuclear recoils with ID data" from
# E. Armengaud, April 2012 (updated 5 July 2012)
# probably after baseline cuts
ENERGY_RESOLUTIONS_ERIC = {
    'ID3' : {
        'Heat' : 0.82,
        'Coll1' : 0.80,
        'Coll2' : 1.60,
        'Veto1' : 1.19,
        'Veto2' : 1.02,
        'Guard1' : 1.90,
        'Guard2' : 1.66,
        'Fiducial' : 0.72
        },
    'ID5' : {
        'Heat' : 1.11,
        'Coll1' : 2.07,
        'Coll2' : 0.90,
        'Veto1' : 2.21,
        'Veto2' : 1.41,
        'Guard1' : 1.31,
        'Guard2' : 999,  # dead
        'Fiducial' : 0.83
        },
    'ID6' : {
        'Heat' : 0.58,
        'Coll1' : 1.02,
        'Coll2' : 1.59,
        'Veto1' : 0.95,
        'Veto2' : 1.76,
        'Guard1' : 0,  # dead
        'Guard2' : 0.92,
        'Fiducial' : 0.86
        },
    'ID401' : {
        'Heat' : 1.31,
        'Coll1' : 0.96,
        'Coll2' : 1.49,
        'Veto1' : 1.40,
        'Veto2' : 1.55,
        'Guard1' : 4.19,
        'Guard2' : 1.33,
        'Fiducial' : 0.81
        },
    'ID404' : {
        'Heat' : 1.08,
        'Coll1' : 1.35,
        'Coll2' : 1.90,
        'Veto1' : 1.37,
        'Veto2' : 1.57,
        'Guard1' : 1.84,
        'Guard2' : 1.25,
        'Fiducial' : 1.10
        }
}


# probably after baseline cuts
AVG_VOLTAGES_ERIC = {
    'ID3' : 6.4,
    'ID5' : 0.0,
    'ID6' : 0.0,
    'ID401' : 7.6,
    'ID404' : 7.2
}
