#!/usr/bin/env python

####################################################################################################
##
## Frequently used parameters like detector masses, low mass paper cross section limit results and
## energy binning choice.
## Lukas Hehn, 2013
##
####################################################################################################


# Choice of range for ionization and recoil energy as well as binsize and number of bins
ENERGY_BINSIZE = 1.0
ENERGY_BINNING = {
    'ion' : {'min' : 0., 'max' : 10., 'binsize' : ENERGY_BINSIZE},\
    'rec' : {'min' : 0., 'max' : 20., 'binsize' : ENERGY_BINSIZE}
}
ENERGY_BINNING['ion']['bins'] = int((ENERGY_BINNING['ion']['max']-ENERGY_BINNING['ion']['min'])/ENERGY_BINNING['ion']['binsize'])
ENERGY_BINNING['rec']['bins'] = int((ENERGY_BINNING['rec']['max']-ENERGY_BINNING['rec']['min'])/ENERGY_BINNING['rec']['binsize'])


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


# Upper limits on SI WIMP-nucleon scattering cross sections taken from technical/internal
# complements to the low mass EDW2 paper: "Search for low-energy nuclear recoils with ID data" from
# E. Armengaud, April 2012 (updated 5 July 2012)
LOW_MASS_RESULTS_ERIC = {
    'ID2' : {
        7 : 1.6E-1,
        8 : 4.20E-3,
        10 : 1.20E-4,
        12 : 2.14E-5,
        15 : 5.47E-6,
        20 : 2.02E-6,
        25 : 1.36E-6,
        30 : 1.15E-6
    },
    'ID3' : {
        7 : 7.90E-4,
        8 : 1.17e-4,
        10 : 1.35e-5,
        12 : 4.17e-6,
        15 : 1.55e-6,
        20 : 7.34e-7,
        25 : 5.47e-7,
        30 : 4.85e-7
    },
    'ID6' : {
        7 : 3.2E-3,
        8 : 2.84E-4,
        10 : 3.64E-5,
        12 : 1.28E-5,
        15 : 5.17E-6,
        20 : 2.20E-6,
        25 : 1.57E-6,
        30 : 1.37E-6
    },
    'ID401' : {
        7 : 3.9E-3,
        8 : 3.85E-4,
        10 : 3.01E-5,
        12 : 7.86E-6,
        15 : 2.61E-6,
        20 : 1.14E-6,
        25 : 8.26E-7,
        30 : 7.25E-7
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
        'Heat' : 0.58, # 0.524
        'Coll1' : 1.02,
        'Coll2' : 1.59,
        'Veto1' : 0.95,
        'Veto2' : 1.76,
        'Guard1' : 999,  # dead
        'Guard2' : 0.92,
        'Fiducial' : 0.86  # 0.834
        },
    'ID401' : {
        'Heat' : 1.31,  # me: 1.05keV baseline
        'Coll1' : 0.96,
        'Coll2' : 1.49,
        'Veto1' : 1.40,
        'Veto2' : 1.55,
        'Guard1' : 4.19,
        'Guard2' : 1.33,
        'Fiducial' : 0.81  # me: 0.76keV baseline
        },
    'ID404' : {
        'Heat' : 1.08,  # 0.956
        'Coll1' : 1.35,
        'Coll2' : 1.90,
        'Veto1' : 1.37,
        'Veto2' : 1.57,
        'Guard1' : 1.84,
        'Guard2' : 1.25,
        'Fiducial' : 1.10  # 1.048
        }
}


# Average voltages found in Eric's internal noted (presumably after baseline cuts)
# Online Threshold in keVee crudely read from internal notes plots
MEASURED_VALUES_ERIC = {
    'ID2' : {'livetime' : 120, 'voltage' : 0.0, 'threshold_ee' : 1.7},
    'ID3' : {'livetime' : 197, 'voltage' : 6.4, 'threshold_ee' : 1.8},
    'ID6' : {'livetime' : 247, 'voltage' : 0.0, 'threshold_ee' : 1.1},
    'ID401' : {'livetime' : 145, 'voltage' : 7.6, 'threshold_ee' : 2.2},
}


# Live time weighted average voltage values
# Online Threshold in keVee from live time weighted values (detector class)
MEASURED_VALUES_LUKAS = {
    'ID3' : {'livetime' : 197.6, 'voltage' : 6.40000009537, 'threshold_ee' : 1.98964249642, 'threshold_nr' : 4.27758077779},
    'ID6' : {'livetime' : 147.1, 'voltage' : 8.0, 'threshold_ee' : 1.09144984615, 'threshold_nr' : 2.64166124022},
    'ID401' : {'livetime' : 145.4, 'voltage' : 7.55912944983, 'threshold_ee' : 2.52587734932, 'threshold_nr' : 5.64077829638},
}


# Parameters for the fitted fiducial efficiency function
FIDUCIAL_EFFICIENCY_PARAMETERS = {
    'ID2' : {'max' : 0.99, 'slope' : -1.73, 'cutoff' : 1.74},
    'ID3' : {'max' : 0.947, 'slope' : -1.876, 'cutoff' : 1.247},  # ID3 has lowest (=best) cutoff
    'ID6' : {'max' : 0.96, 'slope' : -1.23, 'cutoff' : 1.43},
    'ID401' : {'max' : 0.97, 'slope' : -2.386, 'cutoff' : 1.394},
    'ID404' : {'max' : 0.97, 'slope' : -6.34, 'cutoff' : 1.90}
}
