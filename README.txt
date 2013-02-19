This repository contains my Maximum Likelihood analysis with all relevant files:

- DetectorClass.py (the main class that incorporates all detector specific information like calculation of live time periods, efficiency histograms etc)
- Functions.py (a collection of functions to convert dates, energies, etc)
- MLEanalysis.py (the main analysis script for maximum likelihood estimation)
- Parameters.py (some cut parameters)
- FindKDataEvents.py (used to find low energy events with a selection of cuts and save them to a txt-file)

Necessary libraries are:
- ROOT
- RooFit
- pyWIMP
- calender
- datetime
- couchdbkit
