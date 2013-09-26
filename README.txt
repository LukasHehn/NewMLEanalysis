This repository contains all scripts necessary to perform an unbinned maximum likelihood analysis in recoil and ionization energy. NOT contained are the files with event information (data) and the run information (baselines, voltages, etc).

Here's a short description of the more important modules:

- mle.py: the main analysis script for maximum likelihood estimation
- detector.py: holds Detector-class that incorporates all detector specific information like calculation of live time periods, efficiency histograms etc
- functions.py: a collection of functions to convert dates, energies, ER and NR bands, etc.
- parameters.py: holds dictionaries with often used parameters like detector masses, baseline cuts, baseline values, etc.
- find_events.py: used to find low energy events with a selection of cuts, save them to a txt-file and plot them together with WIMP signal and gamma-background.

Necessary libraries are (list might be incomplete):
- ROOT
- RooFit
- pyWIMP
- calender
- datetime
- couchdbkit
