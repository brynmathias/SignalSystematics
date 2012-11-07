SignalSystemacits
=================

Scripts to look at systematic studies on signal scans.

Mostly hacky thrown together stuff to make a few plots in root, these include:

BTag scale uncertianties - a ± one sigma variation around the central value

Jet energy scale uncertianties - again a ± one sigma variation around the central jet energy scale value ~ 1% per jet, more for low PT jets

Cleaning cuts - Take out each of MHT/MET, DeadEcal, Lepton Vetos and plot efficiency of these cuts.

EffMap - Makes simple maps of the efficiencys, greatly improved by sam. Mostly un-needed as this is done in:
http://github.com/elaird/ra1stats along with all the stats.


Dependancies: http://github.com/brynmathias/pyRootUtils for the Print wrapper and pulling TH2's from root files.
