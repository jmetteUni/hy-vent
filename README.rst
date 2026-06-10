Hy-Vent
=======

A package wich aims to provide functions for analyzing CTD Data from hydrothermal vents. Is the basis for the [analysis scripts](https://github.com/jmetteUni/aurora-plume-dispersal-figures) used in the master thesis [Mette, 2025](https://media.suub.uni-bremen.de/handle/elib/22668).

Data is in general handled as pandas dataframes, where data of one station is one dataframe. Several stations (that means dataframes) of the same type are stored in a dicitonary with th key being a unique station identifier containing the cruisename, the station number and the cast number.
Example: "PS137_018_01" means cruise Polarstern 137, station 018 and cast 01. Station and cast number are exspected to have 3 or 2 digits respectively for now.

Preliminary documentation can be found [here](https://hy-vent.readthedocs.io/en/latest/).
