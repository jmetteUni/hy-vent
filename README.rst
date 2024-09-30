UNDER DEVELOPMENT Hy-Vent
=========================

A package wich aims to provide functions for analyzing CTD Data from hydrothermal vents.
Data is in general handled as pandas dataframes, where data of one station is one dataframe. Several stations (that means dataframes) of the same type are stored in a dicitonary with th key being a unique station identifier containing the cruisename, the station number and the cast number.
Example: "PS137_018_01" means cruise Polarstern 137, station 018 and cast 01. Station and cast number are exspected to have 3 or 2 digits respectively for now.


Contains the following modules:

io.py
-----
Input-output functions for reading and writing standardized CTD data and relevent metadata.

misc.py
-------
Miscellaneous functions.

processing.py
-------------
Functions for calculations based on the data.

quality_control.py
------------------
Functions for quality control.

plotting.py
-----------
Basic plot functions.

plotting_maps.py
----------------
Plot functions for 2D overview maps, right now only working for a specific dataset.

quality_control.py
------------------
Functions for quality control of position data from acoustic position trackers

tests
-----------------
Irrelevant test files.

develop.py
----------
Code, which is under heavy development. Working functions will be transferred to other modules.

Dependencies for this code are the following:
- numpy,
- pandas,
- datetime,
- seabird,
- matplotlib,
- lat_lon_parser,
- os
- gsw
