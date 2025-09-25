# EoS_T0 module

This module provides several cold equations of state (EoS) for sgrid.

Supported are:

* Piecewise polytropes (choose: EoS_type = PwP).

* Tabulated cold 1d EoS (choose: EoS_type = tab1d_AtT0).
  See also [EoS_tab1d/Note.txt](EoS_tab1d/Note.txt)

* The script [EoS_tab1d/ConvertEOSTab.py](EoS_tab1d/ConvertEOSTab.py)
  can convert an EoS file from Lorene format to Sgrid's format used in
  this module.
