# MESMER
**MESMER** (Muon Electron Scattering with Multiple Electromagnetic Radiation) is a Monte Carlo event generator for high-precision simulation of muon-electron scattering at low enegies

***The README.md file is under construction.***

## Prerequisites & Compilation
The program is mainly written in Fortran 77 and it has been tested with GCC compilers.

Two external libraries are used, [LoopTools](http://www.feynarts.de/looptools/) and [Collier](https://collier.hepforge.org/), shipped with the code
under LGPL and GPLv3 licenses respectively. The latter requires the availability of the `cmake` command to be compiled.

The interface to Cern ROOT event format requires the [ROOT](https://root.cern/) suite to be installed.

A `Makefile` is provided and the executable `mesmer` is build by simply issuing the command `make`.

### Authors
**MESMER** is developed at INFN, Sezione di Pavia (Italy) by

Carlo M. Carloni Calame (carlo.carloni.calame[THIS IS A AT]pv.infn.it)  
Ettore Budassi  
Clara L. Del Pio  
Syed M. Hasan  
Guido Montagna  
Oreste Nicrosini  
Fulvio Piccinini  