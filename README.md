# MESMER

***THIS README.md FILE IS UNDER CONSTRUCTION.***

**MESMER** (Muon Electron Scattering with Multiple Electromagnetic Radiation) is a Monte Carlo event generator for high-precision simulation of muon-electron scattering at low enegies

### Authors
**MESMER** is developed at INFN, Sezione di Pavia, and Department of Physics, Università di Pavia, by

Carlo M. Carloni Calame      (carlo.carloni.calame[THIS IS A AT]pv.infn.it)  
Ettore Budassi               (ettore.budassi01[THIS IS A AT]universitadipavia.it)  
Mauro Chiesa                 (mauro.chiesa[THIS IS A AT]unipv.it)
Clara L. Del Pio             (claralavinia.delpio01[THIS IS A AT]universitadipavia.it)  
Syed M. Hasan                (syedmehe[THIS IS A AT]pv.infn.it)  
Guido Montagna               (guido.montagna[THIS IS A AT]pv.infn.it)  
Oreste Nicrosini             (oreste.nicrosini[THIS IS A AT]pv.infn.it)  
Fulvio Piccinini             (fulvio.piccinini[THIS IS A AT]pv.infn.it)  

## Prerequisites & Compilation
The program is mainly written in Fortran 77 and it has been tested with GCC compilers.

Three external libraries are used, [`LoopTools`](http://www.feynarts.de/looptools/), [`Collier`](https://collier.hepforge.org/) and
the `C` implementation of [RANLUX](https://luscher.web.cern.ch/luscher/ranlux/), shipped with the code under GPLs licences. The `Collier` library
requires the availability of the `cmake` command to be compiled.

The interface to Cern ROOT event format requires the [ROOT](https://root.cern/) framework to be installed.

A `Makefile` is provided and the executable `mesmer` is build by simply issuing the command `make`.

