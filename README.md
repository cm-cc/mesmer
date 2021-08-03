# MESMER

***THIS README.md FILE IS UNDER CONSTRUCTION.***

**MESMER** (Muon Electron Scattering with Multiple Electromagnetic Radiation) is a Monte Carlo event generator for high-precision simulation of muon-electron scattering at low enegies

### Authors
**MESMER** is developed at INFN, Sezione di Pavia, and Department of Physics, Università di Pavia.

Authors are listed in the file [AUTHORS](AUTHORS).

## Prerequisites & Compilation
The program is mainly written in `Fortran 77` and it has been tested with `GCC` compilers.

Three external libraries are used, [`LoopTools`](http://www.feynarts.de/looptools/), [`Collier`](https://collier.hepforge.org/) and
the `C` implementation of [`RANLUX`](https://luscher.web.cern.ch/luscher/ranlux/), shipped with the code under GPLs licences. The `Collier` library
requires the availability of the `cmake` command to be compiled.

The interface to Cern ROOT event format requires the [`ROOT`](https://root.cern/) framework to be installed.

A `Makefile` is provided and the executable `mesmer` is build by simply issuing the command `make`.

