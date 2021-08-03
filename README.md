# MESMER

***THIS README.md FILE IS UNDER CONSTRUCTION.***

**`MESMER`** (Muon Electron Scattering with Multiple Electromagnetic Radiation) is a Monte Carlo event generator for high-precision simulation of muon-electron scattering at low enegies

### Authors
**`MESMER`** is developed at INFN, Sezione di Pavia, and Department of Physics, Università di Pavia (Italy).  
Authors are listed in the file [`AUTHORS`](AUTHORS.md).

### Citing the code
We'd be grateful if you could cite the following papers when using the **`MESMER`** generator:  
1. [Carloni Calame *et al.*, Towards muon-electron scattering at NNLO, JHEP 11 (2020) 028](https://inspirehep.net/literature/1805205)  
2. [Alacevich *et al.*, Muon-electron scattering at NLO, JHEP 02 (2019) 155](https://inspirehep.net/literature/1703989)

## Prerequisites & Compilation
The program is mainly written in `Fortran 77` and it has been tested with `GCC` compilers.

Three external libraries are used, [`LoopTools`](http://www.feynarts.de/looptools/), [`Collier`](https://collier.hepforge.org/) and
the `C` implementation of [`RANLUX`](https://luscher.web.cern.ch/luscher/ranlux/), shipped with the code under GPL-like licences. To be compiled,
the `Collier` library requires the availability of the `cmake`.

The interface to Cern ROOT event format requires the [`ROOT`](https://root.cern/) framework to be installed.

A `Makefile` is provided and the executable `mesmer` is build by simply issuing the command `make`.

## Running the code


```bash
./mesmer 
  *************************************************************
  ********                                             ********
  ******               Welcome to MESMER                 ******
  ****              ~~~~~~~~~~~~~~~~~~~~~~                 ****
  **       A fully exclusive Monte Carlo event generator     **
  ****          for Muon Electron Scattering  with         ****
  ******        Multiple Electromagnetic Radiation       ******
  ********                                             ********
  *************************************************************
    Questo è quel pezzo / di calamita:
    pietra mesmerica, / ch'ebbe l'origine
    nell'Alemagna, / Che poi sì celebre / là in Francia fu.
           Mozart/Da Ponte, Così Fan Tutte, Act II, Scene XVI
  
 Principal parameters:
   [ type "run" to start generation, "help" for help or "quit" to abort ]
   [ Qmu       ] Incoming muon charge = 1 (in e+ charge units)
   [ Ebeam     ] Muon beam energy = 150. GeV
   [ bspr      ] Beam energy spread = 0. %
   [ Eemin     ] Minimum electron LAB energy = 1. GeV
   [ themax    ] Maximum electron LAB angle = 100. mrad
   [ thmumin   ] Minimum muon LAB angle = 0. mrad
   [ thmumax   ] Maximum muon LAB angle = 100. mrad
   [ acoplcut  ] Apply acoplanarity cut (and value) = no (3.5 mrad)
   [ elastcut  ] Apply "elasticity" cut (and value) = no (0.200000000000000011 mrad)
   [ ord       ] Simulation at "alpha" order
   [ arun      ] alpha running is on
   [ hadoff    ] Hadronic VP is off: no
   [ nev       ] n. of events to generate: 10000000.
   [ store     ] events storage: no
   [ storemode ] mode to store: 0
   [ path      ] files saved in test-run/
   [ seed      ] seed for pseudo-RNG = 42
  
 Internal tweaks:
   [ mode      ] requested evts. are weighted
   [ radchs    ] "radiative" charges: muon ~> 1 electron ~> 1
   [ eps       ] soft photon cutoff in sqrt(s)/2 units = 0.00002
   [ phmass    ] photon mass = 0.0000000001 GeV
   [ nphot     ] max. number of photons mode is -1
   [ nwrite    ] file(s) dumped every 500000 events
   [ nwarmup   ] events for maximum searching = 0
   [ ndistr    ] number of distr. at different orders: 1
   [ sdmax     ] starting "sdifmax" 0.00000000000000001
   [ wnorm     ] normalization cross section = -1.
   [ sync      ] random numbers sequence syncronization: 0 [0/1]
  
 Insert "parameter value" or "run" or "quit": 
```