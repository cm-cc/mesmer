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
The program is mainly written in `Fortran 77` and it has been extensively tested with `GCC` compilers on `GNU/Linux`.

Three external libraries are used, [`LoopTools`](http://www.feynarts.de/looptools/), [`Collier`](https://collier.hepforge.org/) and
the `C` implementation of [`RANLUX`](https://luscher.web.cern.ch/luscher/ranlux/), shipped with the code under GPL-like licences. To be compiled,
the `Collier` library requires the availability of the `cmake`.

The interface to Cern ROOT event format requires the [`ROOT`](https://root.cern/) framework to be installed.

A `Makefile` is provided and the executable `mesmer` is build by simply issuing the command `make`.

## Running the code
Once compiled, **`MESMER`** is run by issuing in the working directory the command `./mesmer`, which displays a command prompt
where parameters for the run can be set.  
Alternatively, an [input data card (`input-example`)](input-example) is provided and it can be fed as input
by piping `./mesmer < input-example`.  
The order in which the input parameters and values are inserted in the data card or at the prompt is indifferent.
If a parameter/value is missing, defaults are used. The only rule is that the last input must be `run`.

The **`MESMER`** prompy looks like:

```
  *************************************************************
  ********                                             ********
  ******               Welcome to MESMER                 ******
  ****              ~~~~~~~~~~~~~~~~~~~~~~                 ****
  **       A fully exclusive Monte Carlo event generator     **
  ****          for Muon Electron Scattering  with         ****
  ******        Multiple Electromagnetic Radiation       ******
  ********                                             ********
  *************************************************************
  
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
   [ elastcut  ] Apply "elasticity" cut (and value) = no (0.2 mrad)
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

By typing `help` at the prompt, a short description of the parameters that can be set is displayed:

```
 Principal parameters:
 Qmu       ---> incoming muon charge in e+ charge units [-1/1]
 Ebeam     ---> nominal muon beam energy (GeV)
 bspr      ---> beam energy spread percentage (%)
 Eemin     ---> minimum electron energy in the LAB (GeV)
 themax    ---> maximum electron angle in the LAB (mrad)
 thmumin   ---> minimum muon angle in the LAB (mrad)
 thmumax   ---> maximum muon angle in the LAB (mrad)
 acoplcut  ---> if applying acoplanarity cut and to which value [yes/no mrad]
 elastcut  ---> if applying distance-from-elasticity-curve cut and to which value [yes/no mrad]
 ord       ---> to which order simulate events [born/alpha/alpha2] for [LO/NLO/NNLO]
 arun      ---> if running of alpha must be used [off/on/hadr5/nsk/knt]
 hadoff    ---> if switching off hadronic VP [yes/no]
 nev       ---> number of events to generate
 store     ---> if events have to be stored [yes/no]
 storemode ---> which mode to use to store events [0/1/2]
             └> [0] plain ascii file 
             └> [1] root file (see README)
             └> [2] on-the-fly xz compressed ascii file
 path      ---> path where to store outputs
 seed      ---> pseudo-RNG seed ("small" int)
  
 Internal tweaks:
 mode      ---> if stored events are weighted or unweighted [weighted/unweighted]
 radchs    ---> "radiative" charges of the muon and electron leg, to switch on/off gauge invariants subsets [0/1 0/1]
 eps       ---> minimum photon CM energy of the emitted photons in sqrt(s)/2 units. A.k.a. soft/hard separator
 phmass    ---> photon mass (IR regulator) (GeV)
 nphot     ---> drives maximum number of hard photons (> eps) and real pair emission [int]
             └> [ < 0 ] maximum possible (up to 2)
             └> [ 0 | 1 | 2 ] exactly this number of photons
             └> [ 1000 ] both real electron and muon pairs
             └> [ 1001 ] only real electron pairs
             └> [ 1002 ] only real muon pairs
             └> [-1000 ] all possible final states
 nwrite    ---> files are dumped every nwrite generated events [int] (negative integer means files written every -nwrite seconds)
 nwarmup   ---> after nwarmup events, also unweighted generation is started [int]. Plays with wnorm.
 ndistr    ---> if writing distributions also at different orders [1/2/3]
 sdmax     ---> maximum integrand for unweightening
 wnorm     ---> typical integrated cross section within applied cuts, used for storage in ROOT format
 sync      ---> syncronization mode for random numbers (see README) [0/1]
```