# MESMER

**`MESMER`** (**M**uon **E**lectron **S**cattering with **M**ultiple **E**lectromagnetic **R**adiation) is a Monte Carlo event generator for
high-precision simulation of muon-electron scattering at low enegies, developed for the MUonE experiment (see [here](https://web.infn.it/MUonE/) and/or [here](https://twiki.cern.ch/twiki/bin/view/MUonE/WebHome)).

### Authors
**`MESMER`** is developed at INFN, Sezione di Pavia, and Department of Physics, Università di Pavia (Italy).  
Authors are listed in the [`AUTHORS`](AUTHORS.md) file.

### Citing the code
We'd be grateful if you could cite the following papers when using the **`MESMER`** generator:  
1. [Carloni Calame *et al.*, Towards muon-electron scattering at NNLO, JHEP 11 (2020) 028](https://inspirehep.net/literature/1805205)  
2. [Alacevich *et al.*, Muon-electron scattering at NLO, JHEP 02 (2019) 155](https://inspirehep.net/literature/1703989)

## Prerequisites & Compilation
The program is mainly written in `Fortran 77` and it has been extensively tested with `GCC` compilers on `GNU/Linux`.

Three external libraries are used, [`LoopTools`](http://www.feynarts.de/looptools/), [`Collier`](https://collier.hepforge.org/) and
the `C` implementation of [`RANLUX`](https://luscher.web.cern.ch/luscher/ranlux/), shipped with the code under GPL-like licences. To be compiled,
the `Collier` library requires the availability of `cmake`.

Furthermore, the program includes three different parameterization for the hadronic vacuum polarization: one by [Fred Jegerlehner](http://www-com.physik.hu-berlin.de/~fjeger/software.html), one by [Keshavarzi-Nomura-Teubner](#footnotes) and one by [Fedor Ignatov](http://cmd.inp.nsk.su/~ignatov/vpl/).

The interface to Cern ROOT event format requires the [`ROOT`](https://root.cern/) framework to be installed.

A `Makefile` is provided and the executable `mesmer` is build by simply issuing the command `make`.

**NOTICE: as of now, the released code runs only at LO, NLO and NNLO. If NNLO extra leptonic pair production are selected, the execution stops with a warning.** 

## Running the code
Once compiled, **`MESMER`** is run by issuing in the working directory the command `./mesmer`, which displays a command prompt
where parameters for the run can be set.  
Alternatively, an [input data card (`input-example`)](input-example) is provided and it can be fed as input
by piping `./mesmer < input-example`.  
The order in which the input parameters and values are inserted in the data card or at the prompt is indifferent.
If a parameter/value is missing, defaults are used. The only mandatory rule is that the last input must be `run`.

The **`MESMER`** prompt	 looks like:

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

By typing `help` at the prompt, a short description of the parameters that can be set is displayed (more details [here](#running-parameters-description)):

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
             └> [1] root file (see README.md)
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
 sync      ---> syncronization mode for random numbers (see README.md) [0/1]
```
## Running parameters description

The parameters that can be set are split into *principal* and *internal* parameters.  
They are set by typing at the prompt `parameter value` (case sensitive).

### Principal parameters

* `Qmu [default 1]`: charge of the incoming muon in e+ units `[1 / -1]`
* `Ebeam [150]`: nominal incoming muon energy in GeV
* `bspr [0]`: Gaussian beam energy spread, in % of `Ebeam`
* `Eemin [1]`: minimum outgoing electron energy in GeV
* `themax [100]`: maximum outgoing electron angle in mrad
* `thmumin [0]`: minimum outgoing muon angle in mrad
* `thmumax [100]`: maximum outgoing muon angle in mrad
* `acoplcut [no 3.5]`: if the acoplanarity cut has to be applied or not `[yes/no mrad]`
* `elastcut [no 0.2]`: if the cut on the geometric distance from the elesticity curve in the [&theta;<sub>&mu;</sub>,&theta;<sub>e</sub>] plane has to be applied or not `[yes/no mrad]`
* `ord [alpha]`: to which order (photonic) corrections are included `[born/alpha/alpha2]` (standing for `[LO/NLO/NNLO]`)
* `arun [on]`: if vacuum polarization (VP) effects must be includer or not. Possible values are
  * `off`: no VP
  * `on`: both leptonic and hadronic VP are included. The hadronic part defaults to `hadr5`
  * `hadr5`: include leptonic VP and for hadronic VP use latest Fred Jegerlehner's parameterization
  * `nsk`: include leptonic VP and for hadronic VP use Fedor Ignatov's parameterization
  * `knt`: include leptonic VP and for hadronic VP use Keshavarzi-Nomura-Teubner's parameterization
* `hadoff [no]`: if hadronic VP has to be switched off `[yes/no]`
* `nev [10000000.]`: number of events to be generated
* `store [no]`: if events have to be stored `[yes/no]`
* `storemode [0]`: if `store yes`, which mode to use to store events `[0/1/2]` (see [below](#description-of-event-format) for a short description of the event format)
  * `0`: plain ASCII text file
  * `1`: `ROOT` format. **`MESMER`** runs in parallel `write-root-events` (a symbolic link to `root-interface/write_MuE_MCevents.exe`), developed by G. Abbiendi, which writes through a named pipe a `.root` file with the events
  * `2`: an `xz` compressed file is saved
* `path [test-run/]`: the directory where all outputs are saved. It will contain some `.txt` files with differential distributions of some variables, the file `events-*.[dat,dat.xz,root]` with saved events if `store yes` and the file `stat_*.txt`, where cross sections and all info of the current run are reported
* `seed [42]`: seed for the pseudo-random-number-generator, it must be set to a "small" integer. Independent generations must use different seeds


### Inner parameters

* `mode [weighted]`: if stored events are weighted or unweighted `[weighted/unweighted]`.  
If `weighted`, to each event is associated a different weight, which must be appropriately accounted for when producing distributions.  
If `unweighted`, all events have the same constant weight, *i.e.* they are directly distributed according to the underlying cross sections.  
**Notice:** unweighted generation can be much slower, due to the unweightening procedure
* `radchs [1 1]`: "radiative" charges of the muon and electron leg, to switch on/off gauge invariants subsets `[0/1 0/1]`
* `eps [1e-5]`: minimum photon CM energy of the emitted photons in sqrt(s)/2 units. Physical observables must not (and do not) depend on `eps`
* `phmass [1e-10]`: photon mass to regularize infra-red divergecies in GeV. Physical observables do not depend on `phmass`
* `nphot [-1]`: it's a code driving which final states are included and the multiplicity of particles in the final state.
Possible values are
  * `< 0 and > -1000`: up to one extra photon in the final state for `ord alpha`, *i.e.* NLO, and up to two photons for `ord alpha2`, *i.e.* NNLO. No extra photons for `ord born`, *i.e.* LO
  * `0 or 1 or 2`: exactly this number of extra photons in the final state
  * `1000`: simulates
  *&mu;<sup>&plusmn;</sup>e<sup>-</sup>&rarr;&mu;<sup>&plusmn;</sup>e<sup>-</sup>e<sup>+</sup>e<sup>-</sup>*
  and
  *&mu;<sup>&plusmn;</sup>e<sup>-</sup>&rarr;&mu;<sup>&plusmn;</sup>e<sup>-</sup>&mu;<sup>+</sup>&mu;<sup>-</sup>*
  together
  * `1001`: simulates *&mu;<sup>&plusmn;</sup>e<sup>-</sup>&rarr;&mu;<sup>&plusmn;</sup>e<sup>-</sup>e<sup>+</sup>e<sup>-</sup>*
  * `1002`: simulates *&mu;<sup>&plusmn;</sup>e<sup>-</sup>&rarr;&mu;<sup>&plusmn;</sup>e<sup>-</sup>&mu;<sup>+</sup>&mu;<sup>-</sup>*
  * `-1000`: all possible final states, *i.e.* maximum number of photons plus extra leptonic pairs
* `nwrite [500000]`: files in `path` are dumped every `nwrite` generated events. If `nwrite < 0`, files are dumped (approximately) every `-nwrite` seconds
* `nwarmup [0]`:  if `mode weighted`, the maximum weight for subsequent unweightening is searched generating first `nwarmup` events, after which unweighted generation is started. Notice that the maximum weight is a guessed value.  
If `store yes` and `nwarmup > 0`, `wnorm` is calculated automatically using the first `nwarmup` events of the generation.  
For unweightening, the maximum weight `sdmax` can be alternatively set by hand (setting `nwarmup 0` at the same time)
* `sdmax [1e-17]`: maximum weight used for unweightening, when `mode unweighted` and `nwarmup 0`
* `wnorm [-1.]`: typical integrated cross section within applied cuts. This value is used for storing events: the true weight of the event *w<sub>t</sub>* is saved as *w = w<sub>t</sub> / wnorm* in order to make the average of the weights *w* over the generated sample close to 1
* `ndistr [1]`: number of distributions at different orders. For example, if running at NNLO (`ord alpha2`), the distribution files saved in `path` can be produced also at NLO and LO with the same run. The defaults `ndistr 1` produces distributions only at the selected order. If `ord alpha` and `ndistr >=2`, distributions are produced at NLO and LO. If `ord alpha2`,  if `ndistr 2` they are produced at NNLO and NLO, if `ndistr 3`, distributions at NNLO, NLO and LO are saved
* `sync [0]`: syncronization mode for random numbers `[0/1]`. This is used to make as correlated as possible two samples with same `seed` but different `Ebeam`. *This feature needs more cross-checks*

## Description of event format
The event format has been developed together with G. Abbiendi.

Version 1 is compliant with the current version of the `ROOT` file writer ([mantained by G. Abbiendi](https://gitlab.cern.ch/muesli/nlo-mc/mue/-/tree/master/writer)), but it is limited to parse only LO, NLO and NNLO events (*events with extra leptonic pairs are not supported yet*).  

**A more general event format will be agreed upon soon.**

### Version 1
The simple ASCII file (`storemode 0`) contains an header section between the tags `<header>` and `</header>` and a closing section between the tags `<footer>` and `</footer>`. After the header, each *event record* is enclosed between the `<event>` and `</event>` tags.  
The header contains useful info about the run and the set parameters. The footer contains some statistics of the generated sample.  
A typical *event record* looks like
```
<event>
          42
                 9998
           3
   6.1762125544134905E-002   6.0976085936283446E-002   6.1756917995296243E-002
   0.0000000000000000     
   150.00000000000000     
   145.44875212920633        5.1471445841701699E-002  -4.2645249941300672E-002   145.44869839336738     
   2.2013510818264361       -2.4148349195201722E-002   1.7863188831557644E-002   2.2011460851234839     
   2.3504077876745284       -2.7323096646499977E-002   2.4782061109743032E-002   2.3501183089789106     
 </event>
```
After the `<event>` tag, the first line is the `seed` of the sample, the second one is the event number and the third is the number of final state particles (in this case a &mu;, an *e* and a &gamma;).  
In the fourth line, three weights *w* are listed: the weight with full VP effect, the one with VP switched off and the one with only leptonic VP effects (without hadronic VP): in this way, with a single run VP effects can be studied.  
The fifth line represents the weight *w<sub>LO</sub>* which can be used to get LO distributions from a NLO sample (in this case it is `0` because it's a 3-body final state).  
The sixth line is the incoming &mu; energy in GeV.  
Finally, the last three lines represent the momenta (in GeV) of the final state &mu;, *e* and &gamma; respectively, in the format E p<sub>x</sub> p<sub>y</sub> p<sub>z</sub>.  
The tag `</event>` closes the *event record*.

The `ROOT` file writer reads through a named pipe this event format and writes all the information into a `.root` file.

### Version 2
Not released yet. To be agreed upon

##### Footnotes
KNT: Available upon request from the authors