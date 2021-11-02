# MESMER

**`MESMER`** (**M**uon **E**lectron **S**cattering with **M**ultiple **E**lectromagnetic **R**adiation) is a Monte Carlo event generator for
high-precision simulation of muon-electron scattering at low enegies, developed for the MUonE experiment (see [here](https://web.infn.it/MUonE/) and/or [here](https://twiki.cern.ch/twiki/bin/view/MUonE/WebHome)).

### Authors
**`MESMER`** is developed at INFN, Sezione di Pavia, and Department of Physics, Università di Pavia (Italy).  
Authors are listed in the [`AUTHORS`](AUTHORS.md) file.

### Citing the code
We'd be grateful if you could cite the following papers when using the **`MESMER`** generator:  
1. [Alacevich *et al.*, Muon-electron scattering at NLO, JHEP 02 (2019) 155](https://inspirehep.net/literature/1703989)
2. [Carloni Calame *et al.*, Towards muon-electron scattering at NNLO, JHEP 11 (2020) 028](https://inspirehep.net/literature/1805205)  
3. [E. Budassi *et al.*, NNLO virtual and real leptonic corrections to muon-electron scattering, 2109.14606 [hep-ph]](https://inspirehep.net/literature/1933852)

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
   [ extmubeam ] Feed generic muon beam externally (beam profile): no
   [ Eemin     ] Minimum electron LAB energy = 1. GeV
   [ themin    ] Minimum electron LAB angle = 0. mrad
   [ themax    ] Maximum electron LAB angle = 100. mrad
   [ thmumin   ] Minimum muon LAB angle = 0. mrad
   [ thmumax   ] Maximum muon LAB angle = 100. mrad
   [ acoplcut  ] Apply acoplanarity cut (and value) = no (3.5 mrad)
   [ elastcut  ] Apply "elasticity" cut (and value) = no (0.2 mrad)
   [ Ethr      ] 'Detectability' energy threshold for extra leptons = 0.2 GeV
   [ ththr     ] Maximum 'detectability' LAB angle for extra leptons = 100. mrad
   [ ord       ] Simulation at "alpha" order
   [ arun      ] alpha running is on
   [ hadoff    ] Hadronic VP is off: no
   [ nev       ] n. of events to generate: 10000000.
   [ store     ] events storage: no old no
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
   [ sdmax     ] starting "sdifmax" 0.0000000001
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
 extmubeam ---> if muon beam 3-momentum is externally fed (beam profile) [yes/no].
                If 'yes', the beam 3-momentum must be fed into the pipe 'path'/beamprofile.fifo
 Eemin     ---> minimum electron energy in the LAB (GeV)
 themin    ---> minimum electron angle in the LAB (mrad)
 themax    ---> maximum electron angle in the LAB (mrad)
 thmumin   ---> minimum muon angle in the LAB (mrad)
 thmumax   ---> maximum muon angle in the LAB (mrad)
 acoplcut  ---> if applying acoplanarity cut and to which value [yes/no mrad]
 elastcut  ---> if applying distance-from-elasticity-curve cut and to which value [yes/no mrad]
 Ethr      ---> 'detectability' threshold energy for extra leptons (GeV)
 ththr     ---> maximum 'detectability' angle for extra leptons  (mrad)
 ord       ---> to which order simulate events [born/alpha/alpha2] for [LO/NLO/NNLO]
 arun      ---> if running of alpha must be used [off/on/hadr5/nsk/knt]
 hadoff    ---> if switching off hadronic VP [yes/no]
 nev       ---> number of events to generate
 store     ---> (3 parameters) if events have to be stored [yes/no],
                which format [new/old] (for v2/v1),
                if writing also coefficients for VP reweighting [yes/no] (active only for v2)
 storemode ---> which mode to use to store events [0/1/2]
             └> [0] plain ascii file 
             └> [1] root file (see README)
             └> [2] on-the-fly xz compressed ascii file
             └> [3] just write to the fifo and wait for an external process to read events from it	     
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
## Running parameters description

In general, the routine `cuts(...)` in the file `cuts.F` can be modified according to the needs of the user. The default one applies the selection criteria described below. The only mandatory generator level cut is `Eemin`, which **must** be stricty larger than *m<sub>e</sub>*. The other cuts can be set at their kinematical limits.

The parameters that can be set are split into *principal* and *internal* parameters and they are set by typing at the prompt `parameter value(s)` (case sensitive).

### Principal parameters

* `Qmu [default 1]`: charge of the incoming muon in e+ units `[1 / -1]`
* `Ebeam [150]`: nominal incoming muon energy in GeV
* `extmubeam [no]`: feed generic muon beam externally (beam profile, see [below](#beam-profile)) `[yes/no]`.  
The muon beam 3-momenta in the format *p<sub>x</sub> p<sub>y</sub> p<sub>z</sub>* must be piped into the file `'path'/beamprofile.fifo`. The three components are generic, not necessarity along the positive z-axis.  
**Notice:** when the muom beam is not along the positive z-axis, all the angles described below, used as generator level cuts, must be intended as relative to the muon beam direction
* `bspr [0]`: Gaussian beam energy spread, in % of `Ebeam`. Active only if `extmubeam no`
* `Eemin [1]`: minimum outgoing electron energy in GeV
* `themax [100]`: maximum outgoing electron angle in mrad
* `thmumin [0]`: minimum outgoing muon angle in mrad
* `thmumax [100]`: maximum outgoing muon angle in mrad
* `acoplcut [no 3.5]`: if the acoplanarity cut has to be applied or not `[yes/no mrad]`
* `elastcut [no 0.2]`: if the cut on the geometric distance from the elesticity curve in the [&theta;<sub>&mu;</sub>,&theta;<sub>e</sub>] plane has to be applied or not `[yes/no mrad]`
* `Ethr [0.2]`: minimum energy above which a lepton possibly triggers the detector (GeV)
* `ththr [100]`: lab. angle above which the detector is blind (mrad).  
  (*i.e.* a lepton with energy > `Ethr` and angle < `ththr` counts as a possible track, otherwise it is considered undetectable. The code requires that events passing the selection criteria have strictly two visible tracks.)

* `ord [alpha]`: to which order (photonic) corrections are included `[born/alpha/alpha2]` (standing for `[LO/NLO/NNLO]`)
* `arun [on]`: if vacuum polarization (VP) effects must be includer or not. Possible values are
  * `off`: no VP
  * `on`: both leptonic and hadronic VP are included. The hadronic part defaults to `hadr5`
  * `hadr5`: include leptonic VP and for hadronic VP use latest Fred Jegerlehner's parameterization
  * `nsk`: include leptonic VP and for hadronic VP use Fedor Ignatov's parameterization
  * `knt`: include leptonic VP and for hadronic VP use Keshavarzi-Nomura-Teubner's parameterization
* `hadoff [no]`: if hadronic VP has to be switched off `[yes/no]`
* `nev [10000000.]`: number of events to be generated
* `store [no old no]`: if events have to be stored `[yes/no]`, which storage version is to be used (`[old]` for v1, `[new]` for v2), if in the v2 event record the coefficients for the VP reweighting have to be written `[yes/no]`. See [below](#description-of-event-format) for a short description of the v1/v2 (old/new) event formats
* `storemode [0]`: if `store yes`, which mode to use to store events `[0/1/2]`
  * `0`: plain ASCII text file
  * `1`: `ROOT` format. **`MESMER`** concurrently runs `write-root-events` (a symbolic link to `root-interface/write_MuE_MCevents.exe`), developed by G. Abbiendi, which writes through a named pipe a `.root` file with the events
  * `2`: an `xz` compressed file is saved
  * `3`: just write events to the fifo file and wait for an external process to read them
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
* `sdmax [1e-10]`: maximum weight used for unweightening, when `mode unweighted` and `nwarmup 0`
* `wnorm [-1.]`: typical integrated cross section within applied cuts. This value is used for storing events: the true weight of the event *w<sub>t</sub>* is saved as *w = w<sub>t</sub> / wnorm* in order to make the average of the weights *w* over the generated sample close to 1
* `ndistr [1]`: number of distributions at different orders. For example, if running at NNLO (`ord alpha2`), the distribution files saved in `path` can be produced also at NLO and LO with the same run. The defaults `ndistr 1` produces distributions only at the selected order. If `ord alpha` and `ndistr >=2`, distributions are produced at NLO and LO. If `ord alpha2`,  if `ndistr 2` they are produced at NNLO and NLO, if `ndistr 3`, distributions at NNLO, NLO and LO are saved
* `sync [0]`: syncronization mode for random numbers `[0/1]`. This is used to make as correlated as possible two samples with same `seed` but different `Ebeam`. *This feature needs more cross-checks*

### Beam profile
If the code is run with `extmubeam yes`, it expects to read muon beam 3-momenta for each event from the named pipe `'path'/beamprofile.fifo`. As an example, a sample of 10<sup>5</sup> incoming muon momenta is provided in the file `beamprofile-example.txt.gz` (thanks to Mateusz Goncerz). To test it, run `mesmer` with `extmubeam yes` (and `nev` < 10<sup>5</sup>), which creates the pipe and waits to read muon beam 3-momenta, and issue in another shell `zcat beamprofile-example.txt.gz > 'path'/beamprofile.fifo`. Any provider of muon beam momenta must write on the same pipe, using the format *p<sub>x</sub> p<sub>y</sub> p<sub>z</sub>* (in GeV). If `store yes new yes|no`, the momenta of the stored particles will account for the generic direction of the incoming muon.

## Description of event format
The event format has been developed together with G. Abbiendi.

Version 1 is compliant with the current version of the `ROOT` file writer ([mantained by G. Abbiendi](https://gitlab.cern.ch/muesli/nlo-mc/mue/-/tree/master/writer)), but it is limited to parse only LO, NLO and NNLO events (*events with extra leptonic pairs are not supported yet*).  

### Version 2 (still under development)

The simple ASCII file (`storemode 0`) contains an header section between the tags `<header>` and `</header>` and a closing section between the tags `<footer>` and `</footer>`. After the header, each *event record* is enclosed between the `<event>` and `</event>` tags.  
The header contains useful info about the run and the set parameters. The footer contains some statistics of the generated sample.  
In the header, the first lines define the `SAMPLE TAG`, independent generations should have different `SAMPLE TAG`.

A typical *event record* looks like
```
 <event>
  42
  3
  4
   7.8095926596468930        7.6978162481410175        7.8080524876985837     
   0.0000000000000000        0.0000000000000000     
  0.0000000000000000 6.5364544453905386E-008 1.4984872179145210E-006 1.0217440564926995E-005 -3.5773619113353054E-004 4.3651960301712784E-004 5.9685156108649662E-004 7.7928861519833872E-002 1.1381502449965999E-003 1.5852826328852689E-002 7.6022089937820363
  -13   0.0000000000000000        0.0000000000000000        149.99996278769049     
  -13 -0.10037474236270096        3.5966016117402458E-005   139.27262373767550     
  11   1.2912320479053441E-002   1.7581244722062357E-003   1.4103729937323426     
  22   6.8958665161184984E-002  -2.9624268532748676E-003   7.4286187467893452     
  22   1.8503756722462530E-002   1.1683363649512295E-003   1.8883473092725849     
 </event>
```
After the `<event>` tag, the first line is the `seed` of the sample, the second line is the event number and the third is the number (`nfs`) of final state particles (in this case a &mu;, an *e* and two &gamma;s).  
In the fourth line, three weights *w* are listed: the weight with full VP effect, the one with VP switched off and the one with only leptonic VP effects (without hadronic VP): in this way, with a single run VP effects can be studied.  
The fifth line represents the weights *w<sub>LO</sub>* and *w<sub>NLO</sub>* which can be used to get LO distributions from a NLO sample and LO and NLO distributions from a NNLO sample (in this case they are `0` because it's a 4-body final state).  
The sixth line (present only if `store yes new yes` is selected, i.e. if the option to store coefficients for VP reweighting is chosen) represents 11 coefficients needed by the analysis tool to reweight the event when changing VP functions.  
The seventh  line is the incoming &mu;, in the format *id p<sub>x</sub> p<sub>y</sub> p<sub>z</sub>*, where *id* is the [PDG Monte Carlo code](https://pdg.lbl.gov/2021/web/viewer.html?file=%2F2021/reviews/rpp2020-rev-monte-carlo-numbering.pdf) for the particle and *p<sub>x</sub>*, *p<sub>y</sub>* and *p<sub>z</sub>* (in GeV) are the three-momentum components of the incoming muon. Its energy can be calculated with the muon mass which is stored in the header section.    
Finally, the last `nfs` lines represent the *id* and three-momenta (in GeV) of the final state &mu;, *e* and two &gamma;s respectively, again in the format *id p<sub>x</sub> p<sub>y</sub> p<sub>z</sub>*.  
The tag `</event>` closes the *event record*.

The upcoming `ROOT` file writer will read through a named pipe this event format and will write all the information into a `.root` file.

### Version 1 (default)
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
Finally, the last three lines represent the momenta (in GeV) of the final state &mu;, *e* and &gamma; respectively, in the format *E p<sub>x</sub> p<sub>y</sub> p<sub>z</sub>*.  
The tag `</event>` closes the *event record*.

The `ROOT` file writer reads through a named pipe this event format and writes all the information into a `.root` file.

##### Footnotes
KNT: Available upon request from the authors