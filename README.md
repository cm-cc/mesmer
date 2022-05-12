# MESMER

**`MESMER`** (**M**uon **E**lectron **S**cattering with **M**ultiple **E**lectromagnetic **R**adiation) is a Monte Carlo event generator for
high-precision simulation of muon-electron scattering at low enegies, developed for the MUonE experiment (see [here](https://web.infn.it/MUonE/) and/or [here](https://twiki.cern.ch/twiki/bin/view/MUonE/WebHome)).

### Authors
**`MESMER`** is developed at INFN, Sezione di Pavia, and Department of Physics, Università di Pavia (Italy).  
Authors are listed in the [`AUTHORS`](AUTHORS.md) file.

### Citing the code
We'd be grateful if you could cite (at least a subset of) the following papers when using the **`MESMER`** generator:  
1. [Alacevich *et al.*, Muon-electron scattering at NLO, JHEP 02 (2019) 155](https://inspirehep.net/literature/1703989)
2. [Carloni Calame *et al.*, Towards muon-electron scattering at NNLO, JHEP 11 (2020) 028](https://inspirehep.net/literature/1805205)  
3. [E. Budassi *et al.*, NNLO virtual and real leptonic corrections to muon-electron scattering, JHEP 11 (2021) 098](https://inspirehep.net/literature/1933852)
4. [E. Budassi *et al.*, Single *&pi;<sup>0</sup>* production in *&mu;e* scattering at MUonE, PLB 829 (2022) 137138](https://inspirehep.net/literature/2044898)


## Prerequisites & Compilation
The program is mainly written in `Fortran 77` and it has been extensively tested with `GCC` compilers on `GNU/Linux`.

Three external libraries are used, [`LoopTools`](http://www.feynarts.de/looptools/), [`Collier`](https://collier.hepforge.org/) and
the `C` implementation of [`RANLUX`](https://luscher.web.cern.ch/luscher/ranlux/), shipped with the code under GPL-like licences. To be compiled,
the `Collier` library requires the availability of `cmake`.

Furthermore, the program includes three different parameterization for the hadronic vacuum polarization: one by [Fred Jegerlehner](http://www-com.physik.hu-berlin.de/~fjeger/software.html), one by [Keshavarzi-Nomura-Teubner](#footnotes) and one by [Fedor Ignatov](http://cmd.inp.nsk.su/~ignatov/vpl/).

Writing events in ROOT format requires the [`ROOT`](https://root.cern/) framework to be installed.

A `Makefile` is provided. The executable `mesmer` and the library `libmesmerfull.a` are build by simply issuing the command `make`.

**NOTICE: as of now, the released code runs only at LO, NLO and NNLO. If NNLO extra leptonic pair production are selected, the execution stops with a warning.** 

## Warning
**When running at NNLO, please understand that a subset of the purely virtual NNLO photonic corrections are approximate as described in 2. of the list of references, i.e. the subset of virtual corrections where at least two photons connect the muon and the electron fermionic lines.** 

## Using the code

Once compiled, **`MESMER`** can be run in [**standalone mode**](STANDALONE.md) or [**embedded mode**](EMBEDDED.md) (i.e. inside a `C/C++` program).

## Running parameters description

In general, the routine `cuts(...)` in the file `cuts.F` can be modified according to the needs of the user. The default one applies the selection criteria described below (generator level cuts). The only mandatory generator level cut is `Eemin`, which **must** be stricty larger than *m<sub>e</sub>*. The other cuts can be set at their kinematical limits.

The parameters that can be set are split into *principal* and *internal* parameters and they are set by typing at the prompt `parameter value(s)` (case sensitive).

The order in which the input parameters and values are inserted in the data card or at the prompt is indifferent.
If a parameter/value is missing, defaults are used. The only mandatory rule is that the last input must be `run`.

### Principal parameters

* `Qmu [default 1]`: charge of the incoming muon in e+ units `[1 / -1]`
* `Ebeam [150]`: nominal incoming muon energy in GeV
* `extmubeam [no]`: feed generic muon beam externally (beam profile, see [here if running in standalone mode](STANDALONE.md#beam-profile) or [here if running in embedded mode](EMBEDDED.md)) `[yes/no]`.  
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
* `store [no yes]`: if events have to be stored `[yes/no]` and <!-- which storage version is to be used (`[old]` for v1, `[new]` for v2),-->
if the coefficients for the VP reweighting have to be written in the event record `[yes/no]`.  
*This parameter is actually used only in standalone mode*, because in embedded mode it's assumed that event storage is managed by the calling program.  
See [here](STANDALONE.md#description-of-event-format) for a description of the event formats
* `storemode [0]`: if `store yes`, which mode to use to store events `[0/1/2/3]`
  * `0`: plain ASCII text file
  * `1`: `ROOT` format. **`MESMER`** concurrently runs `write-root-events`, developed by Giovanni Abbiendi, which writes through a named pipe a `.root` file with the events.  
  If needed, please contact the authors for instuctions
  * `2`: an `xz` compressed file is saved
  * `3`: just write events to the fifo file and wait for an external process to read them
* `path [test-run/]`: the directory where all outputs are saved. It will contain some `.txt` files with differential distributions of some variables, the file `events.[dat,dat.xz,root]` with saved events if `store yes` and the file `stat_*.txt`, where cross sections and all info of the current run are reported
* `seed [42]`: seed for the pseudo-random-number-generator, it must be set to a "small" integer. Independent generations must use different seeds

### Inner parameters

* `mode [weighted]`: if stored events are weighted or unweighted `[weighted/unweighted]`.  
If `weighted`, to each event is associated a different weight, which must be appropriately accounted for when producing distributions.  
If `unweighted`, all events have the same constant weight, *i.e.* they are directly distributed according to the underlying cross sections.  
**Notice:** unweighted generation can be much slower, due to the unweightening procedure
* `nwrite [500000]`: files in `path` are dumped every `nwrite` generated events. If `nwrite < 0`, files are dumped (approximately) every `-nwrite` seconds
* `nwarmup [0]`:  if `mode weighted`, the maximum weight for subsequent unweightening is searched generating first `nwarmup` events, after which unweighted generation is started. Notice that the maximum weight is a guessed value.  
If `store yes` and `nwarmup > 0`, `wnorm` is calculated automatically throwing `nwarmup` events.  
For unweightening, the maximum weight `sdmax` can be alternatively set by hand (setting `nwarmup 0` at the same time)
* `sdmax [1e-10]`: maximum weight used for unweightening, when `mode unweighted` and `nwarmup 0`
* `wnorm [-1.]`: typical integrated cross section within applied cuts. This value is used for storing events: the true weight of the event *w<sub>t</sub>* is saved as *w = w<sub>t</sub> / wnorm* in order to make the average of the weights *w* over the generated sample close to 1
* `sync [0]`: syncronization mode for random numbers `[0/1]`. This is used to make as correlated as possible two samples with same `seed` but different `Ebeam`.
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
  * `1010`: simulates *&mu;<sup>&plusmn;</sup>e<sup>-</sup>&rarr;&mu;<sup>&plusmn;</sup>e<sup>-</sup>&pi;<sup>0</sup>* with *&pi;<sup>0</sup>&rarr;&gamma;&gamma;*
* `ndistr [1]`: number of distributions at different orders.  
For example, if running at NNLO (`ord alpha2`), the distribution files saved in `path` can be produced also at NLO and LO with the same run.  
The defaults `ndistr 1` produces distributions only at the selected order.  
If `ord alpha` and `ndistr >=2`, distributions are produced at NLO and LO.  
If `ord alpha2`,  if `ndistr 2` they are produced at NNLO and NLO, if `ndistr 3`, distributions at NNLO, NLO and LO are saved.  
If `ndistr 0`, no distribution files are saved.

##### Footnotes
KNT: Available upon request from the authors
