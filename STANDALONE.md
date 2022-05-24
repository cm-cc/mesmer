In **standalone mode**, the code is run by issuing in the working directory the command `./mesmer`, which displays a command prompt
where parameters for the run can be set. Alternatively, an [input data card (`input-example`)](input-example) is provided and it can be fed as input
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
   [ store     ] events storage: no yes
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

By typing `help` at the prompt, a short description of the parameters that can be set is displayed (more details [here](README.md#running-parameters-description)):

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
 store     ---> (2 parameters) if events have to be stored [yes/no],
                if writing also coefficients for VP reweighting [yes/no]
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
             └> [-1000 ] all possible final states (except 1010)
             └> [ 1010 ] pi0 production (and decay in to two photons)
 nwrite    ---> files are dumped every nwrite generated events [int] (negative integer means files written every -nwrite seconds)
 nwarmup   ---> after nwarmup events, also unweighted generation is started [int]. Plays with wnorm.
 ndistr    ---> if writing distributions also at different orders [0/1/2/3]  (0 disables writing distributions)
 sdmax     ---> maximum integrand for unweightening
 wnorm     ---> typical integrated cross section within applied cuts, used for storage in ROOT format
 sync      ---> syncronization mode for random numbers (see README) [0/1]
```

The meaning of the parameters is described in the [README](README.md#running-parameters-description)

### Beam profile
In standalone mode, if the code is run with `extmubeam yes`, it expects to read muon beam 3-momenta for each event from the named pipe `'path'/beamprofile.fifo`. As an example, a sample of $10^5$ incoming muon momenta is provided in the file `beamprofile-example.txt.gz` (thanks to Mateusz Goncerz). To test it, run `mesmer` with `extmubeam yes` (and `nev` $< 10^5$), which creates the pipe and waits to read muon beam 3-momenta, and issue in another shell `zcat beamprofile-example.txt.gz > 'path'/beamprofile.fifo`. Any provider of muon beam momenta must write on the same pipe, using the format $p_x\ p_y\ p_z$ (in GeV). If events must be stored, the momenta of the particles will account for the generic direction of the incoming muon.

## Description of event format
The event format has been developed together with Giovanni Abbiendi. The code uses by default Version 2 described in the following. The version flag is hard-wired, but can be easily modified if needed (see [below](#version-1-old-superseded)).

<!--
Version 1 is compliant with the current version of the `ROOT` file writer ([mantained by Giovanni Abbiendi](https://gitlab.cern.ch/muesli/nlo-mc/mue/-/tree/master/writer)), but it is limited to parse only LO, NLO and NNLO events (*events with extra leptonic pairs are not supported yet*).  
-->

### Version 2 (default)

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
After the `<event>` tag, the first line is the `seed` of the sample, the second line is the event number and the third is the number (`nfs`) of final state particles (in this case a $\mu$, an $e$ and two $\gamma$'s).  
In the fourth line, three weights *w* are listed: the weight with full VP effect, the one with VP switched off and the one with only leptonic VP effects (without hadronic VP): in this way, with a single run VP effects can be studied.  
The fifth line represents the weights $w_{LO}$ and $w_{NLO}$ which can be used to get LO distributions from a NLO sample and LO and NLO distributions from a NNLO sample (in this case they are `0` because it's a 4-body final state).  
The sixth line (present only if `store yes yes` is selected, i.e. if the option to store coefficients for VP reweighting is chosen) represents 11 coefficients needed by the analysis tool to reweight the event when changing VP functions.  
The seventh  line is the incoming $\mu$, in the format $id\ p_x\ p_y\ p_z$, where $id$ is the [PDG Monte Carlo code](https://pdg.lbl.gov/2021/web/viewer.html?file=%2F2021/reviews/rpp2020-rev-monte-carlo-numbering.pdf) for the particle and $p_x$, $p_y$ and $p_z$ (in GeV) are the three-momentum components of the incoming muon. Its energy can be calculated with the muon mass which is stored in the header section.    
Finally, the last `nfs` lines represent the $id$ and three-momenta (in GeV) of the final state $\mu$, $e$ and two $\gamma$'s respectively, again in the format $id\ p_x\ p_y\ p_z$.  
The tag `</event>` closes the *event record*.

The `ROOT` file writer reads through a named pipe this event format and writes all the information into a `.root` file.

### `ROOT` file writer

<!-- To download and compile the `.root` file writer, issue the command `make rootwriter`, which clones [Giovanni's repository](https://gitlab.cern.ch/muesli/nlo-mc/mue/), compiles the writer and links it to `./write-root-events`, which in turn is used inside the code to write the `.root` file. -->
<!-- In order to clone Giovanni's repository, a CERN account is needed. In case you don't have one, please ask for a copy of the writer to one of the authors. -->

In order to write in `ROOT` format,  [Giovanni's file writer ](https://gitlab.cern.ch/muesli/nlo-mc/mue/) is needed. Please contact the authors for instructions

### Version 1 (old and superseded)

If needed, modify in `userinterface.F` the line saying `istorver = 2` to `istorver = 1` and recompile.

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
After the `<event>` tag, the first line is the `seed` of the sample, the second one is the event number and the third is the number of final state particles (in this case a $\mu$, an $e$ and a $\gamma$).  
In the fourth line, three weights *w* are listed: the weight with full VP effect, the one with VP switched off and the one with only leptonic VP effects (without hadronic VP): in this way, with a single run VP effects can be studied.  
The fifth line represents the weight $w_{LO}$ which can be used to get LO distributions from a NLO sample (in this case it is `0` because it's a 3-body final state).  
The sixth line is the incoming $\mu$ energy in GeV.  
Finally, the last three lines represent the momenta (in GeV) of the final state $\mu$, $e$ and $\gamma$ respectively, in the format $E\ p_x\ p_y\ p_z$.  
The tag `</event>` closes the *event record*.
