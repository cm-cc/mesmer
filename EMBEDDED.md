In **embedded mode**, the code is supposed to be included into a `C/C++` program, which must be linked against the library `libmesmerfull.a` (together with standard system `libm` and `libgfortran`).

A simple `C` example is provided by the file [`c_driver.c`](c_driver.c); to produce the executable `mesmer_from_c`, it can be compiled for instance with the command  
`gcc c_driver.c -L. -lmesmerfull  -lm -lgfortran -o mesmer_from_c`.

The embedded mode has been tested and developed in parallel with the [fast simulation & analysis code](https://gitlab.cern.ch/muesli/nlo-mc/mue/-/tree/master/) by Giovanni Abbiendi.

The calling program must call the following functions:

* `int init_mesmer(char* input_data_card)`  
Must be called only once, mandatory.  
It reads the input parameters from a data card and initializes **`MESMER`**. It returns `0` if initialization is successful, `1` otherwise

* `void mesmer_setup
    (int* sampletag,
     char* mesmerversion,
     char* hostname,
     char* datetime,
     int* idproc,
     long int* nev,
     int* areweighted,
     char* RCorder,
     int* includedfs,
     int* radmu,
     int* rade,
     int* iseed1,
     int* iseed2,
     double* emulab,
     double* spread,
     int* extmubeam,
     double* Qmu,
     double* mumass,
     double* elmass,
     double* invalpha,
     double* wnorm,
     double* wmax,
     double* eemin,
     double* themin,
     double* themax,
     double* thmumin,
     double* thmumax,
     double* ethr,
     double* ththr,
     int* iacopl,
     double* acopl,
     int* iela,
     double* ela,
     int* ivpwgts,
     int* ihadon,
     int* ivpfl,
     int* nwarmup,
     int* ndistrw,
     int* isync,
     double* eps,
     double* phmass
     )`  
     Optional.  
     It returns the parameters of the run

* `void IncomingMuonMomentum_mesmer(double* pmu)`  
Must be called once for each event, in the do/for loop over events, mandatory.  
It returns the generic incoming muon 4-momentum.  
It can be replaced by any function returning the incoming muon 4-momentum as `double pmu[4]`, with `p[0]` $= E$, `p[1]` $= p_x$, `p[2]` $= p_y$,  `p[3]` $= p_z$

* `void generate_event_mesmer(double *pmu, int *nfs, int *mcids, double (*pmat)[4], double *weight,
				    int *itag, long int *ievtnr, double *wnovp, double *wnohad, double *wLO,
				    double *wNLO, double *cwvp, int *ierr)`  
Must be called once for each event, in the do/for loop over events, mandatory.  
It generate the event and returns the 4-momenta of the particles, their [PDG Monte Carlo codes](https://pdg.lbl.gov/2021/web/viewer.html?file=%2F2021/reviews/rpp2020-rev-monte-carlo-numbering.pdf), various weights of the event and an error code

* `void finalize_mesmer(double *xsw, double *exsw, long int *foohpm, long int *fooh, double *truemax, long int *nabove,
			       long int *nlt0, double *xsbias, double *exsbias, double *xsbiasn,
			       double *exsbiasn, double *sumow, double *sum2ow2,
			       double *sumnow, double *sum2now2)`  
Must be called only once, after the do/for loop over events, mandatory.  
It returns information of the run and finalizes **`MESMER`**

The declarations of all the argument variables can be more easily checked in the `C` example [`c_driver.c`](c_driver.c): they essentially expose the information which is written in the event file when running in standalone mode.  

Otherwise, these functions are wrapped inside functions and classes of the `C++` code by Giovanni Abbiendi, in a transparent way to the user.

Soon a short description of the variables will be included also here.