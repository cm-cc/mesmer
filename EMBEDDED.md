In **embedded mode**, the code is supposed to be included into a `C/C++` program, which must be linked against the library `libmesmerfull.a` (together with system `libm.[a,so]` and `libgfortran.[a,so]`).

An extremely simple `C` example is provided by the file [`c_driver.c`](c_driver.c); it can be compiled for instance with the command
`gcc c_driver.c -L. -lmesmerfull  -lm -lgfortran -o mesmer_from_c`.

The embedded mode has been tested and developed in parallel with the [fast simulation & analysis code](https://gitlab.cern.ch/muesli/nlo-mc/mue/-/tree/master/writer) by Giovanni Abbiendi.
