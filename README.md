This repository contains the source code for the golem95 library. 

It can perform the tensor reduction of one-loop integrals for up to six external legs and also contains the one-loop master integrals. 

For integrals with massless propagators, it contains a method to evaluate tensor integrals directly via a fast one-dimensional numerical integration, this way avoiding the occurrence of inverse Gram determinants. This improves the numerical stability a lot in kinematic regions where the Gram determinant tends to zero.


# INSTALLATION
To install the library without linking LoopTools, type the following command:

```bash
./configure [--prefix=some_path]
  [--with-precision=quadruple] 
  [FC=your_fortran90_compiler]
  [F77=your_fortran77_compiler]
  [LDFLAGS=extra libraries]
make
make install
make clean [optional]
```

This installs a self-contained integral library containing all massive integrals (also for complex masses), 
where the finite scalar massive box/triangle D0/C0 is called from `avh_olo`.

If D0 and C0 should be taken from LoopTools instead, use 

```bash
./configure [--prefix=some_path]
  [--with-looptools=your_path_to_libooptools.a] 
  [FC=your_fortran90_compiler]
  [F77=your_fortran77_compiler]
make
make install
make clean [optional]
```

and set the variable "withlt" to `true` in your main file:
```fortran
withlt = .true. 
```
(this will overwrite the default `withlt = .false.` defined in `src/module/parametre.f90`)


Please note that some_path must be the absolute path of the directory 
where the library will be installed. It can also contain shell variables,
like in `--prefix=$HOME/my_libraries`



If the library should be compiled in quadruple precision, the 
option `--precision=quadruple` should be given when invoking 
configure.

If the option --prefix is not specified, the library will be installed 
to the system default library path. 
If the option compiler is not specified, the installation script
will try to find a suitable compiler.
If the option precision is not specified, the program will use double precision.

# USAGE

The program calculates tensor integrals as a combination 
of form factors and "basis" integrals with or without Feynman 
parameters in the numerator. 
It can be interfaced (by the user) to calculate either full amplitudes 
or only the integrals (with optional Feynman parameters in the numerator).

Examples of main programs, named `demo_npoint.f90`, 
which calculate these integrals, can be found in the `demos/` directory. To run them, go to the demos directory and type: `make`.

This will produce the executables for all included demo programs:

| Demo | Contained Features |
| --- | --- |
| `demo_3point` | three-point functions |
| `demo_4point` | four-point functions |
| `demo_5point` | five-point functions |
| `demo_6point` | six-point functions |
| `demo_4photon` | 4-photon helicity amplitudes |
| `demo_detg` | numerical stability: $\det(G) \rightarrow 0$ |
| `demo_a55_dets_sing` | numerical stability: $\det(S) \rightarrow 0$ |
| `demo_LT` | Golem call with LoopTools conventions |
| `demo_cmplx_masses.f9` | form factors with complex masses |
| `demo_tens_rec/mod.f90` | demo for use of the function `evaluate` |

Running the n-point examples will prompt you for input specifying 
- the number of non-light-like external legs
- the number of Feynman parameters in the numerator and 
  the space-time dimension 
 
Some arbitrary choice of the kinematic point (and of Feynman parameters) 
has been implemented in the demo programs.
The result will be written to the file `test_npoint.txt`.

If the user wants to define the numerical point and the numerator of 
the integral her/himself, she/he should go to the directory `test/` and 
follow the instructions in the `Readme` there.

# COMMON PROBLEMS

* The compilation fails due to lack of write permission:

  Typically, the user has forgotten to specify `--prefix` correctly.


* When compiling, the following error message occurs:   
    _Error in opening the compiled module file.  Check INCLUDE paths. \[PRECISION\]_

   The path to the location of the `.mod` files is not set correctly in the Makefile.
   The `.mod` files are in the subdirectory `/include/golem95`

* When running an example, the following error message occurs:      
   _Error while loading shared libraries: libgolem.so.0: cannot open shared object file: 
   No such file or directory_

    The path to the dynamically linked library is not known. 
    Either link `libgolem.a` statically or 
    add the path to `libgolem.so.0` to the `LD_LIBRARY_PATH` environment variable.
    For example, for bash shell:
    ```bash
    export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:${HOME}/golem95_v1.3.3/lib
    ```


* When running the example `demo_LT` in quadruple precision the result contains `NaN`

  This is the only example program which is not prepared for compilation in 
  quadruple precision. Valid results are only obtained in double precision.

* `make` does not work and complains about not being able to infer `--tag`

  This can happen when the configuration was (re-)generated by autotools
  with a version < 2.62. Please, either run `autoreconf` of a newer version
  or run `make LIBTOOLS_FLAGS=--tag=FC` instead of `make`.
