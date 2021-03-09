How to run swave_sample.f90
===========================

This is a sample FORTRAN90 code to calculate the mean-field solution
for an s-wave pairing state around an impurity.

The code should be fairly self-explanatory. The basic idea is

1. take initial guess of pairing potential

2. construct Hamiltonian

3. solve eigenproblem

4. recalculate pairing potential

5. compare new pairing potential with old pairing potential; if
difference is not smaller than tol, repeat 2-5 until this condition is
satisfied. 

6. output data

On my machine this takes about 9 minutes to run for Nlat=31; for
Nlat=41 it requires about 45 minutes to complete. 


The code consists of 3 parts:

* swave_sample.f90 (the main code)

* diag_mod.f90 (code containing the diagonalization routine)

* nrtype.f90 (code defining some standard numerical objects; I'm not
  sure how necessary this is, but everyone seems to use it. It comes
  from a very influential book "Numerical Recipes")

To construct and run the program, you must compile these three files
together. This is done using the Makefile.

The Makefile compiles the code using gfortran (i.e. the gnu
compiler). MACOSX probably has a different native compiler that you
will need to use, if you wish to run the program itself.

Note that the compilation of the program requires linking to external
libraries "lapack" and "blas", since the subroutine eign calls the
zheevd subroutine of lapack, which also relies on blas. This is very
important, as otherwise the code won't run.

To use the Makefile, you should convert it to an executable. On a unix
system you should be able to do this by

chmod +x ./Makefile

at the command line. Then, to run the Makefile, just type

make

at the command line. 

The code outputs the converged pairing potential as the file
"potential.dat", which can then be plotted. I use the gnuplot plotting
tool, but if you use Julia or Matlab you will have this done
natively. If you use gnuplot, you can plot the data by using the
script 'plot' and typing

load 'plot'

at the gnuplot command line.

