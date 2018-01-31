Program retrieved at the following URL:

http://cpc.cs.qub.ac.uk/summaries/ADZA_v1_0.html

********************************
*   NMscatt program package
*
*   Author: Franci Merzel
*           franc@cmm.ki.si
********************************

The source code including makefile and run_sctipt for the 
test case is given in directory ./source. 
Because of the storage requirements we selected the fragment of 
collagen to provide as the benchmark system in the 
distribution package. 

The input/output data are given in directories
./in and ./out, respectively.

################################################################################
CPC Librarian's Note
################################################################################

Compilation
------------
The makefile should be changed to reflect the name of the Fortran compiler 
you wish to use. It also assumes that the Lapack library is in the same 
directory as the program code. If this is not the case, the makefile must be 
changed to include the Lapack library path. 

In the source directory type
  make
this should produce 4 executables, bead, coh, incoh and photon.

Execution
----------
To run the benchmark test, copy the files in the directory "in", into the 
source directory. Type

  run_script  
  
 The files output can then be checked against those in the "out" directory.

The test run took 15 minutes to complete on a Itanium 2 processor, CPU 1300MHz

