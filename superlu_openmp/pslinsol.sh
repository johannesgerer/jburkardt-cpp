#!/bin/bash
#
#  Compile
#
g++ -c -g -fopenmp -I/$HOME/include pslinsol.cpp >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling pslinsol.cpp"
  exit
fi
rm compiler.txt
#
#  Link and load
#
g++ -fopenmp pslinsol.o -L/$HOME/lib/$ARCH -L/$HOME/libc/$ARCH -lsuperlu_openmp -lm -lblas
if [ $? -ne 0 ]; then
  echo "Errors linking and loading pslinsol.o"
  exit
fi
rm pslinsol.o
mv a.out pslinsol
#
#  Run with 1 processor.
#
export OMP_NUM_THREADS=1
./pslinsol < g20_rua.txt > pslinsol_1_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running pslinsol"
  exit
fi
#
#  Run with 4 processors.
#
export OMP_NUM_THREADS=4
./pslinsol < g20_rua.txt > pslinsol_4_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running pslinsol"
  exit
fi
#
rm pslinsol
#
#  Terminate.
#
echo "Program output written to pslinsol_output.txt"
