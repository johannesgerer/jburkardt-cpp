#!/bin/bash
#
g++ -c multigrid_poisson_1d_prb.cpp
if [ $? -ne 0 ]; then
  echo "Errors compiling multigrid_poisson_1d_prb.cpp."
  exit
fi
#
g++ multigrid_poisson_1d_prb.o /$HOME/libcpp/$ARCH/multigrid_poisson_1d.o -lm
if [ $? -ne 0 ]; then
  echo "Errors linking and loading multigrid_poisson_1d_prb.o."
  exit
fi
#
rm multigrid_poisson_1d_prb.o
#
mv a.out multigrid_poisson_1d_prb
./multigrid_poisson_1d_prb > multigrid_poisson_1d_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running multigrid_poisson_1d_prb."
  exit
fi
rm multigrid_poisson_1d_prb
#
echo "Program output written to multigrid_poisson_1d_prb_output.txt"
