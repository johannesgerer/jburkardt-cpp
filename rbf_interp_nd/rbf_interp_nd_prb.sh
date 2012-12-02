#!/bin/bash
#
g++ -c -g -I/$HOME/include rbf_interp_nd_prb.cpp >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling rbf_interp_nd_prb.cpp"
  exit
fi
rm compiler.txt
#
g++ rbf_interp_nd_prb.o /$HOME/libcpp/$ARCH/rbf_interp_nd.o /$HOME/libcpp/$ARCH/r8lib.o -lm
if [ $? -ne 0 ]; then
  echo "Errors linking and loading rbf_interp_nd_prb.o"
  exit
fi
#
rm rbf_interp_nd_prb.o
#
mv a.out rbf_interp_nd_prb
./rbf_interp_nd_prb > rbf_interp_nd_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running rbf_interp_nd_prb."
  exit
fi
rm rbf_interp_nd_prb
#
echo "Program output written to rbf_interp_nd_prb_output.txt"
