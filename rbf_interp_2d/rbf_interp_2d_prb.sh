#!/bin/bash
#
g++ -c -I/$HOME/include rbf_interp_2d_prb.cpp
if [ $? -ne 0 ]; then
  echo "Errors compiling rbf_interp_2d_prb.cpp"
  exit
fi
#
g++ rbf_interp_2d_prb.o /$HOME/libcpp/$ARCH/rbf_interp_2d.o /$HOME/libcpp/$ARCH/test_interp_2d.o /$HOME/libcpp/$ARCH/r8lib.o -lm
if [ $? -ne 0 ]; then
  echo "Errors linking and loading rbf_interp_2d_prb.o."
  exit
fi
#
rm rbf_interp_2d_prb.o
#
mv a.out rbf_interp_2d_prb
./rbf_interp_2d_prb > rbf_interp_2d_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running rbf_interp_2d_prb."
  exit
fi
rm rbf_interp_2d_prb
#
echo "Program output written to rbf_interp_2d_prb_output.txt"
