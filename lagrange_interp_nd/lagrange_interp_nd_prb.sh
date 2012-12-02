#!/bin/bash
#
g++ -c -g lagrange_interp_nd_prb.cpp >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling lagrange_interp_nd_prb.cpp"
  exit
fi
rm compiler.txt
#
g++ lagrange_interp_nd_prb.o /$HOME/libcpp/$ARCH/lagrange_interp_nd.o -lm
if [ $? -ne 0 ]; then
  echo "Errors linking and loading lagrange_interp_nd_prb.o."
  exit
fi
#
rm lagrange_interp_nd_prb.o
#
mv a.out lagrange_interp_nd_prb
./lagrange_interp_nd_prb > lagrange_interp_nd_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running lagrange_interp_nd_prb."
  exit
fi
rm lagrange_interp_nd_prb
#
echo "Program output written to lagrange_interp_nd_prb_output.txt"
