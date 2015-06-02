#!/bin/bash
#
g++ -c -I/$HOME/include lagrange_nd_prb.cpp
if [ $? -ne 0 ]; then
  echo "Errors compiling lagrange_nd_prb.cpp"
  exit
fi
#
g++ lagrange_nd_prb.o /$HOME/libcpp/$ARCH/lagrange_nd.o -lm
if [ $? -ne 0 ]; then
  echo "Errors linking and loading lagrange_nd_prb.o"
  exit
fi
#
rm lagrange_nd_prb.o
#
mv a.out lagrange_nd_prb
./lagrange_nd_prb > lagrange_nd_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running lagrange_nd_prb."
  exit
fi
rm lagrange_nd_prb
#
echo "Program output written to lagrange_nd_prb_output.txt"
