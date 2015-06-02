#!/bin/bash
#
g++ -c -I/$HOME/include fem1d_lagrange_prb.cpp
if [ $? -ne 0 ]; then
  echo "Errors compiling fem1d_lagrange_prb.cpp"
  exit
fi
#
g++ fem1d_lagrange_prb.o /$HOME/libcpp/$ARCH/fem1d_lagrange.o -lm
if [ $? -ne 0 ]; then
  echo "Errors linking and loading fem1d_lagrange_prb.o"
  exit
fi
#
rm fem1d_lagrange_prb.o
#
mv a.out fem1d_lagrange_prb
./fem1d_lagrange_prb > fem1d_lagrange_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running fem1d_lagrange_prb."
  exit
fi
rm fem1d_lagrange_prb
#
echo "Program output written to fem1d_lagrange_prb_output.txt"
