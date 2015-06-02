#!/bin/bash
#
g++ -c -I/$HOME/include van_der_corput_prb.cpp
if [ $? -ne 0 ]; then
  echo "Errors compiling van_der_corput_prb.cpp"
  exit
fi
#
g++ van_der_corput_prb.o /$HOME/libcpp/$ARCH/van_der_corput.o -lm
if [ $? -ne 0 ]; then
  echo "Errors linking and loading van_der_corput_prb.o."
  exit
fi
#
rm van_der_corput_prb.o
#
mv a.out van_der_corput_prb
./van_der_corput_prb > van_der_corput_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running van_der_corput_prb."
  exit
fi
rm van_der_corput_prb
#
echo "Program output written to van_der_corput_prb_output.txt"
