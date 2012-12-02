#!/bin/bash
#
g++ -c -g -I/$HOME/include wavelet_prb.cpp >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling wavelet_prb.cpp"
  exit
fi
rm compiler.txt
#
g++ wavelet_prb.o /$HOME/libcpp/$ARCH/wavelet.o -lm
if [ $? -ne 0 ]; then
  echo "Errors linking and loading wavelet_prb.o."
  exit
fi
#
rm wavelet_prb.o
#
mv a.out wavelet_prb
./wavelet_prb > wavelet_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running wavelet_prb."
  exit
fi
rm wavelet_prb
#
echo "Program output written to wavelet_prb_output.txt"
