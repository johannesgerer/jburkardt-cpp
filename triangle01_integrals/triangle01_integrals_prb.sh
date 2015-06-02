#!/bin/bash
#
g++ -c -I/$HOME/include triangle01_integrals_prb.cpp
if [ $? -ne 0 ]; then
  echo "Errors compiling triangle01_integrals_prb.cpp"
  exit
fi
#
g++ triangle01_integrals_prb.o /$HOME/libcpp/$ARCH/triangle01_integrals.o -lm
if [ $? -ne 0 ]; then
  echo "Errors linking and loading triangle01_integrals_prb.o."
  exit
fi
#
rm triangle01_integrals_prb.o
#
mv a.out triangle01_integrals_prb
./triangle01_integrals_prb > triangle01_integrals_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running triangle01_integrals_prb."
  exit
fi
rm triangle01_integrals_prb
#
echo "Program output written to triangle01_integrals_prb_output.txt"
