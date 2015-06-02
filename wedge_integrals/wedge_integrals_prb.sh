#!/bin/bash
#
g++ -c -I/$HOME/include wedge_integrals_prb.cpp
if [ $? -ne 0 ]; then
  echo "Errors compiling wedge_integrals_prb.c"
  exit
fi
#
g++ -o wedge_integrals_prb wedge_integrals_prb.o /$HOME/libcpp/$ARCH/wedge_integrals.o -lm
if [ $? -ne 0 ]; then
  echo "Errors linking and loading wedge_integrals_prb.o."
  exit
fi
#
rm wedge_integrals_prb.o
#
./wedge_integrals_prb > wedge_integrals_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running wedge_integrals_prb."
  exit
fi
rm wedge_integrals_prb
#
echo "Program output written to wedge_integrals_prb_output.txt"
