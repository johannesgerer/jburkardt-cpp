#!/bin/bash
#
g++ -c -I/$HOME/include hyperball_integrals_prb.cpp
if [ $? -ne 0 ]; then
  echo "Errors compiling hyperball_integrals_prb.cpp"
  exit
fi
#
g++ hyperball_integrals_prb.o /$HOME/libcpp/$ARCH/hyperball_integrals.o -lm
if [ $? -ne 0 ]; then
  echo "Errors linking and loading hyperball_integrals_prb.o."
  exit
fi
#
rm hyperball_integrals_prb.o
#
mv a.out hyperball_integrals_prb
./hyperball_integrals_prb > hyperball_integrals_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running hyperball_integrals_prb."
  exit
fi
rm hyperball_integrals_prb
#
echo "Program output written to hyperball_integrals_prb_output.txt"
