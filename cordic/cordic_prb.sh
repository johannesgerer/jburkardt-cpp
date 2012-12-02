#!/bin/bash
#
g++ -c -g -I/$HOME/include cordic_prb.cpp >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling cordic_prb.cpp"
  exit
fi
rm compiler.txt
#
g++ cordic_prb.o /$HOME/libcpp/$ARCH/cordic.o -lm
if [ $? -ne 0 ]; then
  echo "Errors linking and loading cordic_prb.o."
  exit
fi
#
rm cordic_prb.o
#
mv a.out cordic_prb
./cordic_prb > cordic_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running cordic_prb."
  exit
fi
rm cordic_prb
#
echo "Program output written to cordic_prb_output.txt"
