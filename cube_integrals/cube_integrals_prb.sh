#!/bin/bash
#
g++ -c -g -I/$HOME/include cube_integrals_prb.cpp >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling cube_integrals_prb.cpp"
  exit
fi
rm compiler.txt
#
g++ cube_integrals_prb.o /$HOME/libcpp/$ARCH/cube_integrals.o -lm
if [ $? -ne 0 ]; then
  echo "Errors linking and loading cube_integrals_prb.o."
  exit
fi
#
rm cube_integrals_prb.o
#
mv a.out cube_integrals_prb
./cube_integrals_prb > cube_integrals_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running cube_integrals_prb."
  exit
fi
rm cube_integrals_prb
#
echo "Program output written to cube_integrals_prb_output.txt"
