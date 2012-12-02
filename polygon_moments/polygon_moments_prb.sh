#!/bin/bash
#
g++ -c -g -I/$HOME/include polygon_moments_prb.cpp >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling polygon_moments_prb.cpp"
  exit
fi
rm compiler.txt
#
g++ polygon_moments_prb.o /$HOME/libcpp/$ARCH/polygon_moments.o -lm
if [ $? -ne 0 ]; then
  echo "Errors linking and loading polygon_moments_prb.o."
  exit
fi
#
rm polygon_moments_prb.o
#
mv a.out polygon_moments_prb
./polygon_moments_prb > polygon_moments_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running polygon_moments_prb."
  exit
fi
rm polygon_moments_prb
#
echo "Program output written to polygon_moments_prb_output.txt"
