#!/bin/bash
#
g++ -c -g -I/$HOME/include sphere_stereograph_prb.cpp >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling sphere_stereograph_prb.cpp"
  exit
fi
rm compiler.txt
#
g++ sphere_stereograph_prb.o /$HOME/libcpp/$ARCH/sphere_stereograph.o -lm
if [ $? -ne 0 ]; then
  echo "Errors linking and loading sphere_stereograph_prb.o."
  exit
fi
#
rm sphere_stereograph_prb.o
#
mv a.out sphere_stereograph_prb
./sphere_stereograph_prb > sphere_stereograph_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running sphere_stereograph_prb."
  exit
fi
rm sphere_stereograph_prb
#
echo "Program output written to sphere_stereograph_prb_output.txt"
