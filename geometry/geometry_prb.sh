#!/bin/bash
#
g++ -c -g -I/$HOME/include geometry_prb.cpp >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling geometry_prb.cpp"
  exit
fi
rm compiler.txt
#
g++ geometry_prb.o /$HOME/libcpp/$ARCH/geometry.o -lm
if [ $? -ne 0 ]; then
  echo "Errors linking and loading geometry_prb.o."
  exit
fi
#
rm geometry_prb.o
#
mv a.out geometry_prb
./geometry_prb > geometry_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running geometry_prb."
  exit
fi
rm geometry_prb
#
echo "Program output written to geometry_prb_output.txt"
