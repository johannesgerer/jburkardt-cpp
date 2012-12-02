#!/bin/bash
#
g++ -c -g -I/$HOME/include grid_prb.cpp >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling grid_prb.cpp"
  exit
fi
rm compiler.txt
#
g++ grid_prb.o /$HOME/libcpp/$ARCH/grid.o -lm
if [ $? -ne 0 ]; then
  echo "Errors linking and loading grid_prb.o."
  exit
fi
#
rm grid_prb.o
#
mv a.out grid_prb
./grid_prb > grid_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running grid_prb."
  exit
fi
rm grid_prb
#
echo "Program output written to grid_prb_output.txt"
