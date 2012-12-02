#!/bin/bash
#
g++ -c -g -I/$HOME/include hex_grid_prb.cpp >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling hex_grid_prb.cpp"
  exit
fi
rm compiler.txt
#
g++ hex_grid_prb.o /$HOME/libcpp/$ARCH/hex_grid.o -lm
if [ $? -ne 0 ]; then
  echo "Errors linking and loading hex_grid_prb.o."
  exit
fi
#
rm hex_grid_prb.o
#
mv a.out hex_grid_prb
./hex_grid_prb > hex_grid_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running hex_grid_prb."
  exit
fi
rm hex_grid_prb
#
echo "Program output written to hex_grid_prb_output.txt"
