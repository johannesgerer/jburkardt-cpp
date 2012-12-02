#!/bin/bash
#
g++ -c -g -I/$HOME/include cell_prb.cpp >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling cell_prb.cpp"
  exit
fi
rm compiler.txt
#
g++ cell_prb.o /$HOME/libcpp/$ARCH/cell.o -lm
if [ $? -ne 0 ]; then
  echo "Errors linking and loading cell_prb.o"
  exit
fi
#
rm cell_prb.o
#
mv a.out cell_prb
./cell_prb > cell_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running cell_prb."
  exit
fi
rm cell_prb
#
echo "Program output written to cell_prb_output.txt"
