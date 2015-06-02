#!/bin/bash
#
g++ -c -I/$HOME/include polygon_grid_prb.cpp
if [ $? -ne 0 ]; then
  echo "Errors compiling polygon_grid_prb.cpp"
  exit
fi
#
g++ -o polygon_grid_prb polygon_grid_prb.o /$HOME/libcpp/$ARCH/polygon_grid.o -lm
if [ $? -ne 0 ]; then
  echo "Errors linking and loading polygon_grid_prb.o."
  exit
fi
#
rm polygon_grid_prb.o
#
./polygon_grid_prb > polygon_grid_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running polygon_grid_prb."
  exit
fi
rm polygon_grid_prb
#
echo "Program output written to polygon_grid_prb_output.txt"
