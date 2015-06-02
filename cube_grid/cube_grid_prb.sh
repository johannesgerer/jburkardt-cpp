#!/bin/bash
#
g++ -c -I/$HOME/include cube_grid_prb.cpp
if [ $? -ne 0 ]; then
  echo "Errors compiling cube_grid_prb.cpp"
  exit
fi
#
g++ -o cube_grid_prb cube_grid_prb.o /$HOME/libcpp/$ARCH/cube_grid.o -lm
if [ $? -ne 0 ]; then
  echo "Errors linking and loading cube_grid_prb.o."
  exit
fi
#
rm cube_grid_prb.o
#
./cube_grid_prb > cube_grid_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running cube_grid_prb."
  exit
fi
rm cube_grid_prb
#
echo "Program output written to cube_grid_prb_output.txt"
