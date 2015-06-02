#!/bin/bash
#
g++ -c -I/$HOME/include pyramid_grid_prb.cpp
if [ $? -ne 0 ]; then
  echo "Errors compiling pyramid_grid_prb.cpp"
  exit
fi
#
g++ -o pyramid_grid_prb pyramid_grid_prb.o /$HOME/libcpp/$ARCH/pyramid_grid.o -lm
if [ $? -ne 0 ]; then
  echo "Errors linking and loading pyramid_grid_prb.o."
  exit
fi
#
rm pyramid_grid_prb.o
#
./pyramid_grid_prb > pyramid_grid_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running pyramid_grid_prb."
  exit
fi
rm pyramid_grid_prb
#
echo "Program output written to pyramid_grid_prb_output.txt"
