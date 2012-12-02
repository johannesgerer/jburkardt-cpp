#!/bin/bash
#
g++ -c -g -I/$HOME/include sphere_grid_prb.cpp >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling sphere_grid_prb.cpp"
  exit
fi
rm compiler.txt
#
g++ sphere_grid_prb.o /$HOME/libcpp/$ARCH/sphere_grid.o -lm
if [ $? -ne 0 ]; then
  echo "Errors linking and loading sphere_grid_prb.o."
  exit
fi
#
rm sphere_grid_prb.o
#
mv a.out sphere_grid_prb
./sphere_grid_prb > sphere_grid_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running sphere_grid_prb."
  exit
fi
rm sphere_grid_prb
#
mv *.xyz ../../datasets/sphere_grid
echo "XYZ data files moved to ../../datasets/sphere_grid directory."
#
echo "Program output written to sphere_grid_prb_output.txt"
