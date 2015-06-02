#!/bin/bash
#
g++ -c -I/$HOME/include sphere_llq_grid_prb.cpp
if [ $? -ne 0 ]; then
  echo "Errors compiling sphere_llq_grid_prb.cpp"
  exit
fi
#
g++ -o sphere_llq_grid_prb sphere_llq_grid_prb.o /$HOME/libcpp/$ARCH/sphere_llq_grid.o -lm
if [ $? -ne 0 ]; then
  echo "Errors linking and loading sphere_llq_grid_prb.o."
  exit
fi
#
rm sphere_llq_grid_prb.o
#
./sphere_llq_grid_prb > sphere_llq_grid_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running sphere_llq_grid_prb."
  exit
fi
rm sphere_llq_grid_prb
#
echo "Program output written to sphere_llq_grid_prb_output.txt"
