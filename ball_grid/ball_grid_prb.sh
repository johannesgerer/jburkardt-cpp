#!/bin/bash
#
g++ -c -g ball_grid_prb.cpp >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling ball_grid_prb.cpp"
  exit
fi
rm compiler.txt
#
g++ ball_grid_prb.o /$HOME/libcpp/$ARCH/ball_grid.o -lm
if [ $? -ne 0 ]; then
  echo "Errors linking and loading ball_grid_prb.o."
  exit
fi
#
rm ball_grid_prb.o
#
mv a.out ball_grid_prb
./ball_grid_prb > ball_grid_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running ball_grid_prb."
  exit
fi
rm ball_grid_prb
#
echo "Program output written to ball_grid_prb_output.txt"
