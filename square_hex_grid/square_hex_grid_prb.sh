#!/bin/bash
#
g++ -c -I/$HOME/include square_hex_grid_prb.cpp
if [ $? -ne 0 ]; then
  echo "Errors compiling square_hex_grid_prb.cpp"
  exit
fi
#
g++ square_hex_grid_prb.o /$HOME/libcpp/$ARCH/square_hex_grid.o -lm
if [ $? -ne 0 ]; then
  echo "Errors linking and loading square_hex_grid_prb.o."
  exit
fi
#
rm square_hex_grid_prb.o
#
mv a.out square_hex_grid_prb
./square_hex_grid_prb > square_hex_grid_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running square_hex_grid_prb."
  exit
fi
rm square_hex_grid_prb
#
echo "Program output written to square_hex_grid_prb_output.txt"
