#!/bin/bash
#
g++ -c -g -I/$HOME/include simplex_coordinates_prb.cpp >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling simplex_coordinates_prb.cpp"
  exit
fi
rm compiler.txt
#
g++ simplex_coordinates_prb.o /$HOME/libcpp/$ARCH/simplex_coordinates.o -lm
if [ $? -ne 0 ]; then
  echo "Errors linking and loading simplex_coordinates_prb.o."
  exit
fi
#
rm simplex_coordinates_prb.o
#
mv a.out simplex_coordinates_prb
./simplex_coordinates_prb > simplex_coordinates_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running simplex_coordinates_prb."
  exit
fi
rm simplex_coordinates_prb
#
echo "Program output written to simplex_coordinates_prb_output.txt"
