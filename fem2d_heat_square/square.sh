#!/bin/bash
#
g++ -c -g square.cpp >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling square.cpp"
  exit
fi
rm compiler.txt
#
#  Link the precompiled main program FEM2D_HEAT with the user routines.
#
g++ ~/libcpp/$ARCH/fem2d_heat.o square.o
if [ $? -ne 0 ]; then
  echo "Errors linking and loading free_fem_heat.o + square.o"
  exit
fi
rm square.o
#
#  Run the program with the user mesh files.
#
chmod ugo+x a.out
mv a.out square
square square_nodes.txt square_elements.txt > square_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running square."
  exit
fi
rm square
#
echo "The square program has been executed."
