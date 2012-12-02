#!/bin/bash
#
cp ball_grid.hpp /$HOME/include
#
g++ -c -g ball_grid.cpp >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling ball_grid.cpp"
  exit
fi
rm compiler.txt
#
mv ball_grid.o ~/libcpp/$ARCH/ball_grid.o
#
echo "Library installed as ~/libcpp/$ARCH/ball_grid.o"
