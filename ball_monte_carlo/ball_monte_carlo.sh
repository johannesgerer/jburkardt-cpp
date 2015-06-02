#!/bin/bash
#
cp ball_monte_carlo.hpp /$HOME/include
#
g++ -c -g -I/$HOME/include ball_monte_carlo.cpp >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling ball_monte_carlo.cpp"
  exit
fi
rm compiler.txt
#
mv ball_monte_carlo.o ~/libcpp/$ARCH/ball_monte_carlo.o
#
echo "Library installed as ~/libcpp/$ARCH/ball_monte_carlo.o"
