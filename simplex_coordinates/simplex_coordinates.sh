#!/bin/bash
#
cp simplex_coordinates.hpp /$HOME/include
#
g++ -c -g -I /$HOME/include simplex_coordinates.cpp >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling simplex_coordinates.cpp"
  exit
fi
rm compiler.txt
#
mv simplex_coordinates.o ~/libcpp/$ARCH/simplex_coordinates.o
#
echo "Library installed as ~/libcpp/$ARCH/simplex_coordinates.o"
