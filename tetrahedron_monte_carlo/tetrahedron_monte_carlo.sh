#!/bin/bash
#
cp tetrahedron_monte_carlo.hpp /$HOME/include
#
g++ -c -g -I /$HOME/include tetrahedron_monte_carlo.cpp >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling tetrahedron_monte_carlo.cpp"
  exit
fi
rm compiler.txt
#
mv tetrahedron_monte_carlo.o ~/libcpp/$ARCH/tetrahedron_monte_carlo.o
#
echo "Library installed as ~/libcpp/$ARCH/tetrahedron_monte_carlo.o"
