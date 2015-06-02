#!/bin/bash
#
cp pyramid_monte_carlo.hpp /$HOME/include
#
g++ -c -I /$HOME/include pyramid_monte_carlo.cpp >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling pyramid_monte_carlo.cpp"
  exit
fi
rm compiler.txt
#
mv pyramid_monte_carlo.o ~/libcpp/$ARCH/pyramid_monte_carlo.o
#
echo "Library installed as ~/libcpp/$ARCH/pyramid_monte_carlo.o"
