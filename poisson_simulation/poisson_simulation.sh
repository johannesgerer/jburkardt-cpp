#!/bin/bash
#
cp poisson_simulation.hpp /$HOME/include
#
g++ -c -g -I /$HOME/include poisson_simulation.cpp >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling poisson_simulation.cpp"
  exit
fi
rm compiler.txt
#
mv poisson_simulation.o ~/libcpp/$ARCH/poisson_simulation.o
#
echo "Library installed as ~/libcpp/$ARCH/poisson_simulation.o"
