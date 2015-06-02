#!/bin/bash
#
cp stochastic_heat2d.hpp /$HOME/include
#
g++ -c -g -I/$HOME/include stochastic_heat2d.cpp >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling stochastic_heat2d.cpp"
  exit
fi
rm compiler.txt
#
mv stochastic_heat2d.o ~/libcpp/$ARCH/stochastic_heat2d.o
#
echo "Library installed as ~/libcpp/$ARCH/stochastic_heat2d.o"
