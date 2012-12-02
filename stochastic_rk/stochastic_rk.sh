#!/bin/bash
#
cp stochastic_rk.hpp /$HOME/include
#
g++ -c -g -I /$HOME/include stochastic_rk.cpp >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling stochastic_rk.cpp"
  exit
fi
rm compiler.txt
#
mv stochastic_rk.o ~/libcpp/$ARCH/stochastic_rk.o
#
echo "Library installed as ~/libcpp/$ARCH/stochastic_rk.o"
