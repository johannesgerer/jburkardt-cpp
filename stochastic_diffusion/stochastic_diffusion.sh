#!/bin/bash
#
cp stochastic_diffusion.hpp /$HOME/include
#
g++ -c -g -I/$HOME/include stochastic_diffusion.cpp >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling stochastic_diffusion.cpp"
  exit
fi
rm compiler.txt
#
mv stochastic_diffusion.o ~/libcpp/$ARCH/stochastic_diffusion.o
#
echo "Library installed as ~/libcpp/$ARCH/stochastic_diffusion.o"
