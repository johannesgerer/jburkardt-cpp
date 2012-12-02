#!/bin/bash
#
cp sobol.hpp /$HOME/include
#
g++ -c -g -I /$HOME/include sobol.cpp >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling sobol.cpp"
  exit
fi
rm compiler.txt
#
mv sobol.o ~/libcpp/$ARCH/sobol.o
#
echo "Library installed as ~/libcpp/$ARCH/sobol.o"
