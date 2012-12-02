#!/bin/bash
#
cp black_scholes.hpp /$HOME/include
#
g++ -c -g -I /$HOME/include black_scholes.cpp >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling black_scholes.cpp."
  exit
fi
rm compiler.txt
#
mv black_scholes.o ~/libcpp/$ARCH/black_scholes.o
#
echo "Library installed as ~/libcpp/$ARCH/black_scholes.o"
