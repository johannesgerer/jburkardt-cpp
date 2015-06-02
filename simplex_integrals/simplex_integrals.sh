#!/bin/bash
#
cp simplex_integrals.hpp /$HOME/include
#
g++ -c -g -I /$HOME/include simplex_integrals.cpp >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling simplex_integrals.cpp"
  exit
fi
rm compiler.txt
#
mv simplex_integrals.o ~/libcpp/$ARCH/simplex_integrals.o
#
echo "Library installed as ~/libcpp/$ARCH/simplex_integrals.o"
