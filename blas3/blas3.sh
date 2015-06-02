#!/bin/bash
#
cp blas3.hpp /$HOME/include
#
g++ -c -O2 -I$HOME/include blas3.cpp >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling blas3.cpp."
  exit
fi
rm compiler.txt
#
mv blas3.o ~/libcpp/$ARCH/blas3.o
#
echo "Library installed as ~/libcpp/$ARCH/blas3.o"
