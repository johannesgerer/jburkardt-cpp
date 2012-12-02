#!/bin/bash
#
cp blas1_d.hpp /$HOME/include
#
g++ -c -g -I /$HOME/include blas1_d.cpp >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling blas1_d.cpp"
  exit
fi
rm compiler.txt
#
mv blas1_d.o ~/libcpp/$ARCH/blas1_d.o
#
echo "Library installed as ~/libcpp/$ARCH/blas1_d.o"
