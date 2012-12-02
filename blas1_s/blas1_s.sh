#!/bin/bash
#
cp blas1_s.hpp /$HOME/include
#
g++ -c -g -I /$HOME/include blas1_s.cpp >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling blas1_s.cpp"
  exit
fi
rm compiler.txt
#
mv blas1_s.o ~/libcpp/$ARCH/blas1_s.o
#
echo "Library installed as ~/libcpp/$ARCH/blas1_s.o"
