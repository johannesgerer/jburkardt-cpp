#!/bin/bash
#
cp blas1_c.hpp /$HOME/include
#
g++ -c -g -I /$HOME/include blas1_c.cpp >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling blas1_c.cpp"
  exit
fi
rm compiler.txt
#
mv blas1_c.o ~/libcpp/$ARCH/blas1_c.o
#
echo "Library installed as ~/libcpp/$ARCH/blas1_c.o"
