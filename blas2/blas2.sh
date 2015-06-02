#!/bin/bash
#
cp blas2.hpp $HOME/include
#
g++ -c -O2 -I$HOME/include blas2.cpp >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling blas2.cpp"
  exit
fi
rm compiler.txt
#
mv blas2.o ~/libcpp/$ARCH/blas2.o
#
echo "Library installed as ~/libcpp/$ARCH/blas2.o"
