#!/bin/bash
#
cp jacobi.hpp /$HOME/include
#
g++ -c -g -I /$HOME/include jacobi.cpp >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling jacobi.cpp"
  exit
fi
rm compiler.txt
#
mv jacobi.o ~/libcpp/$ARCH/jacobi.o
#
echo "Library installed as ~/libcpp/$ARCH/jacobi.o"
