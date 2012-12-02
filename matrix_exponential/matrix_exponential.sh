#!/bin/bash
#
cp matrix_exponential.hpp /$HOME/include
#
g++ -c -g -I/$HOME/include matrix_exponential.cpp >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling matrix_exponential.cpp"
  exit
fi
rm compiler.txt
#
mv matrix_exponential.o ~/libcpp/$ARCH/matrix_exponential.o
#
echo "Library installed as ~/libcpp/$ARCH/matrix_exponential.o"
