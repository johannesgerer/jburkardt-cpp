#!/bin/bash
#
cp test_matrix_exponential.hpp /$HOME/include
#
g++ -c -g -I /$HOME/include test_matrix_exponential.cpp >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling test_matrix_exponential.cpp"
  exit
fi
rm compiler.txt
#
mv test_matrix_exponential.o ~/libcpp/$ARCH/test_matrix_exponential.o
#
echo "Library installed as ~/libcpp/$ARCH/test_matrix_exponential.o"
