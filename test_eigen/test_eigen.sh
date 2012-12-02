#!/bin/bash
#
cp test_eigen.hpp /$HOME/include
#
g++ -c -g -I /$HOME/include test_eigen.cpp >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling test_eigen.cpp."
  exit
fi
rm compiler.txt
#
mv test_eigen.o ~/libcpp/$ARCH/test_eigen.o
#
echo "Library installed as ~/libcpp/$ARCH/test_eigen.o"
