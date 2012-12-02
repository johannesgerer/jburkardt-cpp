#!/bin/bash
#
cp test_zero.hpp /$HOME/include
#
g++ -c -g -I /$HOME/include test_zero.cpp >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling test_zero.cpp"
  exit
fi
rm compiler.txt
#
mv test_zero.o ~/libcpp/$ARCH/test_zero.o
#
echo "Library installed as ~/libcpp/$ARCH/test_zero.o"
