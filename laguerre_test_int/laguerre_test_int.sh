#!/bin/bash
#
cp laguerre_test_int.hpp /$HOME/include
#
g++ -c -g -I /$HOME/include laguerre_test_int.cpp >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling laguerre_test_int.cpp"
  exit
fi
rm compiler.txt
#
mv laguerre_test_int.o ~/libcpp/$ARCH/laguerre_test_int.o
#
echo "Library installed as ~/libcpp/$ARCH/laguerre_test_int.o"
