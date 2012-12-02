#!/bin/bash
#
cp test_int_laguerre.hpp /$HOME/include
#
g++ -c -g -I /$HOME/include test_int_laguerre.cpp >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling test_int_laguerre.cpp"
  exit
fi
rm compiler.txt
#
mv test_int_laguerre.o ~/libcpp/$ARCH/test_int_laguerre.o
#
echo "Library installed as ~/libcpp/$ARCH/test_int_laguerre.o"
