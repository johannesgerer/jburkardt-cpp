#!/bin/bash
#
cp test_int_hermite.hpp /$HOME/include
#
g++ -c -g -I /$HOME/include test_int_hermite.cpp >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling test_int_hermite.cpp"
  exit
fi
rm compiler.txt
#
mv test_int_hermite.o ~/libcpp/$ARCH/test_int_hermite.o
#
echo "Library installed as ~/libcpp/$ARCH/test_int_hermite.o"
