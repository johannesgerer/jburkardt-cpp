#!/bin/bash
#
cp test_int_2d.hpp /$HOME/include
#
g++ -c -g -I /$HOME/include test_int_2d.cpp >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling test_int_2d.cpp"
  exit
fi
rm compiler.txt
#
mv test_int_2d.o ~/libcpp/$ARCH/test_int_2d.o
#
echo "Library installed as ~/libcpp/$ARCH/test_int_2d.o"
