#!/bin/bash
#
cp test_values.hpp /$HOME/include
#
g++ -c -g -I /$HOME/include test_values.cpp >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling test_values.cpp"
  exit
fi
rm compiler.txt
#
mv test_values.o ~/libcpp/$ARCH/test_values.o
#
echo "Library installed as ~/libcpp/$ARCH/test_values.o"
