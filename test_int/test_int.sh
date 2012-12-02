#!/bin/bash
#
cp test_int.hpp /$HOME/include
#
g++ -c -g -I /$HOME/include test_int.cpp >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling test_int.cpp."
  exit
fi
rm compiler.txt
#
mv test_int.o ~/libcpp/$ARCH/test_int.o
#
echo "Library installed as ~/libcpp/$ARCH/test_int.o"
