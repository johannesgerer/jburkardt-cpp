#!/bin/bash
#
cp test_interp.hpp /$HOME/include
#
g++ -c -g -I /$HOME/include test_interp.cpp >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling test_interp.cpp"
  exit
fi
rm compiler.txt
#
mv test_interp.o ~/libcpp/$ARCH/test_interp.o
#
echo "Library installed as ~/libcpp/$ARCH/test_interp.o"
