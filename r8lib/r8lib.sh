#!/bin/bash
#
cp r8lib.hpp /$HOME/include
#
g++ -c -g -I /$HOME/include r8lib.cpp >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling r8lib.cpp"
  exit
fi
rm compiler.txt
#
mv r8lib.o ~/libcpp/$ARCH/r8lib.o
#
echo "Library installed as ~/libcpp/$ARCH/r8lib.o"
