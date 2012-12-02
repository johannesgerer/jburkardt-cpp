#!/bin/bash
#
cp i8lib.hpp /$HOME/include
#
g++ -c -g -I /$HOME/include i8lib.cpp >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling i8lib.cpp"
  exit
fi
rm compiler.txt
#
mv i8lib.o ~/libcpp/$ARCH/i8lib.o
#
echo "Library installed as ~/libcpp/$ARCH/i8lib.o"
