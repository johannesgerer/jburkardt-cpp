#!/bin/bash
#
cp interp.hpp /$HOME/include
#
g++ -c -g -I/$HOME/include interp.cpp >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling interp.cpp"
  exit
fi
rm compiler.txt
#
mv interp.o ~/libcpp/$ARCH/interp.o
#
echo "Library installed as ~/libcpp/$ARCH/interp.o"
