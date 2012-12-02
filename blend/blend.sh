#!/bin/bash
#
cp blend.hpp /$HOME/include
#
g++ -c -g -I /$HOME/include blend.cpp >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling blend.cpp"
  exit
fi
rm compiler.txt
#
mv blend.o ~/libcpp/$ARCH/blend.o
#
echo "Library installed as ~/libcpp/$ARCH/blend.o"
