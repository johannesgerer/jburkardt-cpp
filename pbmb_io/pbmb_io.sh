#!/bin/bash
#
cp pbmb_io.hpp /$HOME/include
#
g++ -c -g -I /$HOME/include pbmb_io.cpp >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling pbmb_io.cpp"
  exit
fi
rm compiler.txt
#
mv pbmb_io.o ~/libcpp/$ARCH/pbmb_io.o
#
echo "Library installed as ~/libcpp/$ARCH/pbmb_io.o"
