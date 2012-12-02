#!/bin/bash
#
cp pbma_io.hpp /$HOME/include
#
g++ -c -g -I /$HOME/include pbma_io.cpp >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling pbma_io.cpp"
  exit
fi
rm compiler.txt
#
mv pbma_io.o ~/libcpp/$ARCH/pbma_io.o
#
echo "Library installed as ~/libcpp/$ARCH/pbma_io.o"
