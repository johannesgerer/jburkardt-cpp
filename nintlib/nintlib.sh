#!/bin/bash
#
cp nintlib.hpp /$HOME/include
#
g++ -c -g -I /$HOME/include nintlib.cpp >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling nintlib.cpp"
  exit
fi
rm compiler.txt
#
mv nintlib.o ~/libcpp/$ARCH/nintlib.o
#
echo "Library installed as ~/libcpp/$ARCH/nintlib.o"
