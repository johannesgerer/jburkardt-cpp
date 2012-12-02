#!/bin/bash
#
cp normal.hpp /$HOME/include
#
g++ -c -g -I /$HOME/include normal.cpp >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling normal.cpp"
  exit
fi
rm compiler.txt
#
mv normal.o ~/libcpp/$ARCH/normal.o
#
echo "Library installed as ~/libcpp/$ARCH/normal.o"
