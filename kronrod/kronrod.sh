#!/bin/bash
#
cp kronrod.hpp /$HOME/include
#
g++ -c -g -I /$HOME/include kronrod.cpp >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling kronrod.cpp"
  exit
fi
rm compiler.txt
#
mv kronrod.o ~/libcpp/$ARCH/kronrod.o
#
echo "Library installed as ~/libcpp/$ARCH/kronrod.o"
