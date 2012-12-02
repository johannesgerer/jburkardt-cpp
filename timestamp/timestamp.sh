#!/bin/bash
#
cp timestamp.hpp /$HOME/include
#
g++ -c -g -I /$HOME/include timestamp.cpp >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling timestamp.cpp"
  exit
fi
rm compiler.txt
#
mv timestamp.o ~/libcpp/$ARCH/timestamp.o
#
echo "Library installed as ~/libcpp/$ARCH/timestamp.o"
