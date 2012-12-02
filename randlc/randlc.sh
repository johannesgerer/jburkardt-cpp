#!/bin/bash
#
cp randlc.hpp /$HOME/include
#
g++ -c -g -I /$HOME/include randlc.cpp >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling randlc.cpp"
  exit
fi
rm compiler.txt
#
mv randlc.o ~/libcpp/$ARCH/randlc.o
#
echo "Library installed as ~/libcpp/$ARCH/randlc.o"
