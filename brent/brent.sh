#!/bin/bash
#
cp brent.hpp /$HOME/include
#
g++ -c -g -I /$HOME/include brent.cpp >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling brent.cpp"
  exit
fi
rm compiler.txt
#
mv brent.o ~/libcpp/$ARCH/brent.o
#
echo "Library installed as ~/libcpp/$ARCH/brent.o"
