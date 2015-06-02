#!/bin/bash
#
cp cdflib.hpp /$HOME/include
#
g++ -c -g -I /$HOME/include cdflib.cpp >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling cdflib.cpp"
  exit
fi
rm compiler.txt
#
mv cdflib.o ~/libcpp/$ARCH/cdflib.o
#
echo "Library installed as ~/libcpp/$ARCH/cdflib.o"
