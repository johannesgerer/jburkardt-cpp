#!/bin/bash
#
cp calpak.hpp /$HOME/include
#
g++ -c -g -I /$HOME/include calpak.cpp >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling calpak.cpp"
  exit
fi
rm compiler.txt
#
mv calpak.o ~/libcpp/$ARCH/calpak.o
#
echo "Library installed as ~/libcpp/$ARCH/calpak.o"
