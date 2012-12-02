#!/bin/bash
#
cp stroud.hpp /$HOME/include
#
g++ -c -g -I /$HOME/include stroud.cpp >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling stroud.cpp"
  exit
fi
rm compiler.txt
#
mv stroud.o ~/libcpp/$ARCH/stroud.o
#
echo "Library installed as ~/libcpp/$ARCH/stroud.o"
