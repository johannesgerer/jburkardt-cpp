#!/bin/bash
#
cp ziggurat_inline.hpp /$HOME/include
#
g++ -c -g -I /$HOME/include ziggurat_inline.cpp >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling ziggurat_inline.cpp."
  exit
fi
rm compiler.txt
#
mv ziggurat_inline.o ~/libcpp/$ARCH/ziggurat_inline.o
#
echo "Library installed as ~/libcpp/$ARCH/ziggurat_inline.o"
