#!/bin/bash
#
cp divdif.hpp /$HOME/include
#
g++ -c -g -I /$HOME/include divdif.cpp >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling divdif.cpp"
  exit
fi
rm compiler.txt
#
mv divdif.o ~/libcpp/$ARCH/divdif.o
#
echo "Library installed as ~/libcpp/$ARCH/divdif.o"
