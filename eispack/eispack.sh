#!/bin/bash
#
cp eispack.hpp /$HOME/include
#
g++ -c -g -I/$HOME/include eispack.cpp >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling eispack.cpp"
  exit
fi
rm compiler.txt
#
mv eispack.o ~/libcpp/$ARCH/eispack.o
#
echo "Library installed as ~/libcpp/$ARCH/eispack.o"
