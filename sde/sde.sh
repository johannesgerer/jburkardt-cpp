#!/bin/bash
#
cp sde.hpp /$HOME/include
#
g++ -c -g -I/$HOME/include sde.cpp >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling sde.cpp"
  exit
fi
rm compiler.txt
#
mv sde.o ~/libcpp/$ARCH/sde.o
#
echo "Library installed as ~/libcpp/$ARCH/sde.o"
