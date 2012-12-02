#!/bin/bash
#
cp walsh.hpp /$HOME/include
#
g++ -c -g -I /$HOME/include walsh.cpp >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling walsh.cpp"
  exit
fi
rm compiler.txt
#
mv walsh.o ~/libcpp/$ARCH/walsh.o
#
echo "Library installed as ~/libcpp/$ARCH/walsh.o"
