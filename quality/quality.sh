#!/bin/bash
#
cp quality.hpp /$HOME/include
#
g++ -c -g -I /$HOME/include quality.cpp >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling quality.cpp"
  exit
fi
rm compiler.txt
#
mv quality.o ~/libcpp/$ARCH/quality.o
#
echo "Library installed as ~/libcpp/$ARCH/quality.o"
