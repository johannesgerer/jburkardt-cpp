#!/bin/bash
#
cp llsq.hpp /$HOME/include
#
g++ -c -g -I /$HOME/include llsq.cpp >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling llsq.cpp"
  exit
fi
rm compiler.txt
#
mv llsq.o ~/libcpp/$ARCH/llsq.o
#
echo "Library installed as ~/libcpp/$ARCH/llsq.o"
