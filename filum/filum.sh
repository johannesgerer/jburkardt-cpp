#!/bin/bash
#
cp filum.hpp /$HOME/include
#
g++ -c -g -I /$HOME/include filum.cpp >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling filum.cpp"
  exit
fi
rm compiler.txt
#
mv filum.o ~/libcpp/$ARCH/filum.o
#
echo "Library installed as ~/libcpp/$ARCH/filum.o"
