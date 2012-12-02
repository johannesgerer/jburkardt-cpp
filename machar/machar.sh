#!/bin/bash
#
cp machar.hpp /$HOME/include
#
g++ -c -g -I /$HOME/include machar.cpp >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling machar.cpp"
  exit
fi
rm compiler.txt
#
mv machar.o ~/libcpp/$ARCH/machar.o
#
echo "Library installed as ~/libcpp/$ARCH/machar.o"
