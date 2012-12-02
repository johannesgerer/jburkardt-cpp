#!/bin/bash
#
cp keast.hpp /$HOME/include
#
g++ -c -g -I /$HOME/include keast.cpp >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling keast.cpp"
  exit
fi
rm compiler.txt
#
mv keast.o ~/libcpp/$ARCH/keast.o
#
echo "Library installed as ~/libcpp/$ARCH/keast.o"
