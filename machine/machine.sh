#!/bin/bash
#
cp machine.hpp /$HOME/include
#
g++ -c -g -I /$HOME/include machine.cpp >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling machine.cpp"
  exit
fi
rm compiler.txt
#
mv machine.o ~/libcpp/$ARCH/machine.o
#
echo "Library installed as ~/libcpp/$ARCH/machine.o"
