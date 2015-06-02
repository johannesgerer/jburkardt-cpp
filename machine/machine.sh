#!/bin/bash
#
cp machine.hpp /$HOME/include
#
g++ -c -I /$HOME/include machine.cpp
if [ $? -ne 0 ]; then
  echo "Errors compiling machine.cpp"
  exit
fi
#
mv machine.o ~/libcpp/$ARCH/machine.o
#
echo "Library installed as ~/libcpp/$ARCH/machine.o"
