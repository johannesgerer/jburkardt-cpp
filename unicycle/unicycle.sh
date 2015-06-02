#!/bin/bash
#
cp unicycle.hpp /$HOME/include
#
g++ -c -I /$HOME/include unicycle.cpp
if [ $? -ne 0 ]; then
  echo "Errors compiling unicycle.cpp"
  exit
fi
#
mv unicycle.o ~/libcpp/$ARCH/unicycle.o
#
echo "Library installed as ~/libcpp/$ARCH/unicycle.o"
