#!/bin/bash
#
cp solve.hpp /$HOME/include
#
g++ -c -I/$HOME/include solve.cpp
if [ $? -ne 0 ]; then
  echo "Errors compiling solve.cpp"
  exit
fi
#
mv solve.o ~/libcpp/$ARCH/solve.o
#
echo "Library installed as ~/libcpp/$ARCH/solve.o"
