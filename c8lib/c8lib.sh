#!/bin/bash
#
cp c8lib.hpp /$HOME/include
#
g++ -c -I /$HOME/include c8lib.cpp
if [ $? -ne 0 ]; then
  echo "Errors compiling c8lib.cpp"
  exit
fi
#
mv c8lib.o ~/libcpp/$ARCH/c8lib.o
#
echo "Library installed as ~/libcpp/$ARCH/c8lib.o"
