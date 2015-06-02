#!/bin/bash
#
cp cpv.hpp /$HOME/include
#
g++ -c -I/$HOME/include cpv.cpp
if [ $? -ne 0 ]; then
  echo "Errors compiling cpv.cpp"
  exit
fi
#
mv cpv.o ~/libcpp/$ARCH/cpv.o
#
echo "Library installed as ~/libcpp/$ARCH/cpv.o"
