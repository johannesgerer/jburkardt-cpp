#!/bin/bash
#
cp fem1d_lagrange.hpp /$HOME/include
#
g++ -c -I /$HOME/include fem1d_lagrange.cpp
if [ $? -ne 0 ]; then
  echo "Errors compiling fem1d_lagrange.cpp"
  exit
fi
#
mv fem1d_lagrange.o ~/libcpp/$ARCH/fem1d_lagrange.o
#
echo "Library installed as ~/libcpp/$ARCH/fem1d_lagrange.o"
