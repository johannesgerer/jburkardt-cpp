#!/bin/bash
#
cp tanh_quad.hpp /$HOME/include
#
g++ -c -g -I /$HOME/include tanh_quad.cpp >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling tanh_quad.cpp"
  exit
fi
rm compiler.txt
#
mv tanh_quad.o ~/libcpp/$ARCH/tanh_quad.o
#
echo "Library installed as ~/libcpp/$ARCH/tanh_quad.o"
