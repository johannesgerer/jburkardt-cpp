#!/bin/bash
#
cp bernstein.hpp /$HOME/include
#
g++ -c -g -I /$HOME/include bernstein.cpp >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling bernstein.cpp"
  exit
fi
rm compiler.txt
#
mv bernstein.o ~/libcpp/$ARCH/bernstein.o
#
echo "Library installed as ~/libcpp/$ARCH/bernstein.o"
