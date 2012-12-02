#!/bin/bash
#
cp floyd.hpp /$HOME/include
#
g++ -c -g -I /$HOME/include floyd.cpp >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling floyd.cpp"
  exit
fi
rm compiler.txt
#
mv floyd.o ~/libcpp/$ARCH/floyd.o
#
echo "Library installed as ~/libcpp/$ARCH/floyd.o"
