#!/bin/bash
#
cp geometry.hpp /$HOME/include
#
g++ -c -g -I /$HOME/include geometry.cpp >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling geometry.cpp"
  exit
fi
rm compiler.txt
#
mv geometry.o ~/libcpp/$ARCH/geometry.o
#
echo "Library installed as ~/libcpp/$ARCH/geometry.o"
