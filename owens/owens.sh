#!/bin/bash
#
cp owens.hpp /$HOME/include
#
g++ -c -g -I /$HOME/include owens.cpp >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling owens.cpp"
  exit
fi
rm compiler.txt
#
mv owens.o ~/libcpp/$ARCH/owens.o
#
echo "Library installed as ~/libcpp/$ARCH/owens.o"
