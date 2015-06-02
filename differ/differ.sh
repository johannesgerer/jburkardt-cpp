#!/bin/bash
#
cp differ.hpp /$HOME/include
#
g++ -c -g -I/$HOME/include differ.cpp >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling differ.cpp"
  exit
fi
rm compiler.txt
#
mv differ.o ~/libcpp/$ARCH/differ.o
#
echo "Library installed as ~/libcpp/$ARCH/differ.o"
