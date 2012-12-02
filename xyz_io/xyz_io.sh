#!/bin/bash
#
cp xyz_io.hpp /$HOME/include
#
g++ -c -g -I /$HOME/include xyz_io.cpp >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling xyz_io.cpp"
  exit
fi
rm compiler.txt
#
mv xyz_io.o ~/libcpp/$ARCH/xyz_io.o
#
echo "Library installed as ~/libcpp/$ARCH/xyz_io.o"
