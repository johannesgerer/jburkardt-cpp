#!/bin/bash
#
cp xy_io.hpp /$HOME/include
#
g++ -c -g -I /$HOME/include xy_io.cpp >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling xy_io.cpp"
  exit
fi
rm compiler.txt
#
mv xy_io.o ~/libcpp/$ARCH/xy_io.o
#
echo "Library installed as ~/libcpp/$ARCH/xy_io.o"
