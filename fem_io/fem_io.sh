#!/bin/bash
#
cp fem_io.hpp /$HOME/include
#
g++ -c -g -I /$HOME/include fem_io.cpp >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling fem_io.cpp"
  exit
fi
rm compiler.txt
#
mv fem_io.o ~/libcpp/$ARCH/fem_io.o
#
echo "Library installed as ~/libcpp/$ARCH/fem_io.o"
