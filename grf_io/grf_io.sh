#!/bin/bash
#
cp grf_io.hpp /$HOME/include
#
g++ -c -g -I /$HOME/include grf_io.cpp >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling grf_io.cpp"
  exit
fi
rm compiler.txt
#
mv grf_io.o ~/libcpp/$ARCH/grf_io.o
#
echo "Library installed as ~/libcpp/$ARCH/grf_io.o"
