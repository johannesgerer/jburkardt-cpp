#!/bin/bash
#
cp vector_read.hpp /$HOME/include
#
g++ -c -g -I /$HOME/include vector_read.cpp >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling vector_read.cpp"
  exit
fi
rm compiler.txt
#
mv vector_read.o ~/libcpp/$ARCH/vector_read.o
#
echo "Library installed as ~/libcpp/$ARCH/vector_read.o"
