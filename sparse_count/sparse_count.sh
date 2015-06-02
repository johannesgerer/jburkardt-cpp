#!/bin/bash
#
cp sparse_count.hpp /$HOME/include
#
g++ -c -g -I /$HOME/include sparse_count.cpp >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling sparse_count.cpp"
  exit
fi
rm compiler.txt
#
mv sparse_count.o ~/libcpp/$ARCH/sparse_count.o
#
echo "Library installed as ~/libcpp/$ARCH/sparse_count.o"
