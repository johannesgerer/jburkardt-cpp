#!/bin/bash
#
cp pgma_io.hpp /$HOME/include
#
g++ -c -g -I /$HOME/include pgma_io.cpp >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling pgma_io.cpp"
  exit
fi
rm compiler.txt
#
mv pgma_io.o ~/libcpp/$ARCH/pgma_io.o
#
echo "Library installed as ~/libcpp/$ARCH/pgma_io.o"
