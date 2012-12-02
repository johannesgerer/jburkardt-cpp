#!/bin/bash
#
cp collatz_recursive.hpp /$HOME/include
#
g++ -c -g collatz_recursive.cpp >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling collatz_recursive.cpp"
  exit
fi
rm compiler.txt
#
mv collatz_recursive.o ~/libcpp/$ARCH/collatz_recursive.o
#
echo "Library installed as ~/libcpp/$ARCH/collatz_recursive.o"
