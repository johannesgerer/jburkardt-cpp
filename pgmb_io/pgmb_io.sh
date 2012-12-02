#!/bin/bash
#
cp pgmb_io.hpp /$HOME/include
#
g++ -c -g -I /$HOME/include pgmb_io.cpp >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling pgmb_io.cpp"
  exit
fi
rm compiler.txt
#
mv pgmb_io.o ~/libcpp/$ARCH/pgmb_io.o
#
echo "Library installed as ~/libcpp/$ARCH/pgmb_io.o"
