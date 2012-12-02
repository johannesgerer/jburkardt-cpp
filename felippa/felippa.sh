#!/bin/bash
#
cp felippa.hpp /$HOME/include
#
g++ -c -g -I /$HOME/include felippa.cpp >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling felippa.cpp"
  exit
fi
rm compiler.txt
#
mv felippa.o ~/libcpp/$ARCH/felippa.o
#
echo "Library installed as ~/libcpp/$ARCH/felippa.o"
