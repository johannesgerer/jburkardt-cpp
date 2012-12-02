#!/bin/bash
#
cp linpack_d.hpp /$HOME/include
#
g++ -c -g -I /$HOME/include linpack_d.cpp >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling linpack_d.cpp"
  exit
fi
rm compiler.txt
#
mv linpack_d.o ~/libcpp/$ARCH/linpack_d.o
#
echo "Library installed as ~/libcpp/$ARCH/linpack_d.o"
