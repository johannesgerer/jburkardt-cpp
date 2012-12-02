#!/bin/bash
#
cp linpack_z.hpp /$HOME/include
#
g++ -c -g -I /$HOME/include linpack_z.cpp >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling linpack_z.cpp"
  exit
fi
rm compiler.txt
#
mv linpack_z.o ~/libcpp/$ARCH/linpack_z.o
#
echo "Library installed as ~/libcpp/$ARCH/linpack_z.o"
