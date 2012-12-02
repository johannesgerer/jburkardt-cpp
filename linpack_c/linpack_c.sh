#!/bin/bash
#
cp linpack_c.hpp /$HOME/include
#
g++ -c -g -I /$HOME/include linpack_c.cpp >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling linpack_c.cpp"
  exit
fi
rm compiler.txt
#
mv linpack_c.o ~/libcpp/$ARCH/linpack_c.o
#
echo "Library installed as ~/libcpp/$ARCH/linpack_c.o"
