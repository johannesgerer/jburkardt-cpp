#!/bin/bash
#
cp ncc_triangle.hpp /$HOME/include
#
g++ -c -g -I /$HOME/include ncc_triangle.cpp >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling ncc_triangle.cpp"
  exit
fi
rm compiler.txt
#
mv ncc_triangle.o ~/libcpp/$ARCH/ncc_triangle.o
#
echo "Library installed as ~/libcpp/$ARCH/ncc_triangle.o"
