#!/bin/bash
#
cp nco_triangle.hpp /$HOME/include
#
g++ -c -g -I /$HOME/include nco_triangle.cpp >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling nco_triangle.cpp"
  exit
fi
rm compiler.txt
#
mv nco_triangle.o ~/libcpp/$ARCH/nco_triangle.o
#
echo "Library installed as ~/libcpp/$ARCH/nco_triangle.o"
