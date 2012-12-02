#!/bin/bash
#
cp polygon_moments.hpp /$HOME/include
#
g++ -c -g -I/$HOME/include polygon_moments.cpp >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling polygon_moments.cpp"
  exit
fi
rm compiler.txt
#
mv polygon_moments.o ~/libcpp/$ARCH/polygon_moments.o
#
echo "Library installed as ~/libcpp/$ARCH/polygon_moments.o"
