#!/bin/bash
#
cp barycentric_interp_1d.hpp /$HOME/include
#
g++ -c -g -I/$HOME/include barycentric_interp_1d.cpp >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling barycentric_interp_1d.cpp"
  exit
fi
rm compiler.txt
#
mv barycentric_interp_1d.o ~/libcpp/$ARCH/barycentric_interp_1d.o
#
echo "Library installed as ~/libcpp/$ARCH/barycentric_interp_1d.o"
