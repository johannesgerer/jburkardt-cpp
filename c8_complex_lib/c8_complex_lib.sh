#!/bin/bash
#
cp c8_complex_lib.hpp /$HOME/include
#
g++ -c -g c8_complex_lib.cpp >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling c8_complex_lib.cpp"
  exit
fi
rm compiler.txt
#
mv c8_complex_lib.o ~/libcpp/$ARCH/c8_complex_lib.o
#
echo "Library installed as ~/libcpp/$ARCH/c8_complex_lib.o"
