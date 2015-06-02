#!/bin/bash
#
cp c4_complex_lib.hpp /$HOME/include
#
g++ -c -g c4_complex_lib.cpp >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling c4_complex_lib.cpp"
  exit
fi
rm compiler.txt
#
mv c4_complex_lib.o ~/libcpp/$ARCH/c4_complex_lib.o
#
echo "Library installed as ~/libcpp/$ARCH/c4_complex_lib.o"
