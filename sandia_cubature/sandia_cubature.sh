#!/bin/bash
#
cp sandia_cubature.hpp /$HOME/include
#
g++ -c -g -I /$HOME/include sandia_cubature.cpp >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling sandia_cubature.cpp"
  exit
fi
rm compiler.txt
#
mv sandia_cubature.o ~/libcpp/$ARCH/sandia_cubature.o
#
echo "Library installed as ~/libcpp/$ARCH/sandia_cubature.o"
