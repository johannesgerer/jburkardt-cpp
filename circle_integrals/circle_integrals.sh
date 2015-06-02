#!/bin/bash
#
cp circle_integrals.hpp /$HOME/include
#
g++ -c -g -I /$HOME/include circle_integrals.cpp >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling circle_integrals.cpp"
  exit
fi
rm compiler.txt
#
mv circle_integrals.o ~/libcpp/$ARCH/circle_integrals.o
#
echo "Library installed as ~/libcpp/$ARCH/circle_integrals.o"
