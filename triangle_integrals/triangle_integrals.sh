#!/bin/bash
#
cp triangle_integrals.hpp /$HOME/include
#
g++ -c -I/$HOME/include triangle_integrals.cpp
if [ $? -ne 0 ]; then
  echo "Errors compiling triangle_integrals.cpp"
  exit
fi
#
mv triangle_integrals.o ~/libcpp/$ARCH/triangle_integrals.o
#
echo "Library installed as ~/libcpp/$ARCH/triangle_integrals.o"
