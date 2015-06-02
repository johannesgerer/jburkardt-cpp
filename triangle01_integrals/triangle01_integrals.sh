#!/bin/bash
#
cp triangle01_integrals.hpp /$HOME/include
#
g++ -c -I /$HOME/include triangle01_integrals.cpp
if [ $? -ne 0 ]; then
  echo "Errors compiling triangle01_integrals.cpp"
  exit
fi
#
mv triangle01_integrals.o ~/libcpp/$ARCH/triangle01_integrals.o
#
echo "Library installed as ~/libcpp/$ARCH/triangle01_integrals.o"
