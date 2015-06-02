#!/bin/bash
#
cp disk_integrals.hpp /$HOME/include
#
g++ -c -g -I/$HOME/include disk_integrals.cpp >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling disk_integrals.cpp"
  exit
fi
rm compiler.txt
#
mv disk_integrals.o ~/libcpp/$ARCH/disk_integrals.o
#
echo "Library installed as ~/libcpp/$ARCH/disk_integrals.o"
