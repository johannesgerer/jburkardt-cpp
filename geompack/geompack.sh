#!/bin/bash
#
cp geompack.hpp /$HOME/include
#
g++ -c -g -I /$HOME/include geompack.cpp >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling geompack.cpp"
  exit
fi
rm compiler.txt
#
mv geompack.o ~/libcpp/$ARCH/geompack.o
#
echo "Library installed as ~/libcpp/$ARCH/geompack.o"
