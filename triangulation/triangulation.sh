#!/bin/bash
#
cp triangulation.hpp /$HOME/include
#
g++ -c -g -I /$HOME/include triangulation.cpp >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling triangulation.cpp"
  exit
fi
rm compiler.txt
#
mv triangulation.o ~/libcpp/$ARCH/triangulation.o
#
echo "Library installed as ~/libcpp/$ARCH/triangulation.o"
