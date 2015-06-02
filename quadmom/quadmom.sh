#!/bin/bash
#
cp quadmom.hpp /$HOME/include
#
g++ -c -g -I/$HOME/include quadmom.cpp >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling quadmom.cpp"
  exit
fi
rm compiler.txt
#
mv quadmom.o ~/libcpp/$ARCH/quadmom.o
#
echo "Library installed as ~/libcpp/$ARCH/quadmom.o"
