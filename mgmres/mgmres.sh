#!/bin/bash
#
cp mgmres.hpp /$HOME/include
#
g++ -c -g -I /$HOME/include mgmres.cpp >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling mgmres.cpp"
  exit
fi
rm compiler.txt
#
mv mgmres.o ~/libcpp/$ARCH/mgmres.o
#
echo "Library installed as ~/libcpp/$ARCH/mgmres.o"
