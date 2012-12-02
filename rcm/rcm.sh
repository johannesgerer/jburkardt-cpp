#!/bin/bash
#
cp rcm.hpp /$HOME/include
#
g++ -c -g -I /$HOME/include rcm.cpp >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling rcm.cpp"
  exit
fi
rm compiler.txt
#
mv rcm.o ~/libcpp/$ARCH/rcm.o
#
echo "Library installed as ~/libcpp/$ARCH/rcm.o"
