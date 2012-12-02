#!/bin/bash
#
cp faure.hpp /$HOME/include
#
g++ -c -g -I /$HOME/include faure.cpp >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling faure.cpp"
  exit
fi
rm compiler.txt
#
mv faure.o ~/libcpp/$ARCH/faure.o
#
echo "Library installed as ~/libcpp/$ARCH/faure.o"
