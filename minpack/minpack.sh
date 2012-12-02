#!/bin/bash
#
cp minpack.hpp /$HOME/include
#
g++ -c -g -I /$HOME/include minpack.cpp >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling minpack.cpp"
  exit
fi
rm compiler.txt
#
mv minpack.o ~/libcpp/$ARCH/minpack.o
#
echo "Library installed as ~/libcpp/$ARCH/minpack.o"
