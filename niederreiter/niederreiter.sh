#!/bin/bash
#
cp niederreiter.hpp /$HOME/include
#
g++ -c -g -I /$HOME/include niederreiter.cpp >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling niederreiter.cpp"
  exit
fi
rm compiler.txt
#
mv niederreiter.o ~/libcpp/$ARCH/niederreiter.o
#
echo "Library installed as ~/libcpp/$ARCH/niederreiter.o"
