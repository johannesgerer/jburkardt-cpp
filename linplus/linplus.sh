#!/bin/bash
#
cp linplus.hpp /$HOME/include
#
g++ -c -g -I /$HOME/include linplus.cpp >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling linplus.cpp"
  exit
fi
rm compiler.txt
#
mv linplus.o ~/libcpp/$ARCH/linplus.o
#
echo "Library installed as ~/libcpp/$ARCH/linplus.o"
