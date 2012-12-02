#!/bin/bash
#
cp cvt.hpp /$HOME/include
#
g++ -c -g -I /$HOME/include cvt.cpp >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling cvt.cpp"
  exit
fi
rm compiler.txt
#
mv cvt.o ~/libcpp/$ARCH/cvt.o
#
echo "Library installed as ~/libcpp/$ARCH/cvt.o"
