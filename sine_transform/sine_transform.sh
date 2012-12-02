#!/bin/bash
#
cp sine_transform.hpp /$HOME/include
#
g++ -c -g -I /$HOME/include sine_transform.cpp >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling sine_transform.cpp"
  exit
fi
rm compiler.txt
#
mv sine_transform.o ~/libcpp/$ARCH/sine_transform.o
#
echo "Library installed as ~/libcpp/$ARCH/sine_transform.o"
