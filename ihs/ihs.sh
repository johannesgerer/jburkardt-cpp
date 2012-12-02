#!/bin/bash
#
cp ihs.hpp /$HOME/include
#
g++ -c -g -I /$HOME/include ihs.cpp >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling ihs.cpp"
  exit
fi
rm compiler.txt
#
mv ihs.o ~/libcpp/$ARCH/ihs.o
#
echo "Library installed as ~/libcpp/$ARCH/ihs.o"
