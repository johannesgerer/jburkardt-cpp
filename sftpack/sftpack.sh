#!/bin/bash
#
cp sftpack.hpp /$HOME/include
#
g++ -c -g -I /$HOME/include sftpack.cpp >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling sftpack.cpp"
  exit
fi
rm compiler.txt
#
mv sftpack.o ~/libcpp/$ARCH/sftpack.o
#
echo "Library installed as ~/libcpp/$ARCH/sftpack.o"
