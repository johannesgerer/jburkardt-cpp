#!/bin/bash
#
cp crc.hpp /$HOME/include
#
g++ -c -g -I /$HOME/include crc.cpp >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling crc.cpp"
  exit
fi
rm compiler.txt
#
mv crc.o ~/libcpp/$ARCH/crc.o
#
echo "Library installed as ~/libcpp/$ARCH/crc.o"
