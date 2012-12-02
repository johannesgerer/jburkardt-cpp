#!/bin/bash
#
cp ice_io.hpp /$HOME/include
#
g++ -c -g -I /usr/local/include ice_io.cpp >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling ice_io.cpp"
  exit
fi
rm compiler.txt
#
mv ice_io.o ~/libcpp/$ARCH/ice_io.o
#
echo "Library installed as ~/libcpp/$ARCH/ice_io.o"
