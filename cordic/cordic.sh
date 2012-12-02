#!/bin/bash
#
cp cordic.hpp /$HOME/include
#
g++ -c -g -I /$HOME/include cordic.cpp >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling cordic.cpp"
  exit
fi
rm compiler.txt
#
mv cordic.o ~/libcpp/$ARCH/cordic.o
#
echo "Library installed as ~/libcpp/$ARCH/cordic.o"
