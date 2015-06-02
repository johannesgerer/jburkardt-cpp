#!/bin/bash
#
cp qwv.hpp /$HOME/include
#
g++ -c -g -I/$HOME/include qwv.cpp >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling qwv.cpp"
  exit
fi
rm compiler.txt
#
mv qwv.o ~/libcpp/$ARCH/qwv.o
#
echo "Library installed as ~/libcpp/$ARCH/qwv.o"
