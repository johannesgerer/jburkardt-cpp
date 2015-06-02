#!/bin/bash
#
cp qwgw.hpp /$HOME/include
#
g++ -c -g -I/$HOME/include qwgw.cpp >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling qwgw.cpp"
  exit
fi
rm compiler.txt
#
mv qwgw.o ~/libcpp/$ARCH/qwgw.o
#
echo "Library installed as ~/libcpp/$ARCH/qwgw.o"
