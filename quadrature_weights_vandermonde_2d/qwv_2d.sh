#!/bin/bash
#
cp qwv_2d.hpp /$HOME/include
#
g++ -c -I/$HOME/include qwv_2d.cpp >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling qwv_2d.cpp"
  exit
fi
rm compiler.txt
#
mv qwv_2d.o ~/libcpp/$ARCH/qwv_2d.o
#
echo "Library installed as ~/libcpp/$ARCH/qwv_2d.o"
