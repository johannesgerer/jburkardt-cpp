#!/bin/bash
#
cp toms097.hpp /$HOME/include
#
g++ -c -g -I/$HOME/include toms097.cpp >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling toms097.cpp"
  exit
fi
rm compiler.txt
#
mv toms097.o ~/libcpp/$ARCH/toms097.o
#
echo "Library installed as ~/libcpp/$ARCH/toms097.o"
