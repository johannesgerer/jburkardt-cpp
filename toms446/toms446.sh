#!/bin/bash
#
cp toms446.hpp /$HOME/include
#
g++ -c -g -I /$HOME/include toms446.cpp >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling toms446.cpp"
  exit
fi
rm compiler.txt
#
mv toms446.o ~/libcpp/$ARCH/toms446.o
#
echo "Library installed as ~/libcpp/$ARCH/toms446.o"
