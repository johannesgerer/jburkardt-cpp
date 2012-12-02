#!/bin/bash
#
cp toms178.hpp /$HOME/include
#
g++ -c -g -I /$HOME/include toms178.cpp >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling toms178.cpp"
  exit
fi
rm compiler.txt
#
mv toms178.o ~/libcpp/$ARCH/toms178.o
#
echo "Library installed as ~/libcpp/$ARCH/toms178.o"
