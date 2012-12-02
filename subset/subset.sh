#!/bin/bash
#
cp subset.hpp /$HOME/include
#
g++ -c -g -I /$HOME/include subset.cpp >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling subset.cpp"
  exit
fi
rm compiler.txt
#
mv subset.o ~/libcpp/$ARCH/subset.o
#
echo "Library installed as ~/libcpp/$ARCH/subset.o"
