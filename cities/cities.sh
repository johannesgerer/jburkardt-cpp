#!/bin/bash
#
cp cities.hpp /$HOME/include
#
g++ -c -g -I /$HOME/include cities.cpp >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling cities.cpp"
  exit
fi
rm compiler.txt
#
mv cities.o ~/libcpp/$ARCH/cities.o
#
echo "Library installed as ~/libcpp/$ARCH/cities.o"
