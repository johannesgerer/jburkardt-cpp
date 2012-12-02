#!/bin/bash
#
cp latinize.hpp /$HOME/include
#
g++ -c -g -I /$HOME/include latinize.cpp >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling latinize.cpp"
  exit
fi
rm compiler.txt
#
mv latinize.o ~/libcpp/$ARCH/latinize.o
#
echo "Library installed as ~/libcpp/$ARCH/latinize.o"
