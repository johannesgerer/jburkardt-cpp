#!/bin/bash
#
cp latin_cover.hpp /$HOME/include
#
g++ -c -g -I /$HOME/include latin_cover.cpp >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling latin_cover.cpp"
  exit
fi
rm compiler.txt
#
mv latin_cover.o ~/libcpp/$ARCH/latin_cover.o
#
echo "Library installed as ~/libcpp/$ARCH/latin_cover.o"
