#!/bin/bash
#
cp latin_center.hpp /$HOME/include
#
g++ -c -g -I /$HOME/include latin_center.cpp >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling latin_center.cpp"
  exit
fi
rm compiler.txt
#
mv latin_center.o ~/libcpp/$ARCH/latin_center.o
#
echo "Library installed as ~/libcpp/$ARCH/latin_center.o"
