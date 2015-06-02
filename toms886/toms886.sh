#!/bin/bash
#
cp toms886.hpp /$HOME/include
#
g++ -c -g -I/$HOME/include toms886.cpp >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling toms886.cpp"
  exit
fi
rm compiler.txt
#
mv toms886.o ~/libcpp/$ARCH/toms886.o
#
echo "Library installed as ~/libcpp/$ARCH/toms886.o"
