#!/bin/bash
#
cp toms179.hpp /$HOME/include
#
g++ -c -g -I /$HOME/include toms179.cpp >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling toms179.cpp"
  exit
fi
rm compiler.txt
#
mv toms179.o ~/libcpp/$ARCH/toms179.o
#
echo "Library installed as ~/libcpp/$ARCH/toms179.o"
