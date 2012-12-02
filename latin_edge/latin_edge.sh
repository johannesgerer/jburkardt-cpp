#!/bin/bash
#
cp latin_edge.hpp /$HOME/include
#
g++ -c -g -I /$HOME/include latin_edge.cpp >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling latin_edge.cpp"
  exit
fi
rm compiler.txt
#
mv latin_edge.o ~/libcpp/$ARCH/latin_edge.o
#
echo "Library installed as ~/libcpp/$ARCH/latin_edge.o"
