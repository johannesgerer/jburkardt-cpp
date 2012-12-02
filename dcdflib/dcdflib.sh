#!/bin/bash
#
cp dcdflib.hpp /$HOME/include
#
g++ -c -g -I /$HOME/include dcdflib.cpp >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling dcdflib.cpp"
  exit
fi
rm compiler.txt
#
mv dcdflib.o ~/libcpp/$ARCH/dcdflib.o
#
echo "Library installed as ~/libcpp/$ARCH/dcdflib.o"
