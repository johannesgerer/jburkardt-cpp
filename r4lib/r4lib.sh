#!/bin/bash
#
cp r4lib.hpp /$HOME/include
#
g++ -c -g -I /$HOME/include r4lib.cpp >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling r4lib.cpp"
  exit
fi
rm compiler.txt
#
mv r4lib.o ~/libcpp/$ARCH/r4lib.o
#
echo "Library installed as ~/libcpp/$ARCH/r4lib.o"
