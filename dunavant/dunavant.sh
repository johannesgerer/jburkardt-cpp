#!/bin/bash
#
cp dunavant.hpp /$HOME/include
#
g++ -c -g -I /$HOME/include dunavant.cpp >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling dunavant.cpp"
  exit
fi
rm compiler.txt
#
mv dunavant.o ~/libcpp/$ARCH/dunavant.o
#
echo "Library installed as ~/libcpp/$ARCH/dunavant.o"
