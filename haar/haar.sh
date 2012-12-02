#!/bin/bash
#
cp haar.hpp /$HOME/include
#
g++ -c -g -I /$HOME/include haar.cpp >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling haar.cpp"
  exit
fi
rm compiler.txt
#
mv haar.o ~/libcpp/$ARCH/haar.o
#
echo "Library installed as ~/libcpp/$ARCH/haar.o"
