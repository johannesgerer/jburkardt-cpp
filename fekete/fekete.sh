#!/bin/bash
#
cp fekete.hpp /$HOME/include
#
g++ -c -g -I /$HOME/include fekete.cpp >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling fekete.cpp"
  exit
fi
rm compiler.txt
#
mv fekete.o ~/libcpp/$ARCH/fekete.o
#
echo "Library installed as ~/libcpp/$ARCH/fekete.o"
