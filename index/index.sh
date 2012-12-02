#!/bin/bash
#
cp index.hpp /$HOME/include
#
g++ -c -g -I/$HOME/include index.cpp >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling index.cpp"
  exit
fi
rm compiler.txt
#
mv index.o ~/libcpp/$ARCH/index.o
#
echo "Library installed as ~/libcpp/$ARCH/index.o"
