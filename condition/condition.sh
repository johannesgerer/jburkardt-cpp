#!/bin/bash
#
cp condition.hpp /$HOME/include
#
g++ -c -g -I /$HOME/include condition.cpp >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling condition.cpp"
  exit
fi
rm compiler.txt
#
mv condition.o ~/libcpp/$ARCH/condition.o
#
echo "Library installed as ~/libcpp/$ARCH/condition.o"
