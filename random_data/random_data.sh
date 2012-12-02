#!/bin/bash
#
cp random_data.hpp /$HOME/include
#
g++ -c -g random_data.cpp >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling random_data.cpp"
  exit
fi
rm compiler.txt
#
mv random_data.o ~/libcpp/$ARCH/random_data.o
#
echo "Library installed as ~/libcpp/$ARCH/random_data.o"
