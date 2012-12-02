#!/bin/bash
#
cp asa159.hpp /$HOME/include
#
g++ -c -g -I /$HOME/include asa159.cpp >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling asa159.cpp"
  exit
fi
rm compiler.txt
#
mv asa159.o ~/libcpp/$ARCH/asa159.o
#
echo "Library installed as ~/libcpp/$ARCH/asa159.o"
