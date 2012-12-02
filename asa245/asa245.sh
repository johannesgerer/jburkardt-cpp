#!/bin/bash
#
cp asa245.hpp /$HOME/include
#
g++ -c -g -I /$HOME/include asa245.cpp >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling asa245.cpp"
  exit
fi
rm compiler.txt
#
mv asa245.o ~/libcpp/$ARCH/asa245.o
#
echo "Library installed as ~/libcpp/$ARCH/asa245.o"
