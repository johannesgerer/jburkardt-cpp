#!/bin/bash
#
cp asa111.hpp /$HOME/include
#
g++ -c -g -I /$HOME/include asa111.cpp >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling asa111.cpp"
  exit
fi
rm compiler.txt
#
mv asa111.o ~/libcpp/$ARCH/asa111.o
#
echo "Library installed as ~/libcpp/$ARCH/asa111.o"
