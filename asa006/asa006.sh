#!/bin/bash
#
cp asa006.hpp /$HOME/include
#
g++ -c -g -I /$HOME/include asa006.cpp >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling asa006.cpp"
  exit
fi
rm compiler.txt
#
mv asa006.o ~/libcpp/$ARCH/asa006.o
#
echo "Library installed as ~/libcpp/$ARCH/asa006.o"
