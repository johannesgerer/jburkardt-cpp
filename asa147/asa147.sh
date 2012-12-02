#!/bin/bash
#
cp asa147.hpp /$HOME/include
#
g++ -c -g -I /$HOME/include asa147.cpp >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling asa147.cpp"
  exit
fi
rm compiler.txt
#
mv asa147.o ~/libcpp/$ARCH/asa147.o
#
echo "Library installed as ~/libcpp/$ARCH/asa147.o"
