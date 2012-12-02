#!/bin/bash
#
cp asa241.hpp /$HOME/include
#
g++ -c -g -I /$HOME/include asa241.cpp >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling asa241.cpp"
  exit
fi
rm compiler.txt
#
mv asa241.o ~/libcpp/$ARCH/asa241.o
#
echo "Library installed as ~/libcpp/$ARCH/asa241.o"
