#!/bin/bash
#
cp asa226.hpp /$HOME/include
#
g++ -c -g -I /$HOME/include asa226.cpp >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling asa226.cpp"
  exit
fi
rm compiler.txt
#
mv asa226.o ~/libcpp/$ARCH/asa226.o
#
echo "Library installed as ~/libcpp/$ARCH/asa226.o"
