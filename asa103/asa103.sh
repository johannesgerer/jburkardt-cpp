#!/bin/bash
#
cp asa103.hpp /$HOME/include
#
g++ -c -g -I /$HOME/include asa103.cpp >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling asa103.cpp"
  exit
fi
rm compiler.txt
#
mv asa103.o ~/libcpp/$ARCH/asa103.o
#
echo "Library installed as ~/libcpp/$ARCH/asa103.o"
