#!/bin/bash
#
cp asa136.hpp /$HOME/include
#
g++ -c -g -I /$HOME/include asa136.cpp >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling asa136.cpp"
  exit
fi
rm compiler.txt
#
mv asa136.o ~/libcpp/$ARCH/asa136.o
#
echo "Library installed as ~/libcpp/$ARCH/asa136.o"
