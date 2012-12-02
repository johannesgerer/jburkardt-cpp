#!/bin/bash
#
cp asa144.hpp /$HOME/include
#
g++ -c -g -I /$HOME/include asa144.cpp >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling asa144.cpp"
  exit
fi
rm compiler.txt
#
mv asa144.o ~/libcpp/$ARCH/asa144.o
#
echo "Library installed as ~/libcpp/$ARCH/asa144.o"
