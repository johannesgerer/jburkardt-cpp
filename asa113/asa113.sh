#!/bin/bash
#
cp asa113.hpp /$HOME/include
#
g++ -c -g -I /$HOME/include asa113.cpp >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling asa113.cpp"
  exit
fi
rm compiler.txt
#
mv asa113.o ~/libcpp/$ARCH/asa113.o
#
echo "Library installed as ~/libcpp/$ARCH/asa113.o"
