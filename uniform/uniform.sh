#!/bin/bash
#
cp uniform.hpp /$HOME/include
#
g++ -c -g -I /$HOME/include uniform.cpp >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling uniform.cpp"
  exit
fi
rm compiler.txt
#
mv uniform.o ~/libcpp/$ARCH/uniform.o
#
echo "Library installed as ~/libcpp/$ARCH/uniform.o"
