#!/bin/bash
#
cp asa152.hpp /$HOME/include
#
g++ -c -g -I /$HOME/include asa152.cpp >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling asa152.cpp"
  exit
fi
rm compiler.txt
#
mv asa152.o ~/libcpp/$ARCH/asa152.o
#
echo "Library installed as ~/libcpp/$ARCH/asa152.o"
