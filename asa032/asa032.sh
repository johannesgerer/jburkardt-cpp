#!/bin/bash
#
cp asa032.hpp /$HOME/include
#
g++ -c -g -I /$HOME/include asa032.cpp >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling asa032.cpp"
  exit
fi
rm compiler.txt
#
mv asa032.o ~/libcpp/$ARCH/asa032.o
#
echo "Library installed as ~/libcpp/$ARCH/asa032.o"
