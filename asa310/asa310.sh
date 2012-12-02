#!/bin/bash
#
cp asa310.hpp /$HOME/include
#
g++ -c -g -I /$HOME/include asa310.cpp >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling asa310.cpp"
  exit
fi
rm compiler.txt
#
mv asa310.o ~/libcpp/$ARCH/asa310.o
#
echo "Library installed as ~/libcpp/$ARCH/asa310.o"
