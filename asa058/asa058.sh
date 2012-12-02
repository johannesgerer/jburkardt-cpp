#!/bin/bash
#
cp asa058.hpp /$HOME/include
#
g++ -c -g -I /$HOME/include asa058.cpp >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling asa058.cpp"
  exit
fi
rm compiler.txt
#
mv asa058.o ~/libcpp/$ARCH/asa058.o
#
echo "Library installed as ~/libcpp/$ARCH/asa058.o"
