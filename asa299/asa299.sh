#!/bin/bash
#
cp asa299.hpp /$HOME/include
#
g++ -c -g -I /$HOME/include asa299.cpp >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling asa299.cpp"
  exit
fi
rm compiler.txt
#
mv asa299.o ~/libcpp/$ARCH/asa299.o
#
echo "Library installed as ~/libcpp/$ARCH/asa299.o"
