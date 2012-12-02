#!/bin/bash
#
cp asa243.hpp /$HOME/include
#
g++ -c -g -I /$HOME/include asa243.cpp >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling asa243.cpp"
  exit
fi
rm compiler.txt
#
mv asa243.o ~/libcpp/$ARCH/asa243.o
#
echo "Library installed as ~/libcpp/$ARCH/asa243.o"
