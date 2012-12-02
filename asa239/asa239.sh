#!/bin/bash
#
cp asa239.hpp /$HOME/include
#
g++ -c -g -I /$HOME/include asa239.cpp >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling asa239.cpp"
  exit
fi
rm compiler.txt
#
mv asa239.o ~/libcpp/$ARCH/asa239.o
#
echo "Library installed as ~/libcpp/$ARCH/asa239.o"
