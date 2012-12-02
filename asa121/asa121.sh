#!/bin/bash
#
cp asa121.hpp /$HOME/include
#
g++ -c -g -I /$HOME/include asa121.cpp >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling asa121.cpp"
  exit
fi
rm compiler.txt
#
mv asa121.o ~/libcpp/$ARCH/asa121.o
#
echo "Library installed as ~/libcpp/$ARCH/asa121.o"
