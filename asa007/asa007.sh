#!/bin/bash
#
cp asa007.hpp /$HOME/include
#
g++ -c -g -I /$HOME/include asa007.cpp >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling asa007.cpp"
  exit
fi
rm compiler.txt
#
mv asa007.o ~/libcpp/$ARCH/asa007.o
#
echo "Library installed as ~/libcpp/$ARCH/asa007.o"
