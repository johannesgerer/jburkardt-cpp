#!/bin/bash
#
cp asa091.hpp /$HOME/include
#
g++ -c -g -I /$HOME/include asa091.cpp >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling asa091.cpp"
  exit
fi
rm compiler.txt
#
mv asa091.o ~/libcpp/$ARCH/asa091.o
#
echo "Library installed as ~/libcpp/$ARCH/asa091.o"
