#!/bin/bash
#
cp asa066.hpp /$HOME/include
#
g++ -c -g -I /$HOME/include asa066.cpp >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling asa066.cpp"
  exit
fi
rm compiler.txt
#
mv asa066.o ~/libcpp/$ARCH/asa066.o
#
echo "Library installed as ~/libcpp/$ARCH/asa066.o"
