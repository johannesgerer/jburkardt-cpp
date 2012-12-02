#!/bin/bash
#
cp sor.hpp /$HOME/include
#
g++ -c -g -I /$HOME/include sor.cpp >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling sor.cpp"
  exit
fi
rm compiler.txt
#
mv sor.o ~/libcpp/$ARCH/sor.o
#
echo "Library installed as ~/libcpp/$ARCH/sor.o"
