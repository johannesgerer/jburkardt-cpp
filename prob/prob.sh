#!/bin/bash
#
cp prob.hpp /$HOME/include
#
g++ -c -g -I /$HOME/include prob.cpp >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling prob.cpp"
  exit
fi
rm compiler.txt
#
mv prob.o ~/libcpp/$ARCH/prob.o
#
echo "Library installed as ~/libcpp/$ARCH/prob.o"
