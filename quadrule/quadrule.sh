#!/bin/bash
#
cp quadrule.hpp /$HOME/include
#
g++ -c -g -I /$HOME/include quadrule.cpp >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling quadrule.cpp"
  exit
fi
rm compiler.txt
#
mv quadrule.o ~/libcpp/$ARCH/quadrule.o
#
echo "Library installed as ~/libcpp/$ARCH/quadrule.o"
