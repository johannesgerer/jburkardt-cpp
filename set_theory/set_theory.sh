#!/bin/bash
#
cp set_theory.hpp /$HOME/include
#
g++ -c -g -I /$HOME/include set_theory.cpp >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling set_theory.cpp"
  exit
fi
rm compiler.txt
#
mv set_theory.o ~/libcpp/$ARCH/set_theory.o
#
echo "Library installed as ~/libcpp/$ARCH/set_theory.o"
