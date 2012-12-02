#!/bin/bash
#
cp lattice_rule.hpp /$HOME/include
#
g++ -c -g -I /$HOME/include lattice_rule.cpp >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling lattice_rule.cpp"
  exit
fi
rm compiler.txt
#
mv lattice_rule.o ~/libcpp/$ARCH/lattice_rule.o
#
echo "Library installed as ~/libcpp/$ARCH/lattice_rule.o"
