#!/bin/bash
#
cp sandia_rules2.hpp /$HOME/include
#
g++ -c -g -I /$HOME/include sandia_rules2.cpp >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling sandia_rules2.cpp"
  exit
fi
rm compiler.txt
#
mv sandia_rules2.o ~/libcpp/$ARCH/sandia_rules2.o
#
echo "Library installed as ~/libcpp/$ARCH/sandia_rules2.o"
