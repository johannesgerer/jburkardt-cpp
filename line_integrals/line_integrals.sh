#!/bin/bash
#
cp line_integrals.hpp /$HOME/include
#
g++ -c -g -I /$HOME/include line_integrals.cpp >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling line_integrals.cpp"
  exit
fi
rm compiler.txt
#
mv line_integrals.o ~/libcpp/$ARCH/line_integrals.o
#
echo "Library installed as ~/libcpp/$ARCH/line_integrals.o"
