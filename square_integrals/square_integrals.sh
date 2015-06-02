#!/bin/bash
#
cp square_integrals.hpp /$HOME/include
#
g++ -c -g -I /$HOME/include square_integrals.cpp >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling square_integrals.cpp"
  exit
fi
rm compiler.txt
#
mv square_integrals.o ~/libcpp/$ARCH/square_integrals.o
#
echo "Library installed as ~/libcpp/$ARCH/square_integrals.o"
