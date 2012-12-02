#!/bin/bash
#
cp van_der_corput.hpp /$HOME/include
#
g++ -c -g -I /$HOME/include van_der_corput.cpp >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling van_der_corput.cpp"
  exit
fi
rm compiler.txt
#
mv van_der_corput.o ~/libcpp/$ARCH/van_der_corput.o
#
echo "Library installed as ~/libcpp/$ARCH/van_der_corput.o"
