#!/bin/bash
#
cp van_der_corput.hpp /$HOME/include
#
g++ -c -I /$HOME/include van_der_corput.cpp
if [ $? -ne 0 ]; then
  echo "Errors compiling van_der_corput.cpp"
  exit
fi
#
mv van_der_corput.o ~/libcpp/$ARCH/van_der_corput.o
#
echo "Library installed as ~/libcpp/$ARCH/van_der_corput.o"
