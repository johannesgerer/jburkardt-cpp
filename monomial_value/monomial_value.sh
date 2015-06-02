#!/bin/bash
#
cp monomial_value.hpp /$HOME/include
#
g++ -c -I/$HOME/include monomial_value.cpp
if [ $? -ne 0 ]; then
  echo "Errors compiling monomial_value.cpp"
  exit
fi
#
mv monomial_value.o ~/libcpp/$ARCH/monomial_value.o
#
echo "Library installed as ~/libcpp/$ARCH/monomial_value.o"
