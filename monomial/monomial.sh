#!/bin/bash
#
cp monomial.hpp /$HOME/include
#
g++ -c -I/$HOME/include monomial.cpp
if [ $? -ne 0 ]; then
  echo "Errors compiling monomial.cpp"
  exit
fi
#
mv monomial.o ~/libcpp/$ARCH/monomial.o
#
echo "Library installed as ~/libcpp/$ARCH/monomial.o"
