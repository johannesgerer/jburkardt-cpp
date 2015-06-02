#!/bin/bash
#
cp hyperball_integrals.hpp /$HOME/include
#
g++ -c -I/$HOME/include hyperball_integrals.cpp
if [ $? -ne 0 ]; then
  echo "Errors compiling hyperball_integrals.cpp"
  exit
fi
#
mv hyperball_integrals.o ~/libcpp/$ARCH/hyperball_integrals.o
#
echo "Library installed as ~/libcpp/$ARCH/hyperball_integrals.o"
