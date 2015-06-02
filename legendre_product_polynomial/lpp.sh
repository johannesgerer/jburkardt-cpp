#!/bin/bash
#
cp lpp.hpp /$HOME/include
#
g++ -c -I /$HOME/include lpp.cpp
if [ $? -ne 0 ]; then
  echo "Errors compiling lpp.cpp"
  exit
fi
#
mv lpp.o ~/libcpp/$ARCH/lpp.o
#
echo "Library installed as ~/libcpp/$ARCH/lpp.o"
