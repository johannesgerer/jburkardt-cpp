#!/bin/bash
#
cp hermite.hpp /$HOME/include
#
g++ -c -I /$HOME/include hermite.cpp
if [ $? -ne 0 ]; then
  echo "Errors compiling hermite.cpp"
  exit
fi
#
mv hermite.o ~/libcpp/$ARCH/hermite.o
#
echo "Library installed as ~/libcpp/$ARCH/hermite.o"
