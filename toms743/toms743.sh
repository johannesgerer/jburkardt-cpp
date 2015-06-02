#!/bin/bash
#
cp toms743.hpp /$HOME/include
#
g++ -c -I/$HOME/include toms743.cpp
if [ $? -ne 0 ]; then
  echo "Errors compiling toms743.cpp"
  exit
fi
#
mv toms743.o ~/libcpp/$ARCH/toms743.o
#
echo "Library installed as ~/libcpp/$ARCH/toms743.o"
