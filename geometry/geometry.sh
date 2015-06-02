#!/bin/bash
#
cp geometry.hpp /$HOME/include
#
g++ -c -I /$HOME/include geometry.cpp
if [ $? -ne 0 ]; then
  echo "Errors compiling geometry.cpp"
  exit
fi
#
mv geometry.o ~/libcpp/$ARCH/geometry.o
#
echo "Library installed as ~/libcpp/$ARCH/geometry.o"
