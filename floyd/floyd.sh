#!/bin/bash
#
cp floyd.hpp /$HOME/include
#
g++ -c -I /$HOME/include floyd.cpp
if [ $? -ne 0 ]; then
  echo "Errors compiling floyd.cpp"
  exit
fi
#
mv floyd.o ~/libcpp/$ARCH/floyd.o
#
echo "Library installed as ~/libcpp/$ARCH/floyd.o"
