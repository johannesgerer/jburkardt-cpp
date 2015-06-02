#!/bin/bash
#
cp wishart.hpp /$HOME/include
#
g++ -c -I/$HOME/include wishart.cpp
if [ $? -ne 0 ]; then
  echo "Errors compiling wishart.cpp"
  exit
fi
#
mv wishart.o ~/libcpp/$ARCH/wishart.o
#
echo "Library installed as ~/libcpp/$ARCH/wishart.o"
