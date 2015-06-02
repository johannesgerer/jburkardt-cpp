#!/bin/bash
#
cp sparse_display.hpp /$HOME/include
#
g++ -c -I/$HOME/include sparse_display.cpp
if [ $? -ne 0 ]; then
  echo "Errors compiling sparse_display.cpp"
  exit
fi
#
mv sparse_display.o ~/libcpp/$ARCH/sparse_display.o
#
echo "Library installed as ~/libcpp/$ARCH/sparse_display.o"
