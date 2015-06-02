#!/bin/bash
#
cp padua.hpp /$HOME/include
#
g++ -c -I/$HOME/include padua.cpp
if [ $? -ne 0 ]; then
  echo "Errors compiling padua.cpp"
  exit
fi
#
mv padua.o ~/libcpp/$ARCH/padua.o
#
echo "Library installed as ~/libcpp/$ARCH/padua.o"
