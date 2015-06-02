#!/bin/bash
#
cp snakes.hpp /$HOME/include
#
g++ -c -I /$HOME/include snakes.cpp
if [ $? -ne 0 ]; then
  echo "Errors compiling snakes.cpp"
  exit
fi
#
mv snakes.o ~/libcpp/$ARCH/snakes.o
#
echo "Library installed as ~/libcpp/$ARCH/snakes.o"
