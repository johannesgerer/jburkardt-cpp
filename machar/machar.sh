#!/bin/bash
#
cp machar.hpp /$HOME/include
#
g++ -c -I /$HOME/include machar.cpp
if [ $? -ne 0 ]; then
  echo "Errors compiling machar.cpp"
  exit
fi
#
mv machar.o ~/libcpp/$ARCH/machar.o
#
echo "Library installed as ~/libcpp/$ARCH/machar.o"
