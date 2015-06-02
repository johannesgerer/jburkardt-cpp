#!/bin/bash
#
cp bellman_ford.hpp /$HOME/include
#
g++ -c -I /$HOME/include bellman_ford.cpp
if [ $? -ne 0 ]; then
  echo "Errors compiling bellman_ford.cpp"
  exit
fi
#
mv bellman_ford.o ~/libcpp/$ARCH/bellman_ford.o
#
echo "Library installed as ~/libcpp/$ARCH/bellman_ford.o"
