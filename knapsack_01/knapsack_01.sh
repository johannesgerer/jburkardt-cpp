#!/bin/bash
#
cp knapsack_01.hpp /$HOME/include
#
g++ -c -I/$HOME/include knapsack_01.cpp
if [ $? -ne 0 ]; then
  echo "Errors compiling knapsack_01.cpp"
  exit
fi
#
mv knapsack_01.o ~/libcpp/$ARCH/knapsack_01.o
#
echo "Library installed as ~/libcpp/$ARCH/knapsack_01.o"
