#!/bin/bash
#
cp subset_sum.hpp /$HOME/include
#
g++ -c -g -I /$HOME/include subset_sum.cpp >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling subset_sum.cpp."
  exit
fi
rm compiler.txt
#
mv subset_sum.o ~/libcpp/$ARCH/subset_sum.o
#
echo "Library installed as ~/libcpp/$ARCH/subset_sum.o"
