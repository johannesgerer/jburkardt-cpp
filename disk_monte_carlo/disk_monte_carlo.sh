#!/bin/bash
#
cp disk_monte_carlo.hpp /$HOME/include
#
g++ -c -g -I/$HOME/include disk_monte_carlo.cpp >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling disk_monte_carlo.cpp"
  exit
fi
rm compiler.txt
#
mv disk_monte_carlo.o ~/libcpp/$ARCH/disk_monte_carlo.o
#
echo "Library installed as ~/libcpp/$ARCH/disk_monte_carlo.o"
