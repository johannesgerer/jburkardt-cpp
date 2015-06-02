#!/bin/bash
#
cp disk_grid.hpp /$HOME/include
#
g++ -c -g disk_grid.cpp >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling disk_grid.cpp"
  exit
fi
rm compiler.txt
#
mv disk_grid.o ~/libcpp/$ARCH/disk_grid.o
#
echo "Library installed as ~/libcpp/$ARCH/disk_grid.o"
