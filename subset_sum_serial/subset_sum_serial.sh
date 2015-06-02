#!/bin/bash
#
cp subset_sum_serial.hpp /$HOME/include
#
g++ -c -g -I/$HOME/include subset_sum_serial.cpp >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling subset_sum_serial.cpp"
  exit
fi
rm compiler.txt
#
mv subset_sum_serial.o ~/libcpp/$ARCH/subset_sum_serial.o
#
echo "Library installed as ~/libcpp/$ARCH/subset_sum_serial.o"
