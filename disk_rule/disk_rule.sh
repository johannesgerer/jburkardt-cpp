#!/bin/bash
#
cp disk_rule.hpp /$HOME/include
#
g++ -c -g -I/$HOME/include disk_rule.cpp >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling disk_rule.cpp"
  exit
fi
rm compiler.txt
#
mv disk_rule.o ~/libcpp/$ARCH/disk_rule.o
#
echo "Library installed as ~/libcpp/$ARCH/disk_rule.o"
