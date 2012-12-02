#!/bin/bash
#
cp point_merge.hpp /$HOME/include
#
g++ -c -g -I /$HOME/include point_merge.cpp >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling point_merge.cpp"
  exit
fi
rm compiler.txt
#
mv point_merge.o ~/libcpp/$ARCH/point_merge.o
#
echo "Library installed as ~/libcpp/$ARCH/point_merge.o"
